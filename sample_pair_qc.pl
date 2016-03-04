#!/usr/bin/env perl

# GOAL: Using allele counts across known SNP positions, generate metrics for sample swaps or
# contamination. Requires a list of TN-paired sample IDs and a folder with allele counts
#
# AUTHOR: Cyriac Kandoth (ckandoth@gmail.com)

use warnings; # Tells Perl to show warnings on anything that might not work as expected
use strict; # Tells Perl to show errors if any of our code is ambiguous
use IO::File; # Helps us read/write files in a safe way
use Getopt::Long qw( GetOptions ); # Helps parse user provided arguments
use Pod::Usage qw( pod2usage ); # Helps us generate nicely formatted help/man content
use JSON::Parse qw( parse_json_safe ); # Helps us parse JSON files
use File::Temp qw( tempfile tempdir ); # Safe way to make temporary files/folders

# Define constants and cutoffs
my ( $min_hom_af, $min_het_af, $max_het_af ) = ( 0.99, 0.20, 0.80 ); # Strict cutoffs to distinguish heterozygous/homozygous sites
my $min_depth = 40; # Don't consider sites with fewer reads than these
my $min_llc_reads = 3; # Minimum variant supporting reads at a site in one sample that could be from low level contamination, where the paired sample is homozygous
my $max_llc_vaf = 0.30; # Maximum VAF at a site in one sample that could be from low level contamination, where the paired sample is homozygous
my $min_llc_sites = 100; # Used for reporting only, minimum number of homozygous sites supporting LLC
my $min_median_llc_vaf = 0.025; # Used for reporting only, minimum median VAF supporting LLC
my $min_t_het_c_hom = 10; # Used for reporting only, minimum number of sites where tumor is het but control is hom i.e. potential T/N swap
use constant { homozyg => 1, var_hom => 2, var_het => 3 };

# Use the CMO JSON to pull paths to tools and data we'll need
my $cmo_cfg_json = `cat /opt/common/CentOS_6-dev/cmo/cmo_resources.json`;
my $cmo_cfg = parse_json_safe( $cmo_cfg_json );
my $r_bin = $cmo_cfg->{programs}->{R}->{default};

# Check for missing or crappy arguments
unless( @ARGV and $ARGV[0]=~m/^-/ ) {
    pod2usage( -verbose => 0, -message => "$0: Missing or invalid arguments!\n", -exitval => 2 );
}

# Parse options and print usage syntax on a syntax error, or if help was explicitly requested
my ( $man, $help ) = ( 0, 0 );
my ( $snp_file, $pairing_file, $readcount_dir, $output_dir );
GetOptions(
    'help!' => \$help,
    'man!' => \$man,
    'snp-file=s' => \$snp_file,
    'pairing-file=s' => \$pairing_file,
    'readcount-dir=s' => \$readcount_dir,
    'output-dir=s' => \$output_dir
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );

# Check that the provided input files exist with a non-zero size, and create the output folder if necessary
die "ERROR: Provided SNP file not found or is empty: $snp_file\n" unless( -s $snp_file );
die "ERROR: Provided T-N pairing file not found or is empty: $pairing_file\n" unless( -s $pairing_file );
mkdir $output_dir unless( -d $output_dir );

# Subroutine to calc median of a given list of numbers
sub median {
    my @values = sort {$a <=> $b} @_;
    my $count = scalar ( @values );
    if( $count % 2 ) {
        return $values[int( $count / 2 )];
    }
    else {
        return( $values[int( $count / 2 ) - 1] + $values[int( $count / 2 )]) / 2;
    }
}

# Load the SNPs whose genotypes will be used for generating our metrics
my %qc_snp = map{chomp; my ($c,$p)=split("\t"); ("$c:$p", 1)}`grep -v ^# $snp_file | cut -f1,2`;

# Locate and read the GATK DepthOfCoverage results for each tumor/control BAM
my %control_readcounts = map {chomp; my ($t,$c)=split("\t"); my ($f)=`find $readcount_dir -name $c*.txt`; chomp($f); ("$t\t$c",$f)} `cat $pairing_file`;
my %tumor_readcounts = map {chomp; my ($t,$c)=split("\t"); my ($f)=`find $readcount_dir -name $t*.txt`; chomp($f); ("$t\t$c",$f)} `cat $pairing_file`;

# Print out column headers into a summary file
mkdir( $output_dir ) unless( -d $output_dir );
my $summary_fh = IO::File->new( "$output_dir/llc_qc_summary.tsv", ">" ) or die "ERROR: Can't open $output_dir/llc_qc_summary.tsv $!";
$summary_fh->print( "Tumor\tControl\tSNP_Sites\tDepth>=$min_depth\tTC_Both_Het\tTC_Matched_Het\tTumor_Hom_Sites\tT_Hom_C_VAF" .
    "\tC_LLC_Sites\tControl_Hom_Sites\tC_Hom_T_VAF\tT_LLC_Sites\tC_Het_T_Hom\tT_Het_C_Hom\tNotes\n" );

# Also create a file where we can store readcounts and VAFs of potentially contaminated positions
my $plot_fh = IO::File->new( "$output_dir/vafs_to_plot.tsv", ">" ) or die "ERROR: Can't open $output_dir/vafs_to_plot.tsv $!";
$plot_fh->print( "SAMPLE\tCHROM\tPOS\tHOM_AF\tVAR_AF\tVAR_READS\n" );

# For each tumor-control pair, print some stats relevant to QC
foreach my $case ( sort keys %control_readcounts ) {
    # Load the readcounts of this bam into a hash
    my ( $control_file, $tumor_file ) = ( $control_readcounts{$case}, $tumor_readcounts{$case} );
    my ( %control_calls, %tumor_calls );
    map {chomp; if( m/^(\w+)[:\t](\d+)\t/ and defined $qc_snp{"$1:$2"} ){ $control_calls{"$1:$2"} = $_ }} `cat $control_file`;
    map {chomp; if( m/^(\w+)[:\t](\d+)\t/ and defined $qc_snp{"$1:$2"} ){ $tumor_calls{"$1:$2"} = $_ }} `cat $tumor_file`;

    my ( %control, %tumor );
    my ( $snp_loci, $loci_checked, $tc_both_het, $matched_het_var, $t_hom_sites, $c_llc_sites, $c_hom_sites, $t_llc_sites, $c_het_t_hom, $t_het_c_hom ) = ( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
    my %counter = ();

    foreach my $locus ( keys %control_calls ) {
        # If not found in both tumor and control bam-readcount files, skip it
        next unless ( defined $tumor_calls{$locus} );
        ++$snp_loci;

        my $line = $control_calls{$locus};
        my ( $chr, $start, $hom_var, $c_depth, $c_depths, $t_depth,  $t_depths );
        # Handle lines from WashU's tool bam-readcount
        if( $line =~ m/^(\S+)\t(\d+)\t\w\t(\d+)\t\S+\t(.*)$/ ){
            ( $chr, $start, $c_depth, $c_depths ) = ( $1, $2, $3, $4 );
            $line = $tumor_calls{$locus};
            ( $t_depth,  $t_depths ) = $line =~ m/^\S+\t\d+\t\w\t(\d+)\t\S+\t(.*)$/;
        }
        # Handle lines from GATK's DepthOfCoverage
        elsif( $line =~ m/^(\S+):(\d+)\t\d+\t\S+\t(\d+)\t(.*)$/ ){
            ( $chr, $start, $c_depth, $c_depths ) = ( $1, $2, $3, $4 );
            $line = $tumor_calls{$locus};
            ( $t_depth,  $t_depths ) = $line =~ m/^\S+:\d+\t\d+\t\S+\t(\d+)\t(.*)$/;
        }
        else { next; }

        next unless( $c_depth >= $min_depth && $t_depth >= $min_depth );
        ++$loci_checked;

        # Get the read-depths of each possible allele at this locus
        my %c_base_depth = map {my @t=split(/:/); (uc($t[0]),$t[1])} split( /\s/, $c_depths );
        my %t_base_depth = map {my @t=split(/:/); (uc($t[0]),$t[1])} split( /\s/, $t_depths );

        # Find variant alleles in the control and tumor, if any
        my ( @c_var, @t_var ) = ((), ());
        foreach my $allele ( sort { $c_base_depth{$b} <=> $c_base_depth{$a} } keys %c_base_depth ) {
            push( @c_var, $allele ) if( $c_base_depth{$allele} > 0 );
        }
        foreach my $allele ( sort { $t_base_depth{$b} <=> $t_base_depth{$a} } keys %t_base_depth ) {
            push( @t_var, $allele ) if( $t_base_depth{$allele} > 0 );
        }

        # The homozygous allele is likely the variant with the highest allele fraction
        $hom_var = ( defined $c_var[0] ? $c_var[0] : ( defined $t_var[0] ? $t_var[0] : '' ));
        # The variant has the next highest allele fraction, that differs from the homozygous allele
        my $c_var = ( defined $c_var[0] and $c_var[0] ne $hom_var ? $c_var[0] : ( defined $c_var[1] and $c_var[1] ne $hom_var ? $c_var[1] : '' ));
        my $t_var = ( defined $t_var[0] and $t_var[0] ne $hom_var ? $t_var[0] : ( defined $t_var[1] and $t_var[1] ne $hom_var ? $t_var[1] : '' ));

        # Find the proportion of reads supporting the variant allele
        my $c_var_reads = (( defined $c_base_depth{$c_var} ) ? $c_base_depth{$c_var} : 0 );
        my $c_vaf = $c_var_reads / $c_depth;
        my $t_var_reads = (( defined $t_base_depth{$t_var} ) ? $t_base_depth{$t_var} : 0 );
        my $t_vaf = $t_var_reads / $t_depth;

        # Find the proportion of reads supporting the homozygous allele
        my $c_raf = (( defined $c_base_depth{$hom_var} ) ? $c_base_depth{$hom_var} : 0 ) / $c_depth;
        my $t_raf = (( defined $t_base_depth{$hom_var} ) ? $t_base_depth{$hom_var} : 0 ) / $t_depth;

        # Now we can figure out whether this site is hom or het in the control and tumor samples
        my ( $c_type, $t_type ) = ( 0, 0 );

        # If either the homozygous or variant allele has the predominant reads, let's call it homozygous
        $c_type = homozyg if( $c_raf >= $min_hom_af );
        $c_type = var_hom if( $c_vaf >= $min_hom_af );
        $t_type = homozyg if( $t_raf >= $min_hom_af );
        $t_type = var_hom if( $t_vaf >= $min_hom_af );

        # Check if this site is heterozygous in either control or tumor
        $c_type = var_het if( $c_vaf > $min_het_af && $c_vaf < $max_het_af );
        $t_type = var_het if( $t_vaf > $min_het_af && $t_vaf < $max_het_af );

        # Count sites with the same heterozygous variant allele in control and tumor
        if( $c_type == var_het && $t_type == var_het ) {
            ++$tc_both_het;
            ++$matched_het_var if( $c_var eq $t_var );
        }

        # If the tumor is homozygous at this site, look for low level contamination in the control
        if( $t_type == homozyg ) {
            ++$t_hom_sites;
            if( $c_var_reads >= $min_llc_reads &&  $c_vaf <= $max_llc_vaf ) {
                push( @{$control{locus}}, $locus );
                push( @{$control{raf}}, $c_raf );
                push( @{$control{vaf}}, $c_vaf );
                push( @{$control{var_reads}}, $c_var_reads );
                ++$c_llc_sites;
            }
        }
        # If the control is homozygous at this site, look for low level contamination in the tumor
        if( $c_type == homozyg ) {
            ++$c_hom_sites;
            if( $t_var_reads >= $min_llc_reads &&  $t_vaf <= $max_llc_vaf ) {
                push( @{$tumor{locus}}, $locus );
                push( @{$tumor{raf}}, $t_raf );
                push( @{$tumor{vaf}}, $t_vaf );
                push( @{$tumor{var_reads}}, $t_var_reads );
                ++$t_llc_sites;
            }
        }

        # Count sites heterozygous in one sample, but homozygous in the matched sample
        # If $t_het_c_hom is high, copy-number losses are happening in the control i.e. a possible TN-swap
        ++$c_het_t_hom if( $c_type == var_het && ( $t_type == homozyg || $t_type == var_hom ));
        ++$t_het_c_hom if( $t_type == var_het && ( $c_type == homozyg || $c_type == var_hom ));
    }

    # Generate stats for this tumor-control pair
    # my $mdn_c_raf  = ( defined $control{raf} ? median( @{$control{raf}} ) * 100 : 0 );
    # my $mdn_t_raf = ( defined $tumor{raf} ? median( @{$tumor{raf}} ) * 100 : 0 );
    my $mdn_c_vaf = ( defined $control{vaf} ? median( @{$control{vaf}} ) : 0 );
    my $mdn_t_vaf = ( defined $tumor{vaf} ? median( @{$tumor{vaf}} ) : 0 );
    $matched_het_var = sprintf( "%.4f", ( $tc_both_het ? ( $matched_het_var / $tc_both_het ) : 0 ));

    # Generate a friendly line of text summarizing any problems with this tumor-control pair
    my $notes = "";
    if( $mdn_c_vaf >= $min_median_llc_vaf and $c_llc_sites >= $min_llc_sites ) {
        $notes .= sprintf( "~%.0f%% LLC in control; ", ( $mdn_c_vaf * 2 * 100 ));
    }
    if( $mdn_t_vaf >= $min_median_llc_vaf and $t_llc_sites >= $min_llc_sites ) {
        $notes .= sprintf( "~%.0f%% LLC in tumor; ", ( $mdn_t_vaf * 2 * 100 ));
    }
    $notes .= "CN losses in control indicate possible sample swap" if( $t_het_c_hom >= $min_t_het_c_hom );

    # Also store these unmatched SNPs with variants into a file for plotting later
    my ( $t_name, $c_name ) = split( "\t", $case );
    for( my $i = 0; $i < ( defined $tumor{locus} ? scalar( @{$tumor{locus}} ) : 0 ); ++$i ) {
        $plot_fh->printf( "%s\t%s\t%i\t%.4f\t%.4f\t%i\n", $t_name, split( ":", ${$tumor{locus}}[$i] ), ${$tumor{raf}}[$i], ${$tumor{vaf}}[$i], ${$tumor{var_reads}}[$i] );
    }
    for( my $i = 0; $i < ( defined $control{locus} ? scalar( @{$control{locus}} ) : 0 ); ++$i ) {
        $plot_fh->printf( "%s\t%s\t%i\t%.4f\t%.4f\t%i\n", $c_name, split( ":", ${$control{locus}}[$i] ), ${$control{raf}}[$i], ${$control{vaf}}[$i], ${$control{var_reads}}[$i] );
    }

    # Add a row to the summary file, with accumulated stats about this TN-pair
    $summary_fh->print( "$case\t$snp_loci\t$loci_checked\t$tc_both_het\t$matched_het_var" );
    $summary_fh->printf( "\t$t_hom_sites\t%.4f\t$c_llc_sites\t$c_hom_sites\t%.4f\t$t_llc_sites", $mdn_c_vaf, $mdn_t_vaf );
    $summary_fh->print( "\t$c_het_t_hom\t$t_het_c_hom\t$notes\n" );
}
$plot_fh->close;
$summary_fh->close;

# Write this R code to a temp file and run it to generate a plot
my $r_code = qq{
    # Make a beanplot of contaminant VAFs per sample:
    library( "beanplot" );
    vafs <- read.delim( "$output_dir/vafs_to_plot.tsv" );
    snps_per_sample <- table( vafs\$SAMPLE );
    plotdata <- split( vafs\$VAR_AF, vafs\$SAMPLE );
    sorted_by_vaf_density <- sort( unlist( lapply( plotdata, sum )));
    plotdata_sorted <- plotdata[names( sorted_by_vaf_density )];
    pdf( "$output_dir/vafs_at_hom_snps.pdf", width=12, height=max( 6, length( plotdata_sorted )/3 ));
    par( mar=c( 5.1, 12.1, 4.1, 2.1 ));
    beanplot( plotdata_sorted, main="Distribution of VAFs at SNPs that are homozygous in a matched sample", xlab="Variant allele fraction (VAF)", what=c( 0, 1, 1, 0 ), horizontal=TRUE, log="", method="stack", grownage=2000, overallline="median", beanlines="median", col=c("darkgrey", NA, NA, "red"), border=c("black"), las=1, ylim=c(0,0.15), xaxt='n' );
    abline( v=c( $min_median_llc_vaf ), lty=2, col="red" );
    axis( side=1, at=seq( 0, 0.2, 0.01 ), cex.lab=0.75 );
    dev.off();
};

my $tmp_dir = tempdir( CLEANUP => 1 );
my ( $r_script_fh, $r_script ) = tempfile( "r_script_XXXX", SUFFIX => '.R', DIR => $tmp_dir );
$r_script_fh->print( $r_code );
my $cmd = "$r_bin/Rscript --vanilla $r_script > /dev/null";
system( $cmd ) == 0 or die "ERROR: Failed to generate vafs_at_hom_snps.pdf! $?";

__DATA__

=head1 NAME

 sample_pair_qc.pl - Look for sample swaps and contamination given pairs of matched samples

=head1 SYNOPSIS

 perl sample_pair_qc.pl --help
 perl sample_pair_qc.pl --snp-file high_het_exac_snps.tsv --pairing-file tn_pairs.tsv --readcount-dir gatk_depth_of_covg --output-dir metrics

=head1 OPTIONS

 --snp-file       TSV file with commonly heterozygous SNP positions. Format: "CHROM POS" or a VCF
 --pairing-file   TSV file with matched sample pairs. Format: "Tumor_Sample_ID Normal_Sample_ID"
 --readcount-dir  Folder with GATK DepthOfCoverage outputs, with file names starting with sample ID
 --output-dir     Path to output folder where metrics and plots will be stored
 --help           Print a brief help message and quit
 --man            Print the detailed manual

=head1 DESCRIPTION

This script generates somatic variant caller commands on each given tumor/normal pair of BAMs, using reasonable defaults for each caller, and prefixed with job scheduler commands (bsub/qsub) using reasonable resource allocation requests. It requires as input a "--pairing-file", which lists the TN-paired sample IDs and their corresponding BAM files in a tab-delimited format like this:

  Tumor_Sample_ID  Normal_Sample_ID

One of the output files is named llc_qc_summary.tsv, containing the following columns:

  The output of this script includes the following columns:
  Tumor, Control - The tumor-vs-control pair whose alleles are being compared for evidence of mismatched samples and/or contamination by exogenous DNA
  SNP_Sites - Number of germline SNP sites (GAF>0.1% in NHLBI exomes) that have at least 1 good read across it, in the sequencing data of both tumor and control
  Depth>=40 - Subset of those germline sites that have at least 40 good reads in both tumor and control
  TC_Both_Het - Sites that are heterozygous in both tumor and control (VAF between 20% and 80%)
  TC_Matched_Het - Percentage of tc_both_het sites that have matching variant alleles (anything <99.9% may indicate that the TC-pair are not from the same patient)
  Tumor_Hom_Sites - Number of sites that are hom in tumor (RAF or VAF >=99%)
  T_Hom_C_VAF - Median control VAF at homozygous sites in tumor (RAF >=99%), that have >=3 reads supporting a variant allele in the control, but <30% control VAF
  C_LLC_Sites - Number of Tumor_Hom_Sites that support low-level contamination in the control as described above
  Control_Hom_Sites - Number of sites that are hom in control (RAF or VAF >=99%)
  C_Hom_T_VAF - Median tumor VAF at homozygous sites in control (RAF >=99%), that have >=3 reads supporting a variant allele in the tumor, but <30% tumor VAF
  T_LLC_Sites - Number of Control_Hom_Sites that support low-level contamination in the tumor as described above
  C_Het_T_Hom - Sites that are het in control (VAF between 20% and 80%) and hom in tumor (RAF or VAF >=99%). A large number of these LOH events in the tumor indicate CN losses, which is not unusual in cancer
  T_Het_C_Hom - Sites that are het in tumor (VAF between 20% and 80%) and hom in control (RAF or VAF >=99%). A large number of these LOH events in the control, relative to the tumor, might indicate an accidental swap between tumor/control
  Notes - Human readable summary of WTF is going on

=head2 Relevant links:

 GATK DepthOfCoverage: https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_coverage_DepthOfCoverage.php

=head1 AUTHORS

 Cyriac Kandoth (ckandoth@gmail.com)

=head1 LICENSE

 Apache-2.0 | Apache License, Version 2.0 | https://www.apache.org/licenses/LICENSE-2.0

=cut
