#!/usr/bin/env perl

# Goal: Measure coverage using Picard's CalculateHSMetrics tool. 
#
# Authors: Cyriac Kandoth (kandothc@mskcc.org), Fanny Dao (daof@mskcc.org)

use warnings; # Tells Perl to show warnings on anything that might not work as expected
use strict; # Tells Perl to show errors if any of our code is ambiguous

# Reference a Perl module that can help us read/write files in a safe way
use IO::File;

# Root path to srv and opt folders
#my $root_dir = "/cbio/cslab/nobackup";
my $root_dir = "/ifs/e63data/schultzlab";

# Hardcode paths to the reference FASTA sequence and Picard dict
my $ref_fa = "$root_dir/srv/Homo_sapiens.GRCh37.73.dna.primary_assembly.reordered.fa";
my $ref_fa_dict = "$root_dir/srv/Homo_sapiens.GRCh37.73.dna.primary_assembly.reordered.dict";

# Make sure we have the necessary number of arguments from the user, and store them into variables for later
( scalar( @ARGV ) == 3 ) or die "Usage: perl $0 <targeted_genes_interval_list> <bam_dir> <output_dir>\n";
my ( $targeted_gene_intervals, $bam_dir, $output_dir   ) = @ARGV;

# Create directories for the resulting logs, interval-list files and metrics files
my $log_dir = "$output_dir/logs/hs_metrics";
my $interval_dir = "$output_dir/gene_intervals";
my $metrics_dir = "$output_dir/metrics";
unless( -e "$output_dir/logs" ) { mkdir "$output_dir/logs" or die "Couldn't create directory! $!"; }
unless( -e $log_dir ) { mkdir $log_dir or die "Couldn't create directory! $!"; }
unless( -e $interval_dir ) { mkdir $interval_dir or die "Couldn't create directory! $!"; }
unless( -e $metrics_dir ) { mkdir $metrics_dir or die "Couldn't create directory! $!"; }

# Create separate interval-list files for each targeted gene. Each interval-list must include the @SQ dictionary entries from the reference fasta
my @genes = map{s/\r?\n$//; $_}`grep -v ^\@ targeted_gene_intervals | cut -f 5 | sort -u`;
foreach my $gene ( @genes ) {
    print `cat $ref_fa_dict > $interval_dir/$gene`;
    print `awk '\$5==\"$gene\"{print}' $targeted_gene_intervals >> $interval_dir/$gene`;
}

# Construct CalculateHsMetrics commands for each BAM and each targeted gene's interval-list
my @bam_files = glob ( "$bam_dir/*.bam" );
chomp ( @bam_files );

# For each BAM, write Picard commands that run sequentially for each gene into a shell script
foreach my $bam ( @bam_files ) {
    my ( $bam_name ) = $bam=~m/\/([^\/]+).bam$/;
    my $cmd_file = "$metrics_dir/run_picard_for_$bam_name.sh";
    my $cmd_fh = IO::File->new( $cmd_file, ">" );
    foreach my $gene ( @genes ) {
        my $gene_interval_file = "$interval_dir/$gene";
        my $output_file = "$metrics_dir/$gene\_$bam_name.tsv";
        $cmd_fh->print( "java -Xmx256m -jar $root_dir/opt/picard/CalculateHsMetrics.jar BI=$gene_interval_file TI=$gene_interval_file I=$bam O=$output_file R=$ref_fa\n" );
    }

    # Write final commands that combines all the per-gene metrics into a single file per BAM, and deletes the per-gene files
    $cmd_fh->print( "egrep -h \"^(#|BAIT_SET|\$)\" $metrics_dir/$genes[0]\_$bam_name.tsv > $metrics_dir/$bam_name.tsv\n" );
    $cmd_fh->print( "egrep -hv \"^(#|BAIT_SET|\$)\" $metrics_dir/*_$bam_name.tsv >> $metrics_dir/$bam_name.tsv\n" );
    $cmd_fh->print( "rm -f $metrics_dir/*_$bam_name.tsv\n" );
    $cmd_fh->close;

    # Submit 1 massive qsub command per bam file
    print "qsub -q all.q -S /bin/bash -V -l mem_free=10G,h_vmem=12G,h_stack=256M -cwd -j y -o $log_dir/$bam_name.log -N hs_metrics $cmd_file\n";
}
