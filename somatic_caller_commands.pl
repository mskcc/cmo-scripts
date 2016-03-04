#!/usr/bin/env perl

# GOAL: Generate somatic variant caller commands on each given tumor/normal pair of BAMs.
# Requires a list of TN-paired sample IDs and their BAM files.
#
# AUTHOR: Cyriac Kandoth (ckandoth@gmail.com)

use warnings; # Tells Perl to show warnings on anything that might not work as expected
use strict; # Tells Perl to show errors if any of our code is ambiguous
use IO::File; # Helps us read/write files in a safe way
use Getopt::Long qw( GetOptions ); # Helps parse user provided arguments
use Pod::Usage qw( pod2usage ); # Helps us generate nicely formatted help/man content
use JSON::Parse qw( parse_json_safe ); # Helps us parse JSON files

# Use the CMO JSON to pull paths to tools and data we'll need
my $cmo_cfg_json = `cat /opt/common/CentOS_6-dev/cmo/cmo_resources.json`;
my $cmo_cfg = parse_json_safe( $cmo_cfg_json );
my $java_bin = $cmo_cfg->{programs}->{java}->{default};
my $samtools_bin = $cmo_cfg->{programs}->{samtools}->{default};
my $varscan_jar = $cmo_cfg->{programs}->{varscan}->{default};
my $mutect_jar = $cmo_cfg->{programs}->{mutect}->{default};
my $gatk_key =  $cmo_cfg->{keys}->{gatk};
my $reference_fa = $cmo_cfg->{genomes}->{GRCh37}->{fasta_vep};
my $dbsnp_vcf = $cmo_cfg->{genomes}->{GRCh37}->{facets_snps};

# Check for missing or crappy arguments
unless( @ARGV and $ARGV[0]=~m/^-/ ) {
    pod2usage( -verbose => 0, -message => "$0: Missing or invalid arguments!\n", -exitval => 2 );
}

# Parse options and print usage syntax on a syntax error, or if help was explicitly requested
my ( $man, $help ) = ( 0, 0 );
my ( $pairing_file, $output_dir );
my ( $normal_purity, $tumor_purity, $job_scheduler ) = ( 0.99, 0.90, "LSF" );
GetOptions(
    'help!' => \$help,
    'man!' => \$man,
    'pairing-file=s' => \$pairing_file,
    'output-dir=s' => \$output_dir,
    'normal-purity' => \$normal_purity,
    'tumor-purity' => \$tumor_purity,
    'job-scheduler' => \$job_scheduler
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );

# Check that the pairing file exists with a non-zero size, and create the output directory if necessary
die "ERROR: Provided T-N pairing file not found or is empty: $pairing_file\n" unless( -s $pairing_file );
mkdir $output_dir unless( -d $output_dir );

# Create directories for the log files and resulting vcf files, unless they already exist
my $log_dir_varscan = "$output_dir/logs/varscan";
my $log_dir_mutect = "$output_dir/logs/mutect";
my $vcf_dir_varscan = "$output_dir/varscan";
my $vcf_dir_mutect = "$output_dir/mutect";
unless( -e "$output_dir/logs" ) { mkdir "$output_dir/logs" or die "Couldn't create directory! $!"; }
unless( -e $log_dir_varscan ) { mkdir $log_dir_varscan or die "Couldn't create directory! $!"; }
unless( -e $log_dir_mutect ) { mkdir $log_dir_mutect or die "Couldn't create directory! $!"; }
unless( -e $vcf_dir_varscan ) { mkdir $vcf_dir_varscan or die "Couldn't create directory! $!"; }
unless( -e $vcf_dir_mutect ) { mkdir $vcf_dir_mutect or die "Couldn't create directory! $!"; }

# For each paired tumor and control, construct and print commands for somatic variant calling
my $paired_sample_fh = IO::File->new( $pairing_file );
while( my $line = $paired_sample_fh->getline ) {
    $line =~ s/\r?\n$//; # This is like chomp, but also works with Windows newlines (\r\n)
    my ( $tumor, $tumor_bam, $normal, $normal_bam ) = split( /\t/, $line );

    # Check if tumor/normal bam files exist
    ( -s "$normal_bam" ) or warn "Normal BAM not found, or is empty!\n$normal_bam\n";
    ( -s "$tumor_bam" ) or warn "Tumor BAM not found, or is empty!\n$tumor_bam\n";
    next unless( -s "$normal_bam" and -s "$tumor_bam" );

    # If the VarScan VCF already exists, print a warning, otherwise print the VarScan command
    if( -e "$vcf_dir_varscan/$tumor\_vs_$normal.snp.vcf" and -e "$vcf_dir_varscan/$tumor\_vs_$normal.indel.vcf" ) {
        warn "VarScan VCF already exists for $tumor\_vs_$normal. Skipping...\n";
    }
    else {
        # Construct samtools mpileup commands for VarScan
        my $samtools_normal = "($samtools_bin mpileup -q 1 -f $reference_fa $normal_bam\)";
        my $samtools_tumor =  "($samtools_bin mpileup -q 1 -f $reference_fa $tumor_bam\)";

        # Delete the log file if it already exists
        my $log_file = "$log_dir_varscan/$tumor\_vs_$normal.log";
        unlink( $log_file ) if( -e $log_file );

        # Generate the VarScan command, using the samtools commands as inputs to the Java command
        my $varscan_cmd = "$java_bin -Xmx6g -Xms512m -jar $varscan_jar somatic <$samtools_normal <$samtools_tumor $vcf_dir_varscan/$tumor\_vs_$normal --min-coverage-normal 8 --min-coverage-tumor 14 --normal-purity $normal_purity --tumor-purity $tumor_purity --strand-filter 1 --output-vcf 1";

        # Write the command into a bash script that we can submit to the cluster
        my $script_file = "$log_dir_varscan/$tumor\_vs_$normal.sh";
        my $cmd_fh = IO::File->new( $script_file, ">" );
        $cmd_fh->print( "#!/bin/bash\n$varscan_cmd\n" );
        $cmd_fh->close;
        chmod 0775, $script_file;

        # Print the job submission command for the user to kick off
        if( $job_scheduler eq "LSF" ) {
            print "bsub -oo $log_file -n 1 -We 96:00 -J varscan $script_file\n";
        }
        elsif( $job_scheduler eq "SGE" ) {
            print "qsub -cwd -S /bin/bash -j y -o $log_file -pe smp 2 -l mem_free=6G,h_vmem=6G,h_stack=512M,virtual_free=6G,h_rt=96:00:00 -N varscan script_file\n";
        }
        elsif( $job_scheduler eq "Torque" ) {
            print "qsub -d . -j oe -o $log_file -l nodes=1:ppn=1,mem=6gb,walltime=96:00:00 -N varscan $script_file\n";
        }
        else {
            print "bash $script_file 2> $log_file\n";
        }
    }

    # If the MuTect VCF already exists, print a warning, otherwise print the MuTect command
    if( -e "$vcf_dir_mutect/$tumor\_vs_$normal.vcf" ) {
        warn "MuTect VCF already exists for $tumor\_vs_$normal. Skipping...\n";
    }
    else {
        # Delete the log file if it already exists
        my $log_file = "$log_dir_mutect/$tumor\_vs_$normal.log";
        unlink( $log_file ) if( -e $log_file );

        # Generate the MuTect command, and write the whole thing into a bash script
        my $mutect_cmd = "$java_bin -Xmx8g -Xms512m -jar $mutect_jar --analysis_type MuTect --phone_home NO_ET --gatk_key $gatk_key --reference_sequence $reference_fa --dbsnp $dbsnp_vcf --input_file:tumor $tumor_bam --input_file:normal $normal_bam --out $vcf_dir_mutect/$tumor\_vs_$normal.txt --vcf $vcf_dir_mutect/$tumor\_vs_$normal.vcf";

        # Write the command into a bash script that we can submit to the cluster
        my $script_file = "$log_dir_mutect/$tumor\_vs_$normal.sh";
        my $cmd_fh = IO::File->new( $script_file, ">" );
        $cmd_fh->print( "#!/bin/bash\n$mutect_cmd\n" );
        $cmd_fh->close;
        chmod 0775, $script_file;

        # Print the job submission command for the user to kick off
        if( $job_scheduler eq "LSF" ) {
            print "bsub -oo $log_file -n 2 -We 72:00 -J mutect $script_file\n";
        }
        elsif( $job_scheduler eq "SGE" ) {
            print "qsub -cwd -S /bin/bash -j y -o $log_file -pe smp 2 -l mem_free=8G,h_vmem=8G,h_stack=512M,virtual_free=8G,h_rt=72:00:00 -N mutect $script_file\n";
        }
        elsif( $job_scheduler eq "Torque" ) {
            print "qsub -d . -j oe -o $log_file -l nodes=1:ppn=2,mem=8gb,walltime=72:00:00 -N mutect $script_file\n";
        }
        else {
            print "bash $script_file 2> $log_file\n";
        }
        
    }
}
$paired_sample_fh->close;

__DATA__

=head1 NAME

 somatic_caller_commands.pl - Generate somatic variant caller commands for given tumor/normal pairs

=head1 SYNOPSIS

 perl somatic_caller_commands.pl --help
 perl somatic_caller_commands.pl --pairing-file tn_pairs.tsv --output-dir somatic_vcfs

=head1 OPTIONS

 --pairing-file   TSV file with: Tumor_Sample_ID Tumor_BAM_Path Normal_Sample_ID Normal_BAM_Path
 --output-dir     Path to output directory where outputs and logs will be kept organized
 --job-scheduler  Job scheduler in use by your cluster (LSF, SGE, Torque, None) [LSF]
 --tumor-purity   Estimated purity of the tumor sample, to use with VarScan [0.90]
 --normal-purity  Estimated purity of the normal sample, to use with VarScan [0.99]
 --help           Print a brief help message and quit
 --man            Print the detailed manual

=head1 DESCRIPTION

This script generates somatic variant caller commands on each given tumor/normal pair of BAMs, using reasonable defaults for each caller, and prefixed with job scheduler commands (bsub/qsub) using reasonable resource allocation requests. It requires as input a "--pairing-file", which lists the TN-paired sample IDs and their corresponding BAM files in a tab-delimited format like this:

  Tumor_Sample_ID  Tumor_BAM_Path Normal_Sample_ID  Normal_BAM_Path

=head2 Relevant links:

 MuTect: https://github.com/mskcc/mutect
 VarScan: http://dkoboldt.github.io/varscan/

=head1 AUTHORS

 Cyriac Kandoth (ckandoth@gmail.com)

=head1 LICENSE

 Apache-2.0 | Apache License, Version 2.0 | https://www.apache.org/licenses/LICENSE-2.0

=cut
