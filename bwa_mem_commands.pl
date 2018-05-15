#!/usr/bin/env perl

# GOAL: Align FASTQs to a Build38 reference sequence using BWA-MEM. A SampleSheet with
# readgroup info per paired FASTQ is expected, as formatted by the IGO/GCL at MSKCC
#
# AUTHOR: Cyriac Kandoth (ckandoth@gmail.com)

use warnings; # Tells Perl to show warnings on anything that might not work as expected
use strict; # Tells Perl to show errors if any of our code is ambiguous
use IO::File; # A Perl module that can help us read/write files in a safe way

# Root path to srv and opt folders containing software/data dependencies
#my $root_dir = "/cbio/cslab/nobackup";
my $root_dir = "/ifs/e63data/schultzlab";

# Path to the reference FASTA file
my $ref_fa = "$root_dir/opt/vep/homo_sapiens/76_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa";

# Paths to the samtools and bwa binaries
my $samtools_bin = "$root_dir/bin/samtools";
my $bwa_bin = "$root_dir/bin/bwa";

# Make sure we have the necessary number of arguments from the user, and store them into variables for later
( scalar( @ARGV ) == 2 ) or die "Usage: perl $0 <fastq_directory> <output_dir>\n";
my ( $fastq_dir, $output_dir ) = @ARGV;

# Create directories for the log files and resulting bam files
my $log_root = "$output_dir/logs";
my $log_dir = "$log_root/bwa_mem";
my $bam_dir = "$output_dir/bams";
unless( -e $log_root ) { mkdir $log_root or die "Couldn't create directory! $!"; }
unless( -e $log_dir ) { mkdir $log_dir or die "Couldn't create directory! $!"; }
unless( -e $bam_dir ) { mkdir $bam_dir or die "Couldn't create directory! $!"; }

# For each subdir containing paired FASTQ files, construct bwa-mem commands
my @dirs = glob( "$fastq_dir/*" );
foreach my $dir ( @dirs ) {

    # Locate any and all paired FASTQs, and skip this subdir if none are found
    my @r1_fqs = glob( "$dir/*{R1_???,_1}.fastq.gz" );
    my @r2_fqs = glob( "$dir/*{R2_???,_2}.fastq.gz" );
    next unless( @r1_fqs and @r2_fqs );

    # Also skip this subdir, if it doesn't contain a SampleSheet
    my ( $sheet ) = glob( "$dir/*{SampleSheet,samplesheet,sample_sheet}.csv" );
    unless( -s $sheet ) {
        warn "Sample sheet with per-readgroup info not found under $dir ...Skipping...\n";
        next;
    }

    # Parse the SampleSheet, and generate bwa-mem commands for each unique FCID:Lane per sample
    my $in_fh = IO::File->new( $sheet );
    while( my $line = $in_fh->getline ) {

        # Skip header lines, and chomp the others
        next if( $line =~ m/^#|^FCID,/ );
        chomp( $line );

        # FORMAT: FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject,Platform[,Library,InsertSize,Date,Center]
        my ( $flowcell, $lane, $sample_id, undef, $index, $description, undef, $recipe, undef, $project,
            $platform, $library, $insert_size, $date, $center ) = split( /,/, $line );

        # Build a unique readgroup ID and a description for it, based on info in the SampleSheet
        my $rg_ID = $flowcell . ( $lane ? ".$lane" : "" );
        $description = ( $description ? "$description;" : "" ) . ( $project ? $project : "" );
        $description = ( $description ? "$description;" : "" ) . ( $recipe ? $recipe : "" );

        # Figure out which version of bwa-mem we're using
        my $rg_PG = `$bwa_bin 2>&1 | grep Version`;
        $rg_PG =~ s/^Version: (\S+).*\n/bwa-mem $1/;

        # Construct a read-group (@RG) line for the BAM header using info from the SampleSheet
        my $rg_line = "\'\@RG\\tID:$rg_ID";
        $rg_line .= "\\tCN:$center" if( $center );
        $rg_line .= "\\tDS:$description" if( $description );
        $rg_line .= "\\tDT:$date" if( $date );
        $rg_line .= "\\tKS:$index" if( $index );
        $rg_line .= "\\tLB:$library" if( $library );
        $rg_line .= "\\tPG:$rg_PG" if( $rg_PG );
        $rg_line .= "\\tPI:$insert_size" if( $insert_size );
        $rg_line .= "\\tPL:$platform" if( $platform );
        $rg_line .= "\\tPU:$rg_ID" if( $rg_ID );
        $rg_line .= "\\tSM:$sample_id" if( $sample_id );
        $rg_line .= "\'";

        # Skip this SampleSheet line if a BAM was already created for this readgroup
        my $bam_file_prefix = "$bam_dir/$sample_id\_$rg_ID";
        if( -s "$bam_file_prefix.bam" ) {
            warn "BAM already exists for $sample_id ...Skipping...\n";
            next;
        }

        # Find the paired FASTQs for this readgroup (Format: FCID.Lane_* or FCID.Lane.Index_*)
        my ( $r1_fq ) = grep{ m/$flowcell.$lane.*_1.fastq.gz/ } @r1_fqs;
        my ( $r2_fq ) = grep{ m/$flowcell.$lane.*_2.fastq.gz/ } @r2_fqs;

        # If not found, parse the read IDs of each FASTQ, to find those from the same flowcell and lane
        unless( $r1_fq and $r2_fq ) {
            ( $r1_fq ) = grep{ my $fc_ln=`gunzip -c $_ | head -1 | cut -f3,4 -d:`; ($fc_ln eq "$flowcell:$lane\n") } @r1_fqs;
            ( $r2_fq ) = grep{ my $fc_ln=`gunzip -c $_ | head -1 | cut -f3,4 -d:`; ($fc_ln eq "$flowcell:$lane\n") } @r2_fqs;
        }

        # Sometimes the FCID:Lane are the first two values in the read IDs
        unless( $r1_fq and $r2_fq ) {
            ( $r1_fq ) = grep{ my $fc_ln=`gunzip -c $_ | head -1 | cut -f1,2 -d:`; ($fc_ln eq "\@$flowcell:$lane\n") } @r1_fqs;
            ( $r2_fq ) = grep{ my $fc_ln=`gunzip -c $_ | head -1 | cut -f1,2 -d:`; ($fc_ln eq "\@$flowcell:$lane\n") } @r2_fqs;
        }

        # Skip this SampleSheet line if corresponding FASTQs could not be found
        unless( defined $r1_fq and defined $r2_fq ) {
            warn "No paired FASTQs found for $flowcell:$lane (FCID:Lane) under $dir ...Skipping...\n";
            next;
        }

        # Delete the bwa-mem runtime log file if it already exists
        my $log_file = "$log_dir/$sample_id\_$rg_ID.log";
        unlink( $log_file ) if( -e $log_file );

        # Generate a bwa-mem command for this readgroup
        my $bwa_mem_cmd = "$bwa_bin mem -t 4 -M -R $rg_line $ref_fa $r1_fq $r2_fq | $samtools_bin view -bhS - | $samtools_bin sort - $bam_file_prefix";

        # Write the command into a bash script that we can submit to the cluster
        my $script_file = "$log_dir/$sample_id\_$rg_ID.sh";
        my $cmd_fh = IO::File->new( $script_file, ">" );
        $cmd_fh->print( "#!/bin/bash\ntime $bwa_mem_cmd\n" );
        $cmd_fh->close;

        # Print the job submission command for the user to kick off
        # Use this with PBS/Torque:
        #print "qsub -d . -j oe -o $log_file -l nodes=1:ppn=4,mem=8gb,walltime=72:00:00 -N bwa_mem $script_file\n";
        # Use this with SGE:
        print "qsub -q all.q -cwd -S /bin/bash -j y -o $log_file -pe smp 4 -l mem_free=8G,h_vmem=8G,h_stack=512M,virtual_free=8G,h_rt=72:00:00 -N bwa_mem $script_file\n";
    }
    $in_fh->close;
}
