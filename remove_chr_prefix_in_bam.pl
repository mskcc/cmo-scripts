#!/usr/bin/env perl

use strict;
use warnings;

# Read SAM format from STDIN or from a provided filename
my $skip_warning = 0;
while (<>) {

    # If this is part of the SAM header, then find and replace chrom names in SN tags
    if( m/^\@/ ) {

        # Skip @SQ lines belonging to unplaced contigs and the MT sequence
        if( !m/^\@SQ/ or m/\tSN:(chr)?([0-9]+|[XY])\t/ ) {
            s/\tSN:chr/\tSN:/;
            print;
        }
        else {
            warn( "Skipping MT and the unplaced contigs." ) unless( $skip_warning );
            $skip_warning = 1;
        }
    }

    # For alignment lines, find and replace chrom names in columns 3, 7, and 12 and after
    else {
        my @cols = split( /\t/ );

        # Skip alignments on unplaced contigs and the MT sequence
        if( $cols[2] =~ m/^(chr)?([0-9]+|[XY])$/ ) {
            @cols[2,6,11..$#cols] = map{ s/chrM/MT/g; s/chr//g; $_ } @cols[2,6,11..$#cols];
            $_ = join( "\t", @cols );
            print;
        }
    }
}
