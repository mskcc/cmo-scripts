#!/usr/bin/env perl

use strict;
use warnings;
use IO::File;

if( scalar( @ARGV ) != 2 ) {
    print STDERR "\nUsage: perl $0 <bed_file> <wig_file>\n";
    exit 1;
}
my $bed_file = $ARGV[0];
my $wig_file = $ARGV[1];

# Create a new WIG file using the re-mapped loci
`joinx sort -s $bed_file -o $bed_file\_sorted`;
my $ifh = IO::File->new( "$bed_file\_sorted" ) or die "Couldn't open $bed_file\_sorted for reading! $!";
my $ofh = IO::File->new( $wig_file, ">" ) or die "Couldn't open $wig_file for reading! $!";
$ofh->print( "track type=wiggle_0 name=SomaticCoverage viewLimits=0:1\n" );
my ( $x_chr, $x_start, $x_stop ) = ( 0, 0, 0 );
while( my $line = $ifh->getline ) {
    chomp( $line );
    my ( $chr, $start, $stop ) = split( /\t/, $line );
    $ofh->print( "fixedStep chrom=$chr start=" . ( $start + 1 ) . " step=1 span=" . ( $stop - $start ) . "\n" );
    $ofh->print( "1\n" );
    ( $x_chr, $x_start, $x_stop ) = ( $chr, $start, $stop );
}
$ofh->close;
$ifh->close;

unlink( "$bed_file\_sorted" );
