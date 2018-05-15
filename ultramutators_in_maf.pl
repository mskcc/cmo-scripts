#!/gsc/bin/perl
## Takes a MAF file, and lists out tumor IDs with more variants than Q3+IQR*4.5
## Author: Cyriac Kandoth
## Last updated: 04.01.2013

use warnings;
use Statistics::Descriptive;

( scalar( @ARGV ) == 1 ) or die "Usage: perl $0 <input_maf>\n";
my ( $maf ) = @ARGV;
( -s $maf ) or die "Provided MAF does not exist or is empty!";

# Count the number of variants per tumor ID
my %var_count;
map{chomp; $var_count{$_}++;}`egrep -v "^(#|Hugo_Symbol)" $maf | cut -f 16`;

# Measure the cutoff for calling a tumor ultramutated
my $stat = Statistics::Descriptive::Full->new();
$stat->add_data( values %var_count );
my $q1 = $stat->quantile( 1 );
my $q3 = $stat->quantile( 3 );
my $cutoff = $q3 + (( $q3 - $q1 ) * 4.5 );

foreach my $tumor_id ( sort keys %var_count ) {
    print "$tumor_id\n" if( $var_count{$tumor_id} > $cutoff );
}

