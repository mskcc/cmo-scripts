#!/usr/bin/env perl

# Author: Cyriac Kandoth
# Date: 02/07/2011
# Desc: Takes a FLAT file with a gene name column and updates it to HUGO names if possible

use strict;
use warnings;
use IO::File;
use LWP::Simple;
use FindBin;
use lib "$FindBin::Bin/lib";
use HugoGene;
use HugoGeneMethods;

if( scalar( @ARGV ) != 2 )
{
  print STDERR "\nUsage: perl $0 <tab_delimited_input_file> <gene_name_column_index>";
  print STDERR "\nThe gene_name_column_index must be a zero-based index\n\n";
  exit 1;
}

#Create hashes that resolve latest HUGO names
my $HugoObjRef = HugoGeneMethods::makeHugoGeneObjects();
my ( $PreviousSymbolToCurrentRef, undef ) = HugoGeneMethods::previousHugoSymbols( $HugoObjRef );
my ( $AliasToCurrentRef, undef ) = HugoGeneMethods::hugoAlias( $HugoObjRef );

my $flat_file = $ARGV[0];
my $gene_column = $ARGV[1];

my $annotFh = IO::File->new( $flat_file ) or die "Can't open $flat_file $!";
my $newFh = IO::File->new( "$flat_file\_hugoified", ">" ) or die "Can't open $flat_file\_hugoified $!";
my $ambigFh = IO::File->new( "$flat_file\_ambiguous", ">" ) or die "Can't open $flat_file\_ambiguous $!";
my ( %undefCnt, %defCnt, %sameCnt ) = ((), (), ());
while( my $line = $annotFh->getline )
{
  if( $line =~ m/^#/ )
  {
    $newFh->print( $line ); #Copy any commented lines to the new file
    next;
  }

  my @segs = split( /\t/, $line );
  chomp( $segs[-1] );
  my $name = $segs[$gene_column];

  #Searching for it as all uppercase, makes it case insensitive
  my $hugo = getHugoSymbol( uc( $name ));
  if( defined $hugo )
  {
    ++$defCnt{$hugo} if ( uc( $hugo ) ne uc( $name ));
    $segs[$gene_column] = $hugo;
  }
  else
  {
    $ambigFh->print( "Couldn't find a unique HUGO name for: $name\n" );
    ++$undefCnt{$name};
  }

  #Write the modified (or unchanged) line to file
  $newFh->print( join( "\t", @segs ), "\n" );
}
$ambigFh->close;
$newFh->close;
$annotFh->close;

print scalar( keys %defCnt ), " gene names needed converting to HUGO\n";
print scalar( keys %undefCnt ), " gene names could not be resolved\n";
print "Unresolved names with more than two variants:\n";
foreach my $key ( sort {$undefCnt{$b} <=> $undefCnt{$a}} ( keys %undefCnt ))
{
  if( $undefCnt{ $key } > 2) { print $undefCnt{ $key }, "\t$key\n"  };
}

sub getHugoSymbol
{
  my $oldSymbol = $_[0];
  #If it's already the latest Hugo name
  if( defined $$HugoObjRef{$oldSymbol} )
  {
    return $$HugoObjRef{$oldSymbol}->symbol();
  }

  #If it's a previous symbol of a unique Hugo name
  if( defined $$PreviousSymbolToCurrentRef{$oldSymbol} &&
         scalar( @{$$PreviousSymbolToCurrentRef{$oldSymbol}} ) == 1 )
  {
    return ${$$PreviousSymbolToCurrentRef{$oldSymbol}}[0];
  }

  #If it's an alias of a unique Hugo name
  if( defined $$AliasToCurrentRef{$oldSymbol} &&
         scalar( @{$$AliasToCurrentRef{$oldSymbol}} ) == 1 )
  {
    return ${$$AliasToCurrentRef{$oldSymbol}}[0];
  }

  #Sorry mate. Couldn't figure it out :-|
  return undef;
}
