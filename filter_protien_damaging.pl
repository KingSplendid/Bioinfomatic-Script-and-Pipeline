#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
my ($inFile, $outFile);

GetOptions(
'i=s'=>\$inFile,
'o=s'=>\$outFile
) or pod2usage (2);

if (!$inFile || !$outFile){
        pod2usage("please input file!");        
}


&filterProtienDamaging( $inFile, $outFile);

sub filterProtienDamaging{
        my ($infile, $outfile) = @_;
        open FILE, $infile || die;
        open OUT, ">$outfile" || die;
        my $lineIndex = -1;
        my %titleIndex;
        while (<FILE>){
                my @field = split /\t/,$_;

                if ( $lineIndex == -1 ){
                        for (my $i=0;$i<@field;$i++){
                                $titleIndex{ $field[$i] } = $i;        
                        }   
                }
                else {
                        if ($field[ $titleIndex{'SIFT_pred;Polyphen2;MutationTaster_pred'} ]=~/disease|damaging|Deleterious/){
                                print OUT $_;                
                        }        
                }
                $lineIndex++;
                
                
        }
        close FILE;
        close OUT;
}

__END__

=head1 NAME

sample - perl $0 [arguments]

=head1 SYNOPSIS

perl $0 [arguments]

  Optional arguments:

        -i    <String>   input file
        -o    <String>   output file
        

  Program function

        filter protien damaging

        Program example

        perl $0 -i inputfile -o outfile

=cut



