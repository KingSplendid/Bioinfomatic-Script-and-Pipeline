#!/usr/bin/env perl

=head2 NAME

Filter vcf file by bed ...

=head2 SYNOPSIS

perl filter_vcf_by_bed.pl -i <vcf_file> -b <bed_file> -o <filter_vcf_file_path>

=head2 AUTHOR

pangguoqiang 2018-06-4

=cut

################################################################################
#                             Options
################################################################################
use 5.010;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Cwd;
use Data::Dumper;

#parameters
my ( $vcfFile, $bedFile, $filterVcfFilePath );


#help
if ( !@ARGV ) {
	pod2usage( -noperldoc => 1, -verbose => 2 );
	exit(0);
}
Getopt::Long::GetOptions(
	"help|h" => sub { pod2usage( -noperldoc => 1, -verbose => 2 ); exit(0); },
	'i=s'   => \$vcfFile,
	'b=s'   => \$bedFile,
	'o=s'   => \$filterVcfFilePath
);

################################################################################
#                              Main
################################################################################
&logMsg("Starting...");

my %bedInfoHR;
&logMsg("Input bed file ...");
&inputBedFile( $bedFile, \%bedInfoHR );

&logMsg("Filter vcf file ...");
&filterVcfFile( $vcfFile, \%bedInfoHR, $filterVcfFilePath );


&logMsg("Finish ...");
################################################################################
#                           Subroutines
################################################################################


sub filterVcfFile {
	my ( $input, $bedInfoHR, $outputPath ) = @_;

	my $infile = $input;
	my @tmp = split /\//, $input;
	my $filename = $tmp[-1];
	my $outfile = $outputPath . "\/" . $filename;


	open IN, "$infile" || die; 
	open OUT, ">$outfile" || die;

	while(<IN>){
		$_ =~ s/[\r\n]//g;
		if(/^#/){
			print OUT "$_\n";
		}else{
			my @field = split /\t/, $_;
			map( $_ =~ s/^\s+|\s+$//g, @field );
			my $chr = $field[0];
			my $pos = $field[1];
			foreach my $start ( keys $bedInfoHR->{$chr} ){
				my $end = $bedInfoHR->{$chr}->{$start};
				if( $pos > $start && $pos < $end ){
					print OUT "$_\n";
				}
			}

		}
	}
	close IN;
	close OUT;
	return 0;
}

sub inputBedFile {
	# body...
	my ( $input, $bedInfoHR ) = @_;

	open FILE, "$input" || die;
	while(<FILE>){
		next if ( $_ !~ /\S+/ );
        $_ =~ s/[\r\n]//g;
        my @field = split /\t/, $_;

        map( $_ =~ s/^\s+|\s+$//g, @field );
        my $chrNum = $field[0];
        my $start = $field[1];
        my $end = $field[2];
        $bedInfoHR->{ $chrNum }->{ $start } = $end;

	}
	close FILE;
	return 0;
}





sub logMsg {
	my ($string) = @_;
	my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime();
	my $format_time = sprintf( "%d-%d-%d %d:%d:%d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec );
	my $memInfo = `pmap $$ |grep total`;
	$memInfo =~ s/^\s*\S+\s+//;
	$memInfo =~ s/\s+//g;
	print STDERR "\n", $string, "\t( Time: ", $format_time, "\tMem: ", $memInfo, " )\n";

}






