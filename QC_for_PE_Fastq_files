#!/usr/bin/perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

my ($fastqFile1, $fastqFile2, $asc_txt, $outpath);

GetOptions(
  'fq1=s'  => \$fastqFile1,
  'fq2=s'  => \$fastqFile2, 
  'asc=s'  => \$asc_txt,
  'out=s'  => \$outpath
) or pod2usage (2);

if(!$fastqFile1 || !$fastqFile2 || !$asc_txt || !$outpath)
{
	pod2usage ("perl $0 --fq1 <input_fq1> --fq2 <input_fq2> --asc <ASCII.txt> --out <outputpath>");
	exit(0);
}

## Read ASCII.txt
my %asc2shi;
open ASC, "<$asc_txt";
while (<ASC>)
{
	$_ =~ s/[\r\n]//g;
	my @lines = split /\t/, $_;
	$asc2shi{$lines[0]} = $lines[1];
}
close ASC;

my $deleted = $outpath . "/deleted.fq";
my $outfq1  = $outpath . "/clean_1.fq";
my $outfq2  = $outpath . "/clean_2.fq";
my ( $title1, $seq1, $ann1, $score1 ) = ( '', '', '', '' );
my ( $title2, $seq2, $ann2, $score2 ) = ( '', '', '', '' );	
open FQ1,  "<$fastqFile1";
open FQ2,  "<$fastqFile2";
open OUT1, ">$outfq1";
open OUT2, ">$outfq2";
open DEL,  ">$deleted";
while ($title1 = <FQ1>)
{
	$seq1   = <FQ1>;
	$ann1   = <FQ1>;
	$score1 = <FQ1>;
	$title2 = <FQ2>;
	$seq2   = <FQ2>;
	$ann2   = <FQ2>;
	$score2 = <FQ2>;
	##确认数据编号对应
	my $seqId1 = $title1;
	$seqId1 =~ /^(\S+)/;
	$seqId1 = $1;
	my $seqId2 = $title2;
	$seqId2 =~ /^(\S+)/;
	$seqId2 = $1;
	if ( $seqId1 ne $seqId2 )
	{
		print STDERR "error:\t$seqId1 ne $seqId2\n\n";
		exit(1);		
	}
	my $judge1 = &filterQC($title1, $seq1, $ann1, $score1);
	my $judge2 = &filterQC($title2, $seq2, $ann2, $score2);
	if ($judge1 >= 1 || $judge2 >= 1)
	{
		print DEL $title1;
		print DEL $seq1;
		print DEL $ann1;
		print DEL $score1;
		print DEL $title2;
		print DEL $seq2;
		print DEL $ann2;
		print DEL $score2;
	}
	else
	{
		print OUT1 $title1;
		print OUT1 $seq1;
		print OUT1 $ann1;
		print OUT1 $score1;
		print OUT2 $title2;
		print OUT2 $seq2;
		print OUT2 $ann2;
		print OUT2 $score2;
	}
}

close FQ1;
close FQ2;
close OUT1;
close OUT2;
close DEL;


###########################################################################
sub filterQC
{
	my ($title, $seq, $ann, $score) = @_;
	my $judge    = 0;
	my $lowcount = 0;
	$title =~ s/[\r\n]//g;
	$seq   =~ s/[\r\n]//g;
	$ann   =~ s/[\r\n]//g;
	$score =~ s/[\r\n]//g;
	my $seq_length = length($seq);
	my $countforN += ($seq =~ tr/[N|n]//);
	my @scores = split//, $score;
	foreach my $i(@scores)
	{
		my $scorenumber = $asc2shi{$i};
		if ($scorenumber <= 38)
		{
			$lowcount += 1;
		}
	}
	my $percentforN = sprintf("%.2f", $countforN / $seq_length);
	my $lowpercent  = sprintf("%.2f", $lowcount / $seq_length);
	if ($percentforN > 0.1 || $lowpercent > 0.5){$judge += 1;}
	return $judge;
}
