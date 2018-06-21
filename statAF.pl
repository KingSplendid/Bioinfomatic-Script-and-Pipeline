#!/usr/bin/env perl

=head2 NAME

Statisitc AF for all samples. v0.1

=head2 SYNOPSIS

perl statAF.pl [vcfPath] [depthPath] [outputPath]

Note:
 1. vcf file with "GT:AD:DP:GQ:PL" info;
 

=head2 AUTHOR

Zhigang Li		 2016-8-17

=cut

################################################################################
#                             Options
################################################################################
use 5.010;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

#parameters
my ( $vcfPath, $depthPath, $outputPath ) = @ARGV;
our $vcfDepthCutoff     = 10;
our $genotypeFreqCutoff = 0.2;    # if seqFreq less the cutoff, this allele is not exists in our consider.
our $depthCutoff        = 10;

# combine some type
# typeOrigin => typeSubject
our %typeCombine = ();

#help
if ( !@ARGV ) {
        pod2usage( -noperldoc => 1, -verbose => 2 );
        exit(0);
}
Getopt::Long::GetOptions( "help|h" => sub { pod2usage( -noperldoc => 1, -verbose => 2 ); exit(0); }, );

################################################################################
#                              Main
################################################################################
&logMsg('Starting...');

&checkFold( $vcfPath, $depthPath );
&createFold($outputPath);


&logMsg("Process VCF file...");

# get all files in vcf path
my @vcfFiles;
&getFileInFold( $vcfPath, \@vcfFiles );

# all allele infor
# e.g. \%alleleInfo = {
#          'ALSCMT' => {
#                        'chr22' => {
#                                     '29885917' => {
#                                                     'rsId' => {
#                                                                 '.' => 1
#                                                               },
#                                                     'altAllele' => {
#                                                                      'T' => 1
#                                                                    },
#                                                     'refAllele' => 'C'
#                                                   },
#                        }
#          }
my %alleleInfo;

# contain vcf file, chr,posi;
# if the position in the sample in exists in vcf file, it's not negative and pass in depth filter step
# e.g. \%positivePosiInSample = {
#          '24:chr15:23086365' => 1,
#          '20:chr2:241664801' => 1,
# }
my %positivePosiInSample;


# default AN number for each type
# e.g. \%defaultANNum = {
#          'ALSCMT' => 102
#        };

my %defaultANNum;

# each sample
# calculate all AF info
my @vcfFileUsed;    # contain all vcf file
my $sampleNum = 0;  # the subscript of @vcfFileUsed
foreach my $fileName (@vcfFiles) {
        next if ( $fileName !~ /\.[vV][cC][fF]$/ );
        push @vcfFileUsed, $fileName;

        my $type = &getType($fileName);
        $defaultANNum{$type} += 2;
        &logStd("\t$sampleNum\t$type\t$fileName\n");
        my $vfile = $vcfPath . '/' . $fileName;
        open VCF, "$vfile" || die;

        while (<VCF>) {
                next if ( $_ =~ /^#/ );
                next if ( $_ !~ /\S/ );
                &processVcfLine( $_, $type, $sampleNum, \%alleleInfo, \%positivePosiInSample );
        }
        close VCF;
        $sampleNum++;
}

############################################################
# filter false negative  based on depth file
&logMsg("Process Depth file...");
$sampleNum = 0;    # the subscript of @vcfFileUsed
foreach my $vf (@vcfFileUsed) {
        my $type    = &getType($vf);
        my $depName = $vf;
        $depName =~ s/\.[^.]+$//;
        $depName .= '.depth';
        my $depFile = $depthPath . '/' . $depName;

        if ( -e $depFile ) {
                my $time = &timer();
                &logStd("\t$sampleNum\t $type\t$depName\t($time)\n");

                # get low quality position
                # e.g %lowQuaPosi = {
                # 'chr15' => {
                #       '44880837' => 1,
                #       '44880846' => 1,
                #       '44880822' => 1,
                # }
                #}
                my %lowQuaPosi;    # contain low quality position info for one sample
                &getLowQualityPosiInDepth( $depFile, \%lowQuaPosi );

                # filter false negative; change AN number
                &filterFN( \%lowQuaPosi, \%alleleInfo, $type, $defaultANNum{$type}, $sampleNum, \%positivePosiInSample );
        }
        else {
                &logWarn( "DepthFileNotExistsFor: $vf", __FILE__, __LINE__ );
        }
        $sampleNum++;
}

############################################################
# output
&logMsg('Output...');
foreach my $type ( keys %alleleInfo ) {
        my $outName = $type . '_af.txt';
        my $outFile = $outputPath . '/' . $outName;

        &logStd("\t$outName\n");
        open OF, ">$outFile" || die;
        print OF "#Chr\tPosi\tRsId\tRef\tAlt\tAC\tAN\tAF\tComment\tType\n";

        # each chr
        foreach my $chr ( keys %{ $alleleInfo{$type} } ) {

                # each posi
                foreach my $posi ( keys %{ $alleleInfo{$type}->{$chr} } ) {
                        my $rsId = join( ',', keys %{ $alleleInfo{$type}->{$chr}->{$posi}->{'rsId'} } );
                        my $ref = $alleleInfo{$type}->{$chr}->{$posi}->{'refAllele'};
                        my $AN =
                            ( exists $alleleInfo{$type}->{$chr}->{$posi}->{'AN'} )
                          ? ( $alleleInfo{$type}->{$chr}->{$posi}->{'AN'} )
                          : ( $defaultANNum{$type} );
                        my %altInfo = %{ $alleleInfo{$type}->{$chr}->{$posi}->{'altAllele'} };

                        my $num     = 0;     # alt sorted allele number
                        my $alt     = '';    # only show the top three alt allele
                        my $AF      = '';    # only show the top three alt allele
                        my $AC      = '';    # only show the top three alt allele
                        my $comment = '';    # show other alt allele info

                        # each alt allele
                        foreach my $allele ( sort { $altInfo{$b} <=> $altInfo{$a} } keys %altInfo ) {

                                #print STDERR "$chr\t$posi\t$allele\t$AN\t", $defaultANNum{$type}, "\t", $type, "\n";

                                #print STDERR Dumper(\%defaultANNum);
                                $num++;
                                if ( $num < 3 ) {
                                        $alt .= $allele . ',';
                                        $AC  .= $altInfo{$allele} . ',';
                                        $AF  .= sprintf( "%.6f", $altInfo{$allele} / $AN ) . ',';
                                }
                                else {
                                        $comment = $allele . ':' . $altInfo{$allele} . ':' . sprintf( "%.6f", $altInfo{$allele} / $AN ) . ',';
                                }
                        }
                        $alt =~ s/,$//;
                        $AF =~ s/,$//;
                        $AC =~ s/,$//;
                        $comment =~ s/,$//;

                        # output
                        print OF $chr, "\t", $posi, "\t", $rsId, "\t", $ref, "\t";
                        print OF $alt, "\t", $AC, "\t", $AN, "\t", $AF, "\t", $comment, "\t", $type, "\n";
                }
        }
        close OF;
}

&logMsg('Finishing...');

################################################################################
#                           Subroutines
################################################################################

# Title   :
# Function:
# Usage   :
# Returns :
# Args    :

sub createFold {
        my ($fold) = @_;
        if ( -d $fold ) {
                system("rm -rf $fold");
        }
        mkdir $fold;
        return 0;
}

# Title   :
# Function:
# Usage   :
# Returns :
# Args    :

sub checkFold {
        my @folds = @_;
        foreach (@folds) {
                if ( !-d $_ ) {
                        &logError( "FoldNotExists: $_", __FILE__, __LINE__ );
                }
        }
        return 0;
}

# Title   :
# Function:
# Usage   :
# Returns :
# Args    :

sub getType {
        my ($name) = @_;
        my @field = split /_/, $name;
        if ( exists $typeCombine{ $field[2] } ) {
                return $typeCombine{ $field[2] };
        }
        else {
                return $field[2];
        }

}

# Title   :
# Function:
# Usage   :
# Returns :
# Args    :

sub filterFN {
        my ( $lowQuaPosiHR, $alleleInfoHR, $type, $anNumDefault, $sampleNum, $positivePosiInSampleHR ) = @_;
        foreach my $chr ( %{$lowQuaPosiHR} ) {
                foreach my $posi ( %{ $lowQuaPosiHR->{$chr} } ) {

                        # if the position in the sample in exists in vcf file, it's not negative and pass
                        next if ( exists $positivePosiInSampleHR->{ $sampleNum . ':' . $chr . ':' . $posi } );

                        # filter false negative
                        if ( exists $alleleInfoHR->{$type}->{$chr}->{$posi} ) {
                                if ( !exists $alleleInfoHR->{$type}->{$chr}->{$posi}->{'AN'} ) {
                                        $alleleInfoHR->{$type}->{$chr}->{$posi}->{'AN'} = $anNumDefault;
                                }

                                # delete FN
                                $alleleInfoHR->{$type}->{$chr}->{$posi}->{'AN'} -= 2;

                        }
                }
        }
        return 0;
}

# Title   :
# Function:
# Usage   :
# Returns :
# Args    :

sub getLowQualityPosiInDepth {
        my ( $file, $lowQuaPosiHR ) = @_;
        open DEP, "$file" || die;
        <DEP>;
        while (<DEP>) {
                $_ =~ /^(\S+):(\d+)\s+(\d+)/;
                if ( $3 < $depthCutoff ) {
                        $lowQuaPosiHR->{$1}->{$2} = 1;
                }
        }
        close DEP;
        return 0;
}

# Title   :
# Function:
# Usage   :
# Returns :
# Args    :

sub processVcfLine {
        my ( $line, $type, $sampleNum, $alleleInfoHR, $positivePosiInSampleHR ) = @_;
        $line =~ s/[\r\n]//g;

        if ( $line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t[\s\S]+GT:AD:DP:GQ:PL\t(\S+)/ ) {
                my ( $chr, $posi, $rsId, $refAllele, $altAllele, $depthInfo ) = ( $1, $2, $3, $4, $5, $6 );

                # all base of allele in this position
                my @baseAllele = ($refAllele);
                push @baseAllele, split( /,/, $altAllele );

                # parsing depth info
                my @depInfoField = split /:/, $depthInfo;
                my $dp           = $depInfoField[2];                    # depth of position, DP
                my @genotype     = split /[\|\/]/, $depInfoField[0];    # GT

                # filter based on sequence depth
                if ( $dp >= $vcfDepthCutoff ) {

                        # modify posi, ref allele, alt allele for each alt allele
                        my %modifyPosiAllele;
                        foreach my $g (@genotype) {
                                next if ( $g == 0 );

                                # modify these value if necessary
                                # check whethe necessary
                                my ( $posiChecked, $refAlleleChecked, $altAlleleChecked );
                                if ( length($refAllele) != 1 || length( $baseAllele[$g] ) != 1 ) {

                                        ( $posiChecked, $refAlleleChecked, $altAlleleChecked ) = &transformAllele( $posi, $refAllele, $baseAllele[$g] );
                                }
                                else {
                                        ( $posiChecked, $refAlleleChecked, $altAlleleChecked ) = ( $posi, $refAllele, $baseAllele[$g] );
                                }

                                # check
                                if ( $refAlleleChecked eq $altAlleleChecked ) {
                                        &logError( "refSameAsAltAllele:\nerf: $refAllele\nalt: $altAllele\n$line\n", __FILE__, __LINE__ );
                                }

                                # after modify
                                $modifyPosiAllele{$g} = [ $posiChecked, $refAlleleChecked, $altAlleleChecked ];

                                # store positive file, chr, posi info
                                $positivePosiInSampleHR->{ $sampleNum . ':' . $chr . ':' . $posiChecked } = 1;
                        }

                        # change genotype based on our criteria
                        my @genotypeModify = &modifyGenotype($depthInfo);
                        foreach my $g (@genotypeModify) {
                                next if ( $g == 0 );

                                # True positive, as AC
                                my ( $posiChecked, $refAlleleChecked, $altAlleleChecked ) = @{ $modifyPosiAllele{$g} };

                                # process and store position info
                                if ( exists $alleleInfoHR->{$type}->{$chr}->{$posiChecked} ) {
                                        my $refAlleleStorage = $alleleInfoHR->{$type}->{$chr}->{$posiChecked}->{'refAllele'};
                                        if ( $refAlleleChecked ne $refAlleleStorage ) {

                                                &logError( "refAlleleDifferent:\n$refAlleleStorage\n$refAllele\n$line\n", __FILE__, __LINE__ );
                                        }
                                        else {
                                                $alleleInfoHR->{$type}->{$chr}->{$posiChecked}->{'altAllele'}->{$altAlleleChecked}++;
                                                $alleleInfoHR->{$type}->{$chr}->{$posiChecked}->{'rsId'}->{$rsId} = 1;
                                        }

                                }
                                else {
                                        $alleleInfoHR->{$type}->{$chr}->{$posiChecked}->{'refAllele'}                      = $refAlleleChecked;
                                        $alleleInfoHR->{$type}->{$chr}->{$posiChecked}->{'altAllele'}->{$altAlleleChecked} = 1;
                                        $alleleInfoHR->{$type}->{$chr}->{$posiChecked}->{'rsId'}->{$rsId}                  = 1;
                                }
                        }
                }
        }
        else {
                &logWarn( "VcfFormatParsingError: $line", __FILE__, __LINE__ );
        }
        return 0;
}

# Title   :
# Function:
# Usage   :
# Returns :
# Args    :

sub modifyGenotype {
        my ($info) = @_;
        my @field    = split /:/, $info;
        my $dp       = $field[2];                    # depth of position, DP
        my @genotype = split /[\|\/]/, $field[0];    # GT
        my @ad       = split /,/, $field[1];         # AD

        &logError( "Genotype more than 2: $info\n", __FILE__, __LINE__ ), if ( @genotype > 2 );
        if ( $genotype[0] == $genotype[1] ) {
                if ( $ad[ $genotype[0] ] / $dp < $genotypeFreqCutoff ) {
                        @genotype = ( 0, 0 );
                }
        }
        else {
                if ( $ad[ $genotype[0] ] / $dp < $genotypeFreqCutoff ) {
                        $genotype[0] = $genotype[1];
                }
                elsif ( $ad[ $genotype[1] ] / $dp < $genotypeFreqCutoff ) {
                        $genotype[1] = $genotype[0];
                }
        }

        return @genotype;
}

# Title   :
# Function: transform allele format. the position is the same as annvoar.
#           But for insertion, the ref allele was included in the alt allele, e.g. C -> CG. Others were the same as annovar.
# Usage   :
# Returns :
# Args    :

sub transformAllele {
        my ( $posi, $refAllele, $altAllele ) = @_;

        #print STDERR $posi,"\t",$refAllele, "\t", $altAllele, "\n";
        my $refLen  = length($refAllele);
        my $altLen  = length($altAllele);
        my $diffLen = abs( $refLen - $altLen );

        # make length between ref and alt same
        if ( $refLen > $altLen ) {
                $altAllele .= '-' x $diffLen;
        }
        elsif ( length($refAllele) < length($altAllele) ) {
                $refAllele .= '-' x $diffLen;
        }

        my @refBase = split //, $refAllele;
        my @altBase = split //, $altAllele;
        my $startBasePosi = 0;
        my $endBasePosi   = $#refBase;

        # delete from left
        my $privousRefBase = '';    # the 5' most previous base
        if ( $refBase[0] eq $altBase[0] ) {
                for ( $startBasePosi = 0 ; $startBasePosi < @refBase ; $startBasePosi++ ) {
                        if ( $refBase[$startBasePosi] eq $altBase[$startBasePosi] ) {
                                $posi++;
                                $privousRefBase = $refBase[$startBasePosi];
                        }
                        else {
                                last;
                        }
                }
        }

        # delete from right
        if ( $refBase[$#refBase] eq $altBase[$#altBase] ) {
                for ( $endBasePosi = $#refBase ; $endBasePosi >= 0 ; $endBasePosi-- ) {
                        last if ( $refBase[$endBasePosi] ne $altBase[$endBasePosi] );
                }
        }

        $refAllele = join( '', @refBase[ $startBasePosi .. $endBasePosi ] );
        $altAllele = join( '', @altBase[ $startBasePosi .. $endBasePosi ] );

        # check whether the event is insertion
        if ( $refAllele =~ /^-/ ) {
                $posi--;
                $refAllele = $privousRefBase;
                $altAllele = $privousRefBase . $altAllele;
        }
        else {
                $refAllele = $refBase[$startBasePosi];
        }

        #print STDERR $posi,"\t", $refAllele, "\t", $altAllele, "\n---------\n";

        return ( $posi, $refAllele, $altAllele );

}

# Title   :
# Function:
# Usage   :
# Returns :
# Args    :

sub getFileInFold {
        my ( $path, $filesAR ) = @_;
        opendir( CPATH, $path );
        @{$filesAR} = readdir(CPATH);
        closedir(CPATH);

        return 0;

}

# some common sub-function
sub timer {
        my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime();
        my $format_time = sprintf( "%d-%d-%d %d:%d:%d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec );
        return $format_time;

}

sub logStd {
        my ($string) = @_;
        print STDERR "$string";
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

sub logError {
        my ( $string, $fileName, $lineNum ) = @_;
        print STDERR "!! [ ERROR ] $fileName\tline: $lineNum\n";
        print STDERR "!! $string\n";
        print STDERR "!! ** Exit 1 **\n";
        exit 1;
}

sub logWarn {
        my ( $string, $fileName, $lineNum ) = @_;
        print STDERR "!! [ WARNING ] $string $fileName\tline: $lineNum\n";
}
