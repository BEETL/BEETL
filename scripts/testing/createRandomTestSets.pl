#!/usr/bin/perl

## Copyright (c) 2011 Illumina, Inc.
##
## 
## This software is covered by the "Illumina Non-Commercial Use Software
## and Source Code License Agreement" and any user of this software or
## source file is bound by the terms therein (see accompanying file
## Illumina_Non-Commercial_Use_Software_and_Source_Code_License_Agreement.pdf)
##
## This file is part of the BEETL software package.
##
## Citation: Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
## Lightweight BWT Construction for Very Large String Collections. 
## Proceedings of CPM 2011, pp.219-231
##

use strict;
use warnings;

die("usage: ./createRandomTestSets.pl <noOfKmers> <lengthOfKmers> <noOfReads>") unless $#ARGV==2;

my ($noOfKmers,$kmerLength,$noOfReads) = ($ARGV[0],$ARGV[1],$ARGV[2]);

my $outputKmersFn= $noOfKmers . "kmers_" . $kmerLength . "length_" . $noOfReads . "reads.KMERS.txt";
my $outputKmersFasta= $noOfKmers . "kmers_" . $kmerLength . "length_" . $noOfReads . "reads.SEQS.fasta";

open(KMERS,">$outputKmersFn");
open(FASTA,">$outputKmersFasta");


my $readLength=100;

my $nucleotidesToSpend = $readLength-$kmerLength;

my @translation=("A","C","G","T");
my @kmersGenerated;

for( my $i=0;$i<$noOfKmers;$i++ )
{
    my $curKmer = "";
    for( my $j=0;$j<$kmerLength;$j++ )
    {
	$curKmer .= $translation[ int( rand(4) ) ];
    }
    push( @kmersGenerated,$curKmer );
    print KMERS $curKmer . "\n";
}

for( my $i=0;$i<$noOfKmers;$i++ )
{

    for( my $k=0;$k<$noOfReads;$k++ )
    {
	print FASTA ">seq_$i\_$k\n";
	
	my $breakpoint = int( rand($nucleotidesToSpend) ); # we have to split 78 nucleotides on either side
	
	my $curSeq="";
	for( my $j=0;$j<$breakpoint;$j++ )
	{
	    $curSeq .= $translation[ int( rand(4) ) ];
	}
	
	$curSeq .= $kmersGenerated[$i];

	for( my $j=$breakpoint;$j<$nucleotidesToSpend;$j++ )
	{
	    $curSeq .= $translation[ int( rand(4) ) ];
	}

	print FASTA $curSeq . "\n";
	
    }
    
}


close(KMERS);
close(FASTA);
