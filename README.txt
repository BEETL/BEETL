
README version: $Id: README.txt,v 1.1 2011/11/18 13:14:50 acox Exp $

BEETL: Burrows-Wheeler Extended Tool Library

Copyright (c) 2011 Illumina, Inc.

This software is covered by the "Illumina Non-Commercial Use Software
and Source Code License Agreement" and any user of this software or
source file is bound by the terms therein 

*** Description ***

BEETL is a suite of applications for building and manipulating the 
Burrows-Wheeler Transform (BWT) of collections of DNA sequences. 
The algorithms employed in BEETL are intended to scale to collections of 
sequences containing one billion entries or more.

The initial release implements two flavours of an algorithm for
building the BWT of a sequence collection - BCR and BCRext.

Subsequent releases will add functionality for efficient inversion and
querying of BWTs.


*** System requirements ***

This program has been compiled and tested using gcc version 4.1.2 under 
the CentOS Linux distribution.

BCR requires around 14 bytes of RAM per input sequence 

BCRext works entirely in external memory but places more demands on your file
system. 

In both cases, for best performance it is recommended to run on a 
dedicated local file system. A solid state hard drive (SSD) has been found
to be particularly suitable.

*** Installation ***

Assumes BEETL tar file is in direction /tar_path
Assumes code should be unpacked into directory /install_path

cd /install_path
tar -jxvf /tar_path/localBEETL_0_0_1_alpha.tar.bz2
cd /install_path/beetl/src
make

*** Usage ***
 

1. BCR algorithm (14 bytes of RAM required per input sequence)

/install_path/beetl/Beetl bcr -i input.fa -o outfile

Input: input.fa
FASTA format file, all sequences must be of same size
Only allowed ambiguity code is 'N'

Output:
outfile
- BWT of entire collection in ASCII format
outfilebwt_0.aux 
- BWT of characters corresponding to suffixes beginning with '$'
outfilebwt_1.aux
- BWT of characters corresponding to suffixes beginning with 'A'
outfilebwt_2.aux
- BWT of characters corresponding to suffixes beginning with 'C'
outfilebwt_3.aux
- BWT of characters corresponding to suffixes beginning with 'G'
outfilebwt_4.aux
- BWT of characters corresponding to suffixes beginning with 'N'
outfilebwt_5.aux
- BWT of characters corresponding to suffixes beginning with 'T'
(so outfile is just the concatenation of outfilebwt_[012345].aux)


2. BCRext algorithm (pure external memory)

/install_path/beetl/Beetl ext -i input.txt -a

Input: input.txt
Raw ASCII format, one sequence per line, each terminated by single carriage 
return 
Only allowed ambiguity code is 'N'

Output: all files are in ASCII format

BCRext-B00
- BWT of characters corresponding to suffixes beginning with '$'
BCRext-B01
- BWT of characters corresponding to suffixes beginning with 'A'
BCRext-B02
- BWT of characters corresponding to suffixes beginning with 'C'
BCRext-B03
- BWT of characters corresponding to suffixes beginning with 'G'
BCRext-B04
- BWT of characters corresponding to suffixes beginning with 'N'
BCRext-B05
- BWT of characters corresponding to suffixes beginning with 'T'

/install_path/beetl/Beetl ext -i input.txt -r

as above, but output files are in run-length-encoded format.


*** Release notes ***

Version 0.0.1 (18th November 2011)

This contains initial implementations of the BCR and BCRext algorithms
as described in our CPM paper.

Short-term goals for the next release are:
i. Huffman encoding as an output format for both BCR and BCRext
ii. BCR and BCR ext both take either raw text or fasta files as input
iii. Consistent naming of BCR and BCRext output files
iv. BCR output files in run-length compressed format

*** Contributors ***

Markus J. Bauer, Illumina UK
Anthony J. Cox, Illumina UK (project lead)
Tobias Jakobi, University of Bielefeld
Giovanna Rosone, University of Palermo
Ole Schulz-Trieglaff, Illumina UK

*** Citation ***

Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
Lightweight BWT Construction for Very Large String Collections. 
Proceedings of CPM 2011, pp.219-231

