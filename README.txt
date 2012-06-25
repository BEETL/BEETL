
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

Install:
cd /install_path
tar -jxvf /tar_path/localBEETL_0_0_1_alpha.tar.bz2
cd /install_path/beetl/src
make

Test:
cd /tmp
[or another folder where you don't mind a small number of files being created]
/install_path/test_v_0_0_2.sh

*** Usage ***
 

1. BCR algorithm (14 bytes of RAM required per input sequence)

/install_path/beetl/Beetl bcr -i input.file -o outfile -m 0

Input: input.file
FASTA/FASTQ/raw sequence format file, all sequences must be of same size
Only allowed ambiguity code is 'N'

Output:
outfile
- BWT of entire collection in ASCII format
mybwt0 
- BWT of characters corresponding to suffixes beginning with '$'
mybwt1
- BWT of characters corresponding to suffixes beginning with 'A'
mybwt2
- BWT of characters corresponding to suffixes beginning with 'C'
mybwt3
- BWT of characters corresponding to suffixes beginning with 'G'
mybwt4
- BWT of characters corresponding to suffixes beginning with 'N'
mybwt5
- BWT of characters corresponding to suffixes beginning with 'T'
(so outfile is just the concatenation of outfilebwt_[012345].aux)

/install_path/beetl/Beetl bcr -i mybwt -o inverted.fa -m 1

Input: BWT files as defined above (i.e. output of "-m 0" mode)

Output: Original sequences in FASTA format 
(NB original sequence names are not retained)


2. BCRext algorithm (pure external memory)

/install_path/beetl/Beetl ext -i input.file -p mybwt [-h -r -a] [-s] [-sap]

Input: input.file
FASTA format file, all sequences must be of same size
Only allowed ambiguity code is 'N'
If -s option is set: input in raw ASCII format, one sequence per line, 
each terminated by single carriage return 

Output: all files are in ASCII format

mybwt-B00
- BWT of characters corresponding to suffixes beginning with '$'
mybwt-B01
- BWT of characters corresponding to suffixes beginning with 'A'
mybwt-B02
- BWT of characters corresponding to suffixes beginning with 'C'
mybwt-B03
- BWT of characters corresponding to suffixes beginning with 'G'
mybwt-B04
- BWT of characters corresponding to suffixes beginning with 'N'
mybwt-B05
- BWT of characters corresponding to suffixes beginning with 'T'

-a: output ASCII encoded files
-r: output runlength encoded files [recommended]
-h: output Hufmann encoded files
-sap: perform implicit permutation of collection to obtain more compressible BWT



Note 1: Configuration parameters for BCR mode in Tools.h

convertFromFasta:	if it is set to 1, it reads the input file (FASTA file: a sequence for line), otherwise it reads the cyc files.
deletePartialBWT:	if it is set to 1, it deletes the BWT-segments files and keeps the entire BWT, otherwise renames them.
deleteCycFile:		if it is set to 1, it deletes the cycs files.
BUILD_SA:		if it is set to 1, it computes the GSA (seqID, position) and the SA (position of the concatenated sequences without a further end-marker).
decodeBackward:		if it is set to 1, it computes the inverse BWT in backward direction, otherwise in forward direction.
BackByVector:		if it is set to 1, it uses the sampling of the BWT segments for inverse BWT. More memory, less time.

Note 2: More detail on inverting BWT by using BCR

	/install_path/beetl/Beetl bcr -i BCR-B0 -o outputInverse -m 1

	Input: outfile of BCR

	1) Backward direction
		Inverse BWT by Backward direction of nText sequences at the same time by lengthRead iterations, in the original order.
		1a)	By using the sampling 
		1b) Without using the sampling
		Output: outputInverse. TDB
			For now, it returns the cyc files.
			TDB: To do the inverse of Transpose class. We should obtain from cyc files the sequences in fasta format.

	2) Forward direction
		Inverse BWT by Forward direction of nText sequences, one sequence at a time, in lexicographic order.
		1a)	By using the sampling 
		1b) Without using the sampling
		Output: outputInverse.
			It contains seqID and the recovered sequence for each row.

(Generalized) Suffix Array by using BCR
	
	/install_path/beetl/Beetl bcr -i input.fa -o BCR-B0

	Input: input.fa (like above)
	
	Output: 
		BCR-B0. It is the concatenation of BCR-B0[012345]
		BCR-B0.pairSA. It is the generalized suffix array (seqID, position) of the collection S.
						It is defined as an array of N pairs (seqID, position), and GSA[i]=(t,j) is the the pair corresponding to the i-th smallest suffix of the strings in S.
		BCR-B0.sa. It is the suffix array of the concatenated strings of the collection (without to append a further end-marker)
		BCR-B0.txt. it produces an information file when verboseEncode == 1. Useful for small input.

Example:
input.fa
	> First sequence
	CGAACAGTTA
	> Second sequence
	ACAGTACGAT 

outfile.txt
Position of end-markers
seqN 	 posN 	 pileN
1	3	1
0	3	2

bwt	pos	numSeq	SA
A	10	0	10
T	10	1	21
T	9	0	9
G	2	0	2
$	0	1	11
A	3	0	3
T	5	1	16
C	2	1	13
C	5	0	5
G	8	1	19
A	1	1	12
A	4	0	4
$	0	0	0
A	6	1	17
C	1	0	1
C	7	1	18
A	3	1	14
A	6	0	6
A	9	1	20
T	8	0	8
G	4	1	15
G	7	0	7

-------------------------------------------------------------

*** Release notes ***

Version 0.0.2 (25th June 2012)

Main new feature is an implementation of the compression boosting strategy
described in our Bioinformatics paper [doi: 10.1093/bioinformatics/bts173].
This is the '-sap' flag in BCRext mode.


Other improvements:

Both BCR and BCRext modes can take FASTA or raw sequence format as input

BCR mode can also take FASTQ as input

BCR mode has a '-m 1' submode for inversion of a BWT back to its original
sequences.

Potential gotcha for '-m 1' mode: if using this mode to invert a BWT created 
via the BCRext mode, you will need to create an additional empty file with the
same name as the input name (i.e. 'touch mybwt', if you are using the same 
names as in the above examples)

BCRext mode can produce output in Huffman or run-length encoded format
as well as ASCII. 

BCR mode can now process reads of any length (both modes still demand that
all reads in a dataset are of the same length, however)


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

