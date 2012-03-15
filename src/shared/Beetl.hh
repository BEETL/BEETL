/**
 ** Copyright (c) 2011 Illumina, Inc.
 **
 ** 
 ** This software is covered by the "Illumina Non-Commercial Use Software
 ** and Source Code License Agreement" and any user of this software or
 ** source file is bound by the terms therein (see accompanying file
 ** Illumina_Non-Commercial_Use_Software_and_Source_Code_License_Agreement.pdf)
 **
 ** This file is part of the BEETL software package.
 **
 ** Citation: Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
 ** Lightweight BWT Construction for Very Large String Collections. 
 ** Proceedings of CPM 2011, pp.219-231
 **
 **/

#include <string>


#ifndef BEETL_HH
#define	BEETL_HH

#define COMMAND_BCR_EXT "ext"
#define COMMAND_BCR "bcr"
#define COMMAND_COUNTWORDS "words"

using namespace std;


//////////////////////////////////////////////////////////////////////////

// global vars

const string bcrModes[] = {("build BCR"),("unBCR"),("search Backward search + Locate SeqID")};

//////////////////////////////////////////////////////////////////////////

// flags for BCRext

bool bcrExtAsciiOutput; // use normal ASCII alphabet as output
bool bcrExtHuffmanOutput; // use huffman encoding as compression
bool bcrExtRunlengthOutput; // use RunLength encoding [default]
bool bcrExtUseSeq; // use fasta input

string bcrExtFileIn; // input file with read set (fasta)
string bcrExtFileOutPrefix; // prefix of the output files
const string bcrExtFileOutPrefixDefault = "BCRext"; // default prefix

//////////////////////////////////////////////////////////////////////////

// flags for BCR by Giovanna

/*unsigned*/ short int bcrMode; // default value
// 0 -> build BCR [default]
// 1 -> unBuild BCR
// 2 -> search BCR
string bcrFileIn; // input file
string bcrFileOut; // output file


//////////////////////////////////////////////////////////////////////////

// flags for countWords

/*unsigned*/ short int minimalOccurencesN; // maximal number of occurences
/*unsigned*/ short int minimalLengthK; // minimal length 

bool ReferenceGenomeInputB; // set A is a reference genome 
bool compressedInputA; // input A is compressed
bool compressedInputB; // input B is compressed
bool compressedBoth; // A & B are compressed

string countWordsInputA; // input set A 
string countWordsInputB; // input set B

//////////////////////////////////////////////////////////////////////////  



void fileIsReadableOrExit (string filename);

void print_usage(char *args);

void isArgumentOrExit(int num, int numArgs);

#endif	/* BEETL_HH */

