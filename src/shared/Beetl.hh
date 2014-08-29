/**
 ** Copyright (c) 2011-2014 Illumina, Inc.
 **
 ** This file is part of the BEETL software package,
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** Citation: Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
 ** Lightweight BWT Construction for Very Large String Collections.
 ** Proceedings of CPM 2011, pp.219-231
 **
 **/

#ifndef BEETL_HH
#define BEETL_HH

#include <string>

using std::string;


#define COMMAND_BCR_EXT "ext"
#define COMMAND_BCR "bcr"
#define COMMAND_COUNTWORDS "words"


//////////////////////////////////////////////////////////////////////////

// global vars

const string bcrModes[] = {( "build BCR" ), ( "unBCR" ), ( "search Backward search + Locate SeqID" )};

//////////////////////////////////////////////////////////////////////////

// flags for BCRext

// TBD these three are better replaced by an enum defined in BCRext.hh
bool bcrExtAsciiOutput; // use normal ASCII alphabet as output
bool bcrExtHuffmanOutput; // use huffman encoding as compression
bool bcrExtRunlengthOutput; // use RunLength encoding [default]
bool bcrExtUseSeq; // use fasta input
bool bcrExtImplicitSort; // do implicit sort of input sequences


string bcrExtFileIn; // input file with read set (fasta)
string bcrExtFileOutPrefix; // prefix of the output files
const string bcrExtFileOutPrefixDefault = "BCRext"; // default prefix

//////////////////////////////////////////////////////////////////////////

// flags for BCR by Giovanna

/*unsigned*/
short int bcrMode; // default value
// 0 -> build BCR [default]
// 1 -> unBuild BCR
// 2 -> search BCR
string bcrFileIn; // input file
string bcrFileOut; // output file
const string bcrFileOutPrefixDefault = "BCR-B0"; // default prefix



//////////////////////////////////////////////////////////////////////////

// flags for countWords

/*unsigned*/
short int minimalOccurencesN; // maximal number of occurences
/*unsigned*/
short int maxLengthK; // max length

//bool ReferenceGenomeInputB; // set A is a reference genome
char whichHandler;
bool compressedInputA; // input A is compressed
bool compressedInputB; // input B is compressed

string countWordsInputA; // input set A
string countWordsInputB; // input set B

string mergedCFiles; //fileNumbers of merging (one filenumber for one filePosition)
string ncbiInformation; //fileNumber To NCBI Information
bool testDatabase; // flag to indicate that singletons of different lengths in the database should be found, along with the normal metagenome comparison
//////////////////////////////////////////////////////////////////////////



void fileIsReadableOrExit ( string filename );

void print_usage( char *args, const char *command = 0 );

void isArgumentOrExit( int num, int numArgs );

#endif /* BEETL_HH */

