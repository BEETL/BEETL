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

#ifndef DEFINED_CONFIG_HH
#define DEFINED_CONFIG_HH


/* Control flags - modify and recompile to change program behaviour */

//#define DEBUG 1
#define USE_STANDARD_LEXICOGRAPHIC_ORDER 1

// USE_4_BITS_PER_BASE: if set, convert the input sequences from ASCII into
// a 4-bits-per-base format and use that from then on
#define USE_4_BITS_PER_BASE 1


// USE_PREFIX_ONLY: if set, don't copy the whole of the sequences but
// only the prefix that's yet to be processed
// NB - USE_PREFIX_ONLY must only be used with USE_4_BITS_PER_BASE
#define USE_PREFIX_ONLY 1

// TRACK_SEQUENCE_NUMBER: if set, store an unsigned integer with each
// sequence that reports its originating position
//#define TRACK_SEQUENCE_NUMBER 1


// REPORT_COMPRESSION_RATIO: if set, output the compression achieved in the
// partial BWTs (only has an effect if compression is enabled during runtime)
#define REPORT_COMPRESSION_RATIO 1

// REMOVE_TEMPORARY_FILES: if set, remove all files except BWT once creation
// is complete
#define REMOVE_TEMPORARY_FILES 1

/* Control parameters - modify and recompile to change program behaviour */

const int maxSeqSize( 32 * 1024 ); // 32k to fit Moleculo reads
const int bwtBufferSize( 16384 ); // 1<<20=1048576

// Read this many sequences into RAM at once
static const unsigned int ReadBufferSize( 1024 );

// LJ options
//#define SEND_DATA_TO_FILES_FOR_DEBUGGING
//#define READ_DATA_FROM_FILES_FOR_DEBUGGING



// ################## Huffman encoding parameters ###########################

#define ACTIVATE_HUFFMAN
#ifdef ACTIVATE_HUFFMAN

//#define USE_EXTRA_CHARACTER_Z

const int huffmanBufferSize( 128 );
const int huffmanWriterBufferSize( 2097152 );

// number of double codes used to introduce runs
const unsigned int numDoubleCodes( 6 );

// EOT character 'Z'
const unsigned int numSingleCodes( numDoubleCodes + 1 );

// maximal size of one token is 8 bit
const unsigned int maxTokenSize( 8 );

// we will have 256 different tokens
const unsigned int numTokens( 1 << maxTokenSize );

// 255 in binary -> 1111 1111
// used to extract the last 8 bits during processing
const unsigned int tokenMask( numTokens - 1 );

// for once a typedef here, otherwise types.hh would have to include config.hh
// this is a lookup table for fast 8bit mask -> orginal alphabet letter
// i.e. tokentable assigns several 8bit binary tokens to specific letters in
// the input alphabet. A lookup using this table supersedes a slow bit by bit
// comparison deciding which letter is the next.
typedef unsigned int TokenTable[numTokens];

// single letter codes
// hexadecimal codes for each letter used by Huffman encoder
// for reference / debugging below the corresponding binary values
// ACGTN$Z= 011 101 110 100 1110001 11100000

static const unsigned long long singleCharCode[numSingleCodes]
    = {0x17, 0x6, 0x5, 0x3, 0x47, 0x1, 0x7}; // Tony's choice
// = {0x14, 0x6, 0x4, 0x1, 0x2B, 0x3, 0xA8};   // Tobias' try
//  $     A    C    G    N     T    Z
// 10111 011 101   110 1110001 100 11100000
// length in bit for each token
static const unsigned int singleCharLength[numSingleCodes]

//= {5,3,3,3,6,4,8}; // Tobias
    = {5, 3, 3, 3, 7, 3, 8}; // Tony

// without 'Z', runlength encoding codes
// hexadecimal codes for each letter used by Huffman encoder
// for reference / debugging below the corresponding binary values
// ACGTN$ = 111001 1111 001 010 11100001 000

static const unsigned long long doubleCharCode[numDoubleCodes]
//   $     A    C    G     N    T // no Z, double codes
    = { 0x27, 0x0, 0x4, 0x2, 0x87, 0xF}; // Tony's choice
// = { 0xA9, 0x7, 0x2, 0xB, 0x55, 0x0}; // Tobias' choice


// length in bit for each token
static const unsigned int doubleCharLength[numDoubleCodes]
//    $ A C G N T
// = {8,3,3,4,7,3}; // Tobias
    = {6, 3, 3, 3, 8, 4}; // Tony

// position of the notinalphabet char in the baseNames array
static const unsigned int finalCharCode( 6 );

#endif //ifdef ACTIVATE_HUFFMAN

// ################## END Huffman parameters ###########################

#endif
