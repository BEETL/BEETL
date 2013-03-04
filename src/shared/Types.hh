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

#ifndef DEFINED_TYPES_HH
#define DEFINED_TYPES_HH

#include <string>

using std::string;


/* Standard data types */

// defined in Tool.h by giovanna, commented to
// catch compile time error for umtiple typedefs

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long ulong;

#ifdef USE_TYPEDEF_NOT_HASH_DEFINE
#ifndef uchar
#define uchar unsigned char
#endif
#ifndef uint
#define uint unsigned int
#endif
#ifndef ulong
#define ulong unsigned long
#endif
#endif

// used for huffman encoder
union BitBuffer
{
    unsigned int ui;
    unsigned long long ull;
};

// below limits to 4 billion reads max - change to unsigned long for more
typedef unsigned int SequenceNumberType;

// Should work for BWT of up to 2^64 characters in size
typedef unsigned long long LetterCountType;

const LetterCountType maxLetterCountType( static_cast<LetterCountType>( -1 ) );

const LetterCountType matchFlag( ( ( LetterCountType )1 )<<( ( 8*sizeof( LetterCountType ) )-1 ) );
const LetterCountType matchMask( ~matchFlag );

const string taxLevel[] = {"superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"};
const uint taxLevelSize = 8;

#endif
