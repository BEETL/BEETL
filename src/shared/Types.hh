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

#include <stdint.h>
#include <string>

using std::string;


// Standard data types
typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long ulong;


// Type to represent: Number of sequences
// below limits to 4 billion reads max - change to uint64_t for more
typedef uint32_t SequenceNumber;

// Type to represent: Sequence length (in biologic case 100)
typedef uint32_t SequenceLength;

// Type to represent: character position or number of characters in BWT
// Should work for BWT of up to 2^64 characters in size
typedef uint64_t LetterNumber;
const LetterNumber maxLetterNumber( static_cast<LetterNumber>( -1 ) );

// Type to represent number of letters in an indexed chunk of the BWT
// Must be the case that max run size per token * max tokens per chunk
// fits into this number without overflow. Might be OK to replace with
// a 16-bit type but we'll leave this at 32 bit initially.
typedef uint32_t LetterNumberCompact;

// For Countwords
const LetterNumber matchFlag( ( ( LetterNumber )1 ) << ( ( 8 * sizeof( LetterNumber ) ) - 1 ) );
const LetterNumber matchMask( ~matchFlag );

// For Metagenomics
const string taxLevelNames[] = {"superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"};
const uint taxLevelSize = 8;
typedef uint32_t MetagFileNumRefType;

// For generalized suffix array (GSA): Definition of each element
struct ElementType
{
    SequenceLength sa;          //It is the position in the sequence, so it goes from 0 a length read
    SequenceNumber numSeq;  //It is the number of the sequence.
};

// For Huffman encoder
union BitBuffer
{
    unsigned int ui;
    unsigned long long ull;
};

// USE_ATTRIBUTE_PACKED: if set, set __attribute__((packed)) for the
// struct sortElement in Sorting.hh; reduces the memory consumption from
// 24 to 17 bytes per input base
#define USE_ATTRIBUTE_PACKED 1


#endif
