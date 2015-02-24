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

const SequenceNumber maxSequenceNumber( static_cast<SequenceNumber>( -1 ) );
const SequenceLength maxSequenceLength( static_cast<SequenceLength>( -1 ) );
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
//#define METABEETL_FOR_HIV
const string taxLevelNames[] = {"superkingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies", "unnamed rank 9", "unnamed rank 10", "unnamed rank 11"
#ifdef METABEETL_FOR_HIV
, "rank12"
, "rank13"
, "rank14"
, "rank15"
, "rank16"
, "rank17"
, "rank18"
, "rank19"
, "rank20"
, "rank21"
, "rank22"
, "rank23"
, "rank24"
, "rank25"
, "rank26"
, "rank27"
, "rank28"
, "rank29"
, "rank30"
, "rank31"
, "rank32"
, "rank33"
, "rank34"
, "rank35"
, "rank36"
, "rank37"
, "rank38"
, "rank39"
, "rank40"
, "rank41"
, "rank42"
, "rank43"
, "rank44"
, "rank45"
, "rank46"
, "rank47"
, "rank48"
, "rank49"
, "rank50"
, "rank51"
, "rank52"
, "rank53"
, "rank54"
, "rank55"
, "rank56"
, "rank57"
, "rank58"
, "rank59"
, "rank60"
, "rank61"
, "rank62"
, "rank63"
, "rank64"
, "rank65"
, "rank66"
, "rank67"
, "rank68"
, "rank69"
, "rank70"
, "rank71"
, "rank72"
, "rank73"
, "rank74"
, "rank75"
, "rank76"
, "rank77"
, "rank78"
, "rank79"
, "rank80"
, "rank81"
, "rank82"
, "rank83"
, "rank84"
, "rank85"
, "rank86"
, "rank87"
, "rank88"
, "rank89"
, "rank90"
, "rank91"
, "rank92"
, "rank93"
, "rank94"
, "rank95"
, "rank96"
, "rank97"
, "rank98"
, "rank99"
, "rank100"
, "rank101"
, "rank102"
, "rank103"
, "rank104"
, "rank105"
, "rank106"
, "rank107"
, "rank108"
, "rank109"
, "rank110"
, "rank111"
, "rank112"
, "rank113"
, "rank114"
, "rank115"
, "rank116"
, "rank117"
, "rank118"
, "rank119"
, "rank120"
, "rank121"
, "rank122"
#endif //ifdef METABEETL_FOR_HIV
};
#ifdef METABEETL_FOR_HIV
const uint taxLevelSize = 122;
#else
const uint taxLevelSize = 11; // Change this to 122 for HIV
#endif

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
