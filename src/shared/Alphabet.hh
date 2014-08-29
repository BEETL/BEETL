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

#ifndef DEFINED_ALPHABET_HH
#define DEFINED_ALPHABET_HH

#include "Config.hh"
#include "Types.hh"


// Convention is that first char is the string terminator char and
// is lexicographically less than the rest

#ifdef USE_STANDARD_LEXICOGRAPHIC_ORDER
# ifndef USE_EXTRA_CHARACTER_Z

static const char alphabet[] = "$ACGNT";
static const int whichPileInverse[] =
{
    '$', 'A', 'C', 'G', 'N', 'T'
};
enum
{
    alphabetSize = 6
};

# else //ifndef USE_EXTRA_CHARACTER_Z

static const char alphabet[] = "$ACGNTZ";
static const int whichPileInverse[] =
{
    '$', 'A', 'C', 'G', 'N', 'T', 'Z'
};
enum
{
    alphabetSize = 7
};

# endif //ifndef USE_EXTRA_CHARACTER_Z

#else //ifdef USE_STANDARD_LEXICOGRAPHIC_ORDER

# ifndef USE_EXTRA_CHARACTER_Z

static const char alphabet[] = "$ACGTN";
static const int whichPileInverse[] =
{
    '$', 'A', 'C', 'G', 'T', 'N'
};
enum
{
    alphabetSize = 6
};

# else //ifndef USE_EXTRA_CHARACTER_Z

static const char alphabet[] = "$ACGTNZ";
static const int whichPileInverse[] =
{
    '$', 'A', 'C', 'G', 'T', 'N', 'Z'
};
enum
{
    alphabetSize = 7
};

# endif //ifndef USE_EXTRA_CHARACTER_Z
#endif //ifdef USE_STANDARD_LEXICOGRAPHIC_ORDER


// Type able to contain symbols of the alphabet (in biologic case 6 ($,A,C,G,N,T))
typedef uchar AlphabetSymbol;

// Next is a character that should not be in the alphabet
// because it gets used as a marker in a couple of places in the code
static const char notInAlphabet( 'Z' );

// One of the characters has special status: marks end of a string
static const char terminatorChar( '$' );

// One of the characters has special status: 'don't know' character
static const char dontKnowChar( 'N' );

const int nv( -1 );


#ifdef USE_STANDARD_LEXICOGRAPHIC_ORDER
// encoding is ($,A,C,G,N,T)=(0,1,2,3,4,5)
//static const char baseNames [] = "ACGNT";


static const int whichPile[] =
{
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, 0, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, /* next is 'A' */
    1, nv, 2, nv, nv, nv, 3, nv, nv, nv, nv, nv, nv, 4, nv, /* next is 'P' */
    nv, nv, nv, nv, 5, nv, nv, nv, nv, nv,
# ifdef USE_EXTRA_CHARACTER_Z
    6, // Z for Huffman
# else //ifdef USE_EXTRA_CHARACTER_Z
    nv,
# endif //ifdef USE_EXTRA_CHARACTER_Z
    nv, nv, nv, nv, nv,
    nv, /* next is 'a' */
    1, nv, 2, nv, nv, nv, 3, nv, nv, nv, nv, nv, nv, 4, nv, /* next is 'p' */
    nv, nv, nv, nv, 5, nv, nv, nv, nv, nv,
# ifdef USE_EXTRA_CHARACTER_Z
    6, // Z for Huffman
# else //ifdef USE_EXTRA_CHARACTER_Z
    nv,
# endif //ifdef USE_EXTRA_CHARACTER_Z
    nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv
};
#else
// encoding is ($,A,C,G,T,N)=(0,1,2,3,4,5) + optional Z=6
//static const char baseNames [] = "ACGTN";

static const int whichPile[] =
{
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, 0, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, /* next is 'A' */
    1, nv, 2, nv, nv, nv, 3, nv, nv, nv, nv, nv, nv, 5, nv, /* next is 'P' */
    nv, nv, nv, nv, 4, nv, nv, nv, nv, nv,
# ifdef USE_EXTRA_CHARACTER_Z
    6, // Z for Huffman
# else //ifdef USE_EXTRA_CHARACTER_Z
    nv,
# endif //ifdef USE_EXTRA_CHARACTER_Z
    nv, nv, nv, nv, nv,
    nv, /* next is 'a' */
    1, nv, 2, nv, nv, nv, 3, nv, nv, nv, nv, nv, nv, 5, nv, /* next is 'p' */
    nv, nv, nv, nv, 4, nv, nv, nv, nv, nv,
# ifdef USE_EXTRA_CHARACTER_Z
    6, // Z for Huffman
# else //ifdef USE_EXTRA_CHARACTER_Z
    nv,
# endif //ifdef USE_EXTRA_CHARACTER_Z
    nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv,
    nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv, nv
};
#endif


// Used by BackTracker classes
typedef bool AlphabetFlag[alphabetSize];


#endif

