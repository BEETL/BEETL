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

#ifndef DEFINED_ALPHABET_HH
#define DEFINED_ALPHABET_HH

#include "Config.hh"


// convention is that first char is the string terminator char and
// is lexicographically less than the rest

#ifdef USE_STANDARD_LEXICOGRAPHIC_ORDER 
static const char alphabet[]="$ACGNTZ";
#else
static const char alphabet[]="$ACGTN";
#endif

//static const int alphabetSize(strlen(alphabet));
enum
{
  alphabetSize=7
};

// Next is a character that should not be in the alphabet
// because it gets used as a marker in a couple of places in the code
static const char notInAlphabet('Z'); 

// One of the characters has special status: marks end of a string
static const char terminatorChar('$');

// One of the characters has special status: 'don't know' character
static const char dontKnowChar('N');

const int nv(-1);


#ifdef USE_STANDARD_LEXICOGRAPHIC_ORDER 
// encoding is ($,A,C,G,N,T)=(0,1,2,3,4,5)
//static const char baseNames [] = "ACGNT";


static const int whichPile[] = 
{
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,0, nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv, /* next is 'A' */ 
1,nv,2,nv,nv,nv,3,nv,nv,nv,nv,nv,nv,4,nv, /* next is 'P' */
nv,nv,nv,nv,5,nv,nv,nv,nv,nv,6,nv,nv,nv,nv,nv, // added Z for Huffman
nv, /* next is 'a' */ 
1,nv,2,nv,nv,nv,3,nv,nv,nv,nv,nv,nv,4,nv, /* next is 'p' */
nv,nv,nv,nv,5,nv,nv,nv,nv,nv,6,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv
};
#else
// encoding is ($,A,C,G,T,N)=(0,1,2,3,4,5)
//static const char baseNames [] = "ACGTN";

static const int whichPile[] = 
{
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,0, nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv, /* next is 'A' */ 
1,nv,2,nv,nv,nv,3,nv,nv,nv,nv,nv,nv,5,nv, /* next is 'P' */
nv,nv,nv,nv,4,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv, /* next is 'a' */ 
1,nv,2,nv,nv,nv,3,nv,nv,nv,nv,nv,nv,5,nv, /* next is 'p' */
nv,nv,nv,nv,4,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv
};
#endif

#endif

