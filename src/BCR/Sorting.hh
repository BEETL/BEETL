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

#ifndef SORTED_INCLUDED
#define SORTED_INCLUDED

#include "Alphabet.hh"
#include "Types.hh"

//2020-12-04
#include "Tools.hh"  // it is included for BUILD_LCP

#include <vector>

using std::vector;

struct sortElement
{
#if BUILD_LCP == 0

    sortElement() {}

    sortElement( AlphabetSymbol z, LetterNumber x, SequenceNumber y )
        : pileN( z )
        , posN( x )
        , seqN( y )
    {}

    SequenceLength getLcpCurN() const
    {
        return 0;
    }
    SequenceLength getLcpSucN() const
    {
        return 0;
    }
    void setLcpCurN( const SequenceLength val ) { }
    void setLcpSucN( const SequenceLength val ) { }

#else

    sortElement() : pileN( 0 ), lcpCurN( 0 ), lcpSucN( 0 ) {}

    sortElement( AlphabetSymbol z, LetterNumber x, SequenceNumber y, SequenceLength l1 = 0, SequenceLength l2 = 0 )
        : pileN( z )
        , posN( x )
        , seqN( y )
        , lcpCurN( l1 )
        , lcpSucN( l2 )
    {}

    SequenceLength getLcpCurN() const
    {
        return lcpCurN;
    }
    SequenceLength getLcpSucN() const
    {
        return lcpSucN;
    }
    void setLcpCurN( const SequenceLength val )
    {
        lcpCurN = val;
    }
    void setLcpSucN( const SequenceLength val )
    {
        lcpSucN = val;
    }

#endif
    ~sortElement() {};
    AlphabetSymbol pileN;
    LetterNumber posN;
    SequenceNumber seqN;
#if BUILD_LCP == 1
    SequenceLength lcpCurN;
    SequenceLength lcpSucN;
#endif

#if USE_ATTRIBUTE_PACKED == 1
} __attribute__ ( ( packed ) );
#else
};
#endif

void quickSort( vector< sortElement > &v );


#endif
