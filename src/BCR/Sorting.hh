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

#ifndef SORTED_INCLUDED
#define SORTED_INCLUDED

#include "Tools.hh"

#include <vector>

using std::vector;


struct sortElement
{
#if BUILD_LCP == 0

    sortElement() {}

    sortElement( dataTypedimAlpha z, dataTypeNChar x, dataTypeNSeq y )
        : pileN( z )
        , posN( x )
        , seqN( y )
    {}

    dataTypelenSeq getLcpCurN() const
    {
        return 0;
    }
    dataTypelenSeq getLcpSucN() const
    {
        return 0;
    }
    void setLcpCurN( const dataTypelenSeq val ) { }
    void setLcpSucN( const dataTypelenSeq val ) { }

#else

    sortElement() : lcpCurN( 0 ), lcpSucN( 0 ) {}

    sortElement( dataTypedimAlpha z, dataTypeNChar x, dataTypeNSeq y, dataTypelenSeq l1 = 0, dataTypelenSeq l2 = 0 )
        : pileN( z )
        , posN( x )
        , seqN( y )
        , lcpCurN( l1 )
        , lcpSucN( l2 )
    {}

    dataTypelenSeq getLcpCurN() const
    {
        return lcpCurN;
    }
    dataTypelenSeq getLcpSucN() const
    {
        return lcpSucN;
    }
    void setLcpCurN( const dataTypelenSeq val )
    {
        lcpCurN = val;
    }
    void setLcpSucN( const dataTypelenSeq val )
    {
        lcpSucN = val;
    }

#endif
    ~sortElement() {};
    dataTypedimAlpha pileN;
    dataTypeNChar posN;
    dataTypeNSeq seqN;
#if BUILD_LCP == 1
    dataTypelenSeq lcpCurN;
    dataTypelenSeq lcpSucN;
#endif

#if USE_ATTRIBUTE_PACKED == 1
} __attribute__ ( ( packed ) );
#else
};
#endif

void quickSort( vector< sortElement > &v );


#endif
