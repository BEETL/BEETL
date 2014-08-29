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

#include "KmerSearchRange.hh"

using namespace std;


bool KmerSearchRange::writeTo( TemporaryFile *pFile, RangeState &currentState ) const
{
    Range::writeTo( pFile, currentState );

    writeCompressedNum( pFile, static_cast<LetterNumber>( data_.start ) );
    writeCompressedNum( pFile, static_cast<LetterNumber>( data_.end ) );

    return true;
}

bool KmerSearchRange::readFrom( TemporaryFile *pFile, RangeState &currentState )
{
    if ( Range::readFrom( pFile, currentState ) == false )
        return false;

    LetterNumber start;
    readCompressedNum( pFile, start );
    data_.start = start;

    LetterNumber end;
    readCompressedNum( pFile, end );
    data_.end = end;

    return true;
}
