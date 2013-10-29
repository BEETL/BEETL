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
