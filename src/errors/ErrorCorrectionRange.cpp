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

#include "ErrorCorrectionRange.hh"

using namespace std;


bool ErrorCorrectionRange::writeTo( TemporaryFile *pFile, RangeState &currentState ) const
{
    Range::writeTo( pFile, currentState );

    writeCompressedNum( pFile, static_cast<LetterNumber>( intervalType_ ) );
    writeCompressedNum( pFile, bwtPosns_.size() );
    for ( int i = 0; i < bwtPosns_.size(); i++ )
        writeCompressedNum( pFile, bwtPosns_[i] );
    writeCompressedNum( pFile, errBwtPosns_.size() );
    for ( int i = 0; i < errBwtPosns_.size(); i++ )
        writeCompressedNum( pFile, errBwtPosns_[i] );

    return true;
}

bool ErrorCorrectionRange::readFrom( TemporaryFile *pFile, RangeState &currentState )
{
    if ( Range::readFrom( pFile, currentState ) == false )
        return false;

    LetterNumber intervalTypeNo;
    readCompressedNum( pFile, intervalTypeNo );
    intervalType_ = static_cast<IntervalType>( intervalTypeNo );

    assert( bwtPosns_.empty() );
    LetterNumber numBwtPosns;
    readCompressedNum( pFile, numBwtPosns );
    for ( int i = 0; i < numBwtPosns; i++ )
    {
        LetterNumber newBwtPos;
        readCompressedNum( pFile, newBwtPos );
        bwtPosns_.push_back( newBwtPos );
    }

    assert( errBwtPosns_.empty() );
    LetterNumber numErrBwtPosns;
    readCompressedNum( pFile, numErrBwtPosns );
    for ( int i = 0; i < numErrBwtPosns; i++ )
    {
        LetterNumber newErrBwtPos;
        readCompressedNum( pFile, newErrBwtPos );
        errBwtPosns_.push_back( newErrBwtPos );
    }

    return true;
}
