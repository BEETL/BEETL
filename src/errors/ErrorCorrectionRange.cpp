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

#include "ErrorCorrectionRange.hh"

using namespace std;


bool ErrorCorrectionRange::writeTo( TemporaryFile *pFile, RangeState &currentState ) const
{
    Range::writeTo( pFile, currentState );

    writeCompressedNum( pFile, static_cast<LetterNumber>( data_.errorIntervalType ) );
    writeCompressedNum( pFile, data_.correctionForBwtPosns.size() );
    for ( uint i = 0; i < data_.correctionForBwtPosns.size(); i++ )
        writeCompressedNum( pFile, data_.correctionForBwtPosns[i] );
    writeCompressedNum( pFile, data_.errorsForBwtPosns.size() );
    for ( uint i = 0; i < data_.errorsForBwtPosns.size(); i++ )
        writeCompressedNum( pFile, data_.errorsForBwtPosns[i] );

    return true;
}

bool ErrorCorrectionRange::readFrom( TemporaryFile *pFile, RangeState &currentState )
{
    if ( Range::readFrom( pFile, currentState ) == false )
        return false;

    LetterNumber intervalTypeNo;
    readCompressedNum( pFile, intervalTypeNo );
    data_.errorIntervalType = static_cast<IntervalType>( intervalTypeNo );

    assert( data_.correctionForBwtPosns.empty() );
    LetterNumber numBwtPosns;
    readCompressedNum( pFile, numBwtPosns );
    for ( uint i = 0; i < numBwtPosns; i++ )
    {
        LetterNumber newBwtPos;
        readCompressedNum( pFile, newBwtPos );
        data_.correctionForBwtPosns.push_back( newBwtPos );
    }

    assert( data_.errorsForBwtPosns.empty() );
    LetterNumber numErrBwtPosns;
    readCompressedNum( pFile, numErrBwtPosns );
    for ( uint i = 0; i < numErrBwtPosns; i++ )
    {
        LetterNumber newErrBwtPos;
        readCompressedNum( pFile, newErrBwtPos );
        data_.errorsForBwtPosns.push_back( newErrBwtPos );
    }

    return true;
}
