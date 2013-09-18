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

#ifndef INCLUDED_ERROR_CORRECTION_RANGE_HH
#define INCLUDED_ERROR_CORRECTION_RANGE_HH

#include "Config.hh"
#include "countWords/Range.hh"

#include <string>


enum IntervalType
{
    INTERVAL_TYPE_DEFAULT = 0,
    INTERVAL_TYPE_CORRECTOR = 1,
    INTERVAL_TYPE_ERROR = 2
};

class ErrorCorrectionRange : public Range
{
public:
    ErrorCorrectionRange(
        const string &word,
        const LetterNumber pos,
        const LetterNumber num,
        const bool isBkptExtension,
        const IntervalType intervalType,
        const vector<LetterNumber> &bwtPosns,
        const vector<LetterNumber> &errBwtPosns
    ) :
        Range( word, pos, num, isBkptExtension ),
        intervalType_( intervalType ),
        bwtPosns_( bwtPosns ),
        errBwtPosns_( errBwtPosns )
    {}

    ErrorCorrectionRange(
        const string &word,
        const LetterNumber pos,
        const LetterNumber num,
        const bool isBkptExtension,
        const IntervalType intervalType
    ) :
        Range( word, pos, num, isBkptExtension ),
        intervalType_( intervalType ),
        bwtPosns_( vector<LetterNumber>() ),
        errBwtPosns_( vector<LetterNumber>() )
    {}

    ErrorCorrectionRange(
        const string &word,
        const LetterNumber pos,
        const LetterNumber num,
        const bool isBkptExtension = false
    ):
        Range( word, pos, num, isBkptExtension ),
        intervalType_( INTERVAL_TYPE_DEFAULT ),
        bwtPosns_( vector<LetterNumber>() ),
        errBwtPosns_( vector<LetterNumber>() )
    {}

    ErrorCorrectionRange( void ) :
        Range(),
        intervalType_( INTERVAL_TYPE_DEFAULT ),
        bwtPosns_( vector<LetterNumber>() ),
        errBwtPosns_( vector<LetterNumber>() )
    {}

    virtual ~ErrorCorrectionRange() {}

    virtual void clear()
    {
        *this = ErrorCorrectionRange();
    }

    virtual bool writeTo( TemporaryFile *pFile, RangeState &currentState ) const;
    virtual bool readFrom( TemporaryFile *pFile, RangeState &currentState );

    IntervalType intervalType_;
    vector<LetterNumber> bwtPosns_;
    vector<LetterNumber> errBwtPosns_;

};


#endif
