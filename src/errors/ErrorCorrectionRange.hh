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

#ifndef INCLUDED_ERROR_CORRECTION_RANGE_HH
#define INCLUDED_ERROR_CORRECTION_RANGE_HH

#include "Config.hh"
#include "Range.hh"

#include <string>


enum IntervalType
{
    INTERVAL_TYPE_DEFAULT = 0,
    INTERVAL_TYPE_CORRECTOR = 1,
    INTERVAL_TYPE_ERROR = 2
};


struct BwtCorrectorDataForInterval//: public DataForSubIntervals
{
    BwtCorrectorDataForInterval()
        : errorIntervalType( INTERVAL_TYPE_DEFAULT )
    {}

    //    virtual ~BwtCorrectorDataForInterval() {}
    //    virtual void clear() {}

    IntervalType errorIntervalType;
    vector<LetterNumber> correctionForBwtPosns; //correctionBwtPosns;
    vector<LetterNumber> errorsForBwtPosns; //errBwtPosns;

};


class ErrorCorrectionRange : public Range
{
public:
    ErrorCorrectionRange(
        const string &word,
        const LetterNumber pos,
        const LetterNumber num,
        const bool isBkptExtension,
        Range &parentRange,
        const int subIntervalNum
    ) :
        Range( word, pos, num, isBkptExtension ),
        data_( dynamic_cast< ErrorCorrectionRange & >( parentRange ).getDataForSubInterval( subIntervalNum ) )
    {}

    ErrorCorrectionRange(
        const string &word,
        const LetterNumber pos,
        const LetterNumber num,
        const bool isBkptExtension = false
    ):
        Range( word, pos, num, isBkptExtension )
    {}

    ErrorCorrectionRange( void ) :
        Range()
    {}

    virtual ~ErrorCorrectionRange() {}

    virtual void clear()
    {
        *this = ErrorCorrectionRange();
    }

    virtual bool writeTo( TemporaryFile *pFile, RangeState &currentState ) const;
    virtual bool readFrom( TemporaryFile *pFile, RangeState &currentState );

    BwtCorrectorDataForInterval data_;

    void clearDataForSubIntervals()
    {
        dataForSubIntervals_.resize( alphabetSize );
        for ( int i = 0; i < alphabetSize; i++ )
        {
            dataForSubIntervals_[i].errorIntervalType = INTERVAL_TYPE_DEFAULT;
            dataForSubIntervals_[i].correctionForBwtPosns = vector<LetterNumber>();
            dataForSubIntervals_[i].errorsForBwtPosns = vector<LetterNumber>();
        }
    }
    BwtCorrectorDataForInterval &getDataForSubInterval( int l )
    {
        return dataForSubIntervals_[l];
    }

private:
    vector<BwtCorrectorDataForInterval> dataForSubIntervals_;
};


#endif
