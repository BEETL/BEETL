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

#ifndef KMER_SEARCH_RANGE_HH
#define KMER_SEARCH_RANGE_HH

#include "Config.hh"
#include "Range.hh"

#include <string>


struct KmerSearchDataForInterval//: public DataForSubIntervals
{
    KmerSearchDataForInterval()
        : start( 0 ), end( 0 )
    {}

    //    virtual ~BwtCorrectorDataForInterval() {}
    //    virtual void clear() {}

    int start, end;
};


class KmerSearchRange : public Range
{
public:
    KmerSearchRange(
        const string &word,
        const LetterNumber pos,
        const LetterNumber num,
        const bool isBkptExtension = false
    ):
        Range( word, pos, num, isBkptExtension )
    {}

    KmerSearchRange(
        const string &word,
        const LetterNumber pos,
        const LetterNumber num,
        const bool isBkptExtension,
        int start, int end
    ):
        Range( word, pos, num, isBkptExtension )
    {
        data_.start = start;
        data_.end = end;
    }

    KmerSearchRange(
        const string &word,
        const LetterNumber pos,
        const LetterNumber num,
        const bool isBkptExtension,
        Range &parentRange,
        const int subIntervalNum
    ):
        Range( word, pos, num, isBkptExtension ),
        data_( dynamic_cast< KmerSearchRange & >( parentRange ).getDataForSubInterval( subIntervalNum ) )
    {
        assert( data_.start < data_.end );
    }

    KmerSearchRange( void ) :
        Range()
    {}

    virtual ~KmerSearchRange() {}

    virtual void clear()
    {
        *this = KmerSearchRange();
    }

    virtual bool writeTo( TemporaryFile *pFile, RangeState &currentState ) const;
    virtual bool readFrom( TemporaryFile *pFile, RangeState &currentState );

    KmerSearchDataForInterval data_;

    void clearDataForSubIntervals()
    {
        dataForSubIntervals_.resize( alphabetSize );
        for ( int i = 0; i < alphabetSize; i++ )
        {
            dataForSubIntervals_[i].start = 0;
            dataForSubIntervals_[i].end = 0;
        }
    }
    KmerSearchDataForInterval &getDataForSubInterval( int l )
    {
        return dataForSubIntervals_[l];
    }

private:
    vector<KmerSearchDataForInterval> dataForSubIntervals_;
};


#endif // KMER_SEARCH_RANGE_HH
