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

#ifndef INCLUDED_BWTCORRECTORINTERVALHANDLER_HH
#define INCLUDED_BWTCORRECTORINTERVALHANDLER_HH

#include "Config.hh"
#include "ErrorInfo.hh"
#include "ErrorCorrectionRange.hh"
#include "IntervalHandlerBase.hh"
#include "RangeStore.hh"

#include <iostream>
#include <map>

using namespace std;


struct BwtCorrectorIntervalHandler : public IntervalHandlerBase
{
    BwtCorrectorIntervalHandler(
        ErrorStore &inErrorStore,
        int minWitnessLength,
        int minOccurrences,
        int cycle
    ):
        errorStore_( inErrorStore ),
        minWitnessLength_( minWitnessLength ),
        minOccurrences_( minOccurrences )
    {
        needSubstring = ( cycle >= minWitnessLength_ );
    }

    virtual ~BwtCorrectorIntervalHandler() {}
    virtual void foundInBoth
    ( const int pileNum,
      const LetterCount &countsThisRangeA,
      const LetterCount &countsThisRangeB,
      const Range &thisRangeA,
      const Range &thisRangeB,
      AlphabetFlag &propagateIntervalA,
      AlphabetFlag &propagateIntervalB,
      bool &isBreakpointDetected,
      const int cycle
    )
    {
        assert( false );
    }

    /*
        virtual bool needSubstring()
        {
            return needSubstring_;
    //        return true;    // Determines which of the following 2 callbacks to use
        }
    */
    virtual void foundInAOnly(
        const int pileNum,
        const LetterCount &countsSoFarA,
        const LetterCount &countsThisRangeA,
        const char *bwtSubstring,
        Range &thisRangeA,
        AlphabetFlag &propagateIntervalA,
        const int cycle
    );

    virtual void foundInBOnly
    ( const int pileNum,
      const LetterCount &countsSoFarB,
      const LetterCount &countsThisRangeB,
      const char *bwtSubstring,
      Range &thisRangeB,
      AlphabetFlag &propagateIntervalB,
      const int cycle
    )
    {
        assert( false );
    }

    virtual Range &getSubIntervalRange (
        const string &word,
        const LetterNumber pos,
        const LetterNumber num,
        const bool isBkptExtension,
        Range &parentRange,
        const int subIntervalNum
    )
    {
        static ErrorCorrectionRange r;
        r = ErrorCorrectionRange( word, pos, num, isBkptExtension, parentRange, subIntervalNum );
        return r;
    }

    //errorStore_ - a reference to the ErrorStore owned by the BwtCorrector algorithm, so that in processing this interval, we can
    //add/update error objects for BWT positions (bwtposn's) should we come across new/existing bwtposn's which are likely errors
    //in this interval
    ErrorStore &errorStore_;

    //the smallest length of substring Q for which the Q-interval should be inspected for errors.
    int minWitnessLength_;

    //the minimum number of times a letter has to occur in an interval for it to be considered the 'correct' letter
    unsigned int minOccurrences_;

    // which cycle we are at
    //    const bool needSubstring_;

    //method to see if we can discover errors in an interval
    //  TRUE if we can - i.e. there is exactly one letter with at least T occurrences, and there is at least one occurrence of some other
    //      letter which is not a dollar.
    //  FALSE if we can't, and hence we don't then scan along the bwtSubstring to find the locations of the putative errors within the range.
    bool defaultDetermineErrors( LetterCount intervalLetterCount, int &correct );
};

#endif
