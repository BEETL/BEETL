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

#ifndef INCLUDED_BWTCORRECTORINTERVALHANDLER_HH
#define INCLUDED_BWTCORRECTORINTERVALHANDLER_HH

using namespace std;

#include <iostream>
#include <map>
#include "Config.hh"
#include "ErrorInfo.hh"
#include "ErrorCorrectionRange.hh"
#include "countWords/IntervalHandlerBase.hh"
#include "countWords/RangeStore.hh"

struct BwtCorrectorIntervalHandler : public IntervalHandlerBase
{
    BwtCorrectorIntervalHandler(
        ErrorStore &inErrorStore,
        int minWitnessLength,
        int minOccurrences
    ):
        errorStore_( inErrorStore ),
        minWitnessLength_( minWitnessLength ),
        minOccurrences_( minOccurrences )
    {}

    virtual ~BwtCorrectorIntervalHandler() {}
    virtual void foundInBoth
    ( const int pileNum,
      const LetterCount &countsThisRangeA, const LetterCount &countsThisRangeB,
      const Range &thisRangeA, const Range &thisRangeB,
      AlphabetFlag &propagateIntervalA, AlphabetFlag &propagateIntervalB,
      bool &isBreakpointDetected );

    virtual void foundInAOnly
    ( const int pileNum,
      const LetterCount &countsSoFarA,
      const LetterCount &countsThisRangeA,
      const Range &thisRangeA,
      AlphabetFlag &propagateIntervalA );

    virtual void foundInAOnly(
        const int pileNum,
        const LetterCount &countsSoFarA,
        //the actual characters of the bwt corresponding to this range
        char *bwtSubstring,
        LetterCount &countsThisRangeA,
        const Range &thisRangeA,
        AlphabetFlag &propagateIntervalA,
        //a vector of interval types - either CORRECTOR, ERROR or DEFAULT for the descending intervals
        //(i.e. the vector is the same length as the alphabet, one element per letter/possible extension)
        IntervalType ( &errorIntervalType )[alphabetSize],
        //a vector of lists of BWT positions with which descending correction intervals should be tagged
        //  this object is passed to OneBwtBacktracker, which creates the descending intervals
        vector<LetterNumber> ( &correctionForBwtPosns )[alphabetSize],
        //a vector of lists of BWT positions with which descending error intervals should be tagged
        vector<LetterNumber> ( &errorsForBwtPosns )[alphabetSize],
        const int cycle
    );

    virtual void foundInBOnly
    ( const int pileNum,
      const LetterCount &countsSoFarB,
      const LetterCount &countsThisRangeB,
      const Range &thisRangeB,
      AlphabetFlag &propagateIntervalB );

    //errorStore_ - a reference to the ErrorStore owned by the BwtCorrector algorithm, so that in processing this interval, we can
    //add/update error objects for BWT positions (bwtposn's) should we come across new/existing bwtposn's which are likely errors
    //in this interval
    ErrorStore &errorStore_;

    //the smallest length of substring Q for which the Q-interval should be inspected for errors.
    int minWitnessLength_;

    //the minimum number of times a letter has to occur in an interval for it to be considered the 'correct' letter
    int minOccurrences_;

    //method to see if we can discover errors in an interval
    //  TRUE if we can - i.e. there is exactly one letter with at least T occurrences, and there is at least one occurrence of some other
    //      letter which is not a dollar.
    //  FALSE if we can't, and hence we don't then scan along the bwtSubstring to find the locations of the putative errors within the range.
    bool defaultDetermineErrors( LetterCount intervalLetterCount, int &correct );

};

#endif
