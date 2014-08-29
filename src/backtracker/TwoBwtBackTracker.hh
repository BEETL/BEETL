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

#ifndef TWO_BWT_BACKTRACKER_HH
#define TWO_BWT_BACKTRACKER_HH

#include "BackTrackerBase.hh"


class TwoBwtBackTracker : public BackTrackerBase
{
public:
    TwoBwtBackTracker( BwtReaderBase *inBwtA, BwtReaderBase *inBwtB,
                       LetterNumber &currentPosA, LetterNumber &currentPosB,
                       RangeStoreExternal &rA, RangeStoreExternal &rB,
                       LetterCount &countsSoFarA, LetterCount &countsSoFarB,
                       int minOcc, const int maxLength, const string &subset, const int cycle,
                       const bool doesPropagateBkptToSeqNumInSetA,
                       const bool doesPropagateBkptToSeqNumInSetB,
                       const bool noComparisonSkip,
                       const bool propagateSequence );


    void process (
        int pileNum
        , string &thisWord
        , IntervalHandlerBase &intervalHandler_
    );

    BwtReaderBase *inBwtA_;
    BwtReaderBase *inBwtB_;

    LetterNumber &currentPosA_;
    LetterNumber &currentPosB_;

    RangeStoreExternal &rA_;
    RangeStoreExternal &rB_;
    LetterCount &countsSoFarA_;
    LetterCount &countsSoFarB_;

    int minOcc_;
    int maxLength_;
    const string &subset_;
    const int cycle_;

    //    LetterNumber numRanges_;
    //    LetterNumber numSingletonRanges_;
    LetterNumber numNotSkippedA_;
    LetterNumber numNotSkippedB_;
    LetterNumber numSkippedA_;
    LetterNumber numSkippedB_;

    AlphabetFlag propagateIntervalA_;
    AlphabetFlag propagateIntervalB_;

    const bool doesPropagateBkptToSeqNumInSetA_;
    const bool doesPropagateBkptToSeqNumInSetB_;
    //    const bool noComparisonSkip_;
};

#endif // TWO_BWT_BACKTRACKER_HH
