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

#ifndef INCLUDED_ONEBWTBACKTRACKER_HH
#define INCLUDED_ONEBWTBACKTRACKER_HH

#include "BackTrackerBase.hh"
#include "BwtReader.hh"
#include "EndPosFile.hh"
#include "IntervalHandlerBase.hh"
#include "LetterCount.hh"
#include "RangeStore.hh"
#include "Types.hh"
#include "libzoo/util/Logger.hh"

#include <string>

using namespace std;


// OneBwtBackTracker class implements the backward search and takes its
// cue from IntervalHandler as to whether to continue down a particular
// branch of the search tree
class OneBwtBackTracker: public BackTrackerBase
{
public:
    OneBwtBackTracker(
        BwtReaderBase *inBwt,
        LetterNumber &currentPos,
        RangeStoreExternal &r,
        LetterCount &countsSoFar,
        const string &subset,
        const int cycle,
        const bool doesPropagateBkptToSeqNumInSet,
        const bool noComparisonSkip,
        const bool propagateSequence,
        EndPosFile &endPosFile
    );

    void skipIfNecessary( const Range &thisRange,
                          LetterNumber &currentPos,
                          BwtReaderBase &inBwt,
                          LetterCount &countsSoFar );

    void process(
        int i,
        string &thisWord,
        IntervalHandlerBase &intervalHandler,
        Range &rangeDerivedObject
    );

    BwtReaderBase *inBwt_;
    LetterNumber &currentPos_;

    RangeStoreExternal &r_;
    LetterCount &countsSoFar_;

    const string &subset_;
    const int cycle_;

    //    LetterNumber numRanges_;
    //    LetterNumber numSingletonRanges_;

    AlphabetFlag propagateInterval_;

    const bool doesPropagateBkptToSeqNumInSet_;
    //    const bool noComparisonSkip_;
    EndPosFile &endPosFile_;
};

#endif
