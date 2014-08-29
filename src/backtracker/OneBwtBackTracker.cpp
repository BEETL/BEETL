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

#include "OneBwtBackTracker.hh"

#include "libzoo/util/Logger.hh"

using namespace std;


OneBwtBackTracker::OneBwtBackTracker(
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
)
    : BackTrackerBase( subset, cycle, noComparisonSkip, propagateSequence )
    , inBwt_( inBwt ),
    currentPos_( currentPos ),
    r_( r ),
    countsSoFar_( countsSoFar ),
    subset_( subset ),
    cycle_( cycle ),
    //    numRanges_( 0 ),
    //    numSingletonRanges_( 0 ),
    doesPropagateBkptToSeqNumInSet_( doesPropagateBkptToSeqNumInSet )
    //    noComparisonSkip_( noComparisonSkip )
    , endPosFile_( endPosFile )
{
    for ( int l( 0 ); l < alphabetSize; ++l )
        propagateInterval_[l] = false;
}

void OneBwtBackTracker::process (
    int pileNum,
    string &thisWord,
    IntervalHandlerBase &intervalHandler,
    Range &thisRange
)
{
    LetterCount countsThisRange;
    bool notAtLast( true );

    processSingletons(
        pileNum
        , notAtLast
        , r_
        , thisRange
        , currentPos_
        , inBwt_
        , countsSoFar_
        , countsThisRange
        , intervalHandler
        , propagateInterval_
        , thisWord
        , doesPropagateBkptToSeqNumInSet_
        , ( IntervalHandler_FoundCallbackPtr )( &IntervalHandlerBase::foundInAOnly )
        , endPosFile_
        , 1
    );

} // ~OneBwtBackTracker::operator()


