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

#ifndef INCLUDED_BACKTRACKER_BASE_HH
#define INCLUDED_BACKTRACKER_BASE_HH

#include "BwtReader.hh"
#include "EndPosFile.hh"
#include "IntervalHandlerBase.hh"
#include "LetterCount.hh"
#include "RangeStore.hh"
#include "Types.hh"
#include "libzoo/util/Logger.hh"

#include <string>

using namespace std;


// BackTrackerBase class implements the backward search and takes its
// cue from IntervalHandler as to whether to continue down a particular
// branch of the search tree
struct BackTrackerBase
{
    BackTrackerBase( const string &subset, const int cycle, const bool noComparisonSkip, const bool propagateSequence );

    void skipIfNecessary( const Range &thisRange,
                          LetterNumber &currentPos,
                          BwtReaderBase &inBwt,
                          LetterCount &countsSoFar );

    void operator() ( int i, string &thisWord, IntervalHandlerBase &intervalHandler );

    const string &subset_;
    const int cycle_;
    const bool noComparisonSkip_;

    LetterNumber numRanges_;
    LetterNumber numSingletonRanges_;

    void processSingletons(
        const int pileNum
        , bool &notAtLast
        , RangeStoreExternal &rA_
        , Range &thisRange
        , LetterNumber &currentPosA_
        , BwtReaderBase *inBwtA_
        , LetterCount &countsSoFar
        , LetterCount &countsThisRange
        , IntervalHandlerBase &intervalHandler
        , AlphabetFlag &propagateInterval
        , string &thisWord
        , const bool doesPropagateToEnd
        , IntervalHandler_FoundCallbackPtr foundCallbackPtr
        , EndPosFile &endPosFile
        , int sampleId
    );

    void prepareCallbackArgs(
        Range &thisRange
        , LetterNumber &currentPos
        , BwtReaderBase *inBwt
        , LetterCount &countsSoFar
        , LetterCount &countsThisRange
        , IntervalHandlerBase &intervalHandler
        , char *&bwtSubstring
    );

    void updatePropagatedSuffixWord( bool &hasChild, const Range &thisRange, string &thisWord, const AlphabetSymbol l );

    vector<char> bwtSubstringStore_;
    const bool propagateSequence_;
};

#endif

