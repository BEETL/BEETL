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

#ifndef INCLUDED_INTERVALHANDLER_BASE_HH
#define INCLUDED_INTERVALHANDLER_BASE_HH

#include "Alphabet.hh"
#include "LetterCount.hh"
#include "RangeStore.hh"
#include "Types.hh"


//
// IntervalHandler
//
// Idea here is that different algorithms can be implemented by defining
// new subclasses of IntervalHandler

struct IntervalHandlerBase
{
    IntervalHandlerBase() : needSubstring( false ) {}
    virtual ~IntervalHandlerBase() {}

    void countString( char *bwtSubstring, int length, LetterCount &countsThisRange );

    bool needSubstring;
    /*
        virtual bool needSubstring()
        {
            return false;    // Determines if bwtSubstring is filled in the following callbacks
        }
    */

    virtual void foundInAOnly(
        const int pileNum,
        const LetterCount &countsSoFarA,
        const LetterCount &countsThisRangeA,
        const char *bwtSubstring, // = NULL if needSubtring() method returns false OR if it would be too large (see code for threshold)
        Range &thisRangeA,
        AlphabetFlag &propagateIntervalA,
        const int cycle
    ) = 0;

    virtual void foundInBOnly(
        const int pileNum,
        const LetterCount &countsSoFarB,
        const LetterCount &countsThisRangeB,
        const char *bwtSubstring, // = NULL if needSubtring() method returns false OR if it would be too large (see code for threshold)
        Range &thisRangeB,
        AlphabetFlag &propagateIntervalB,
        const int cycle
    ) = 0;

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
    ) = 0;

    // Override this method when using derived Range types
    virtual Range &getSubIntervalRange (
        const string &word,
        const LetterNumber pos,
        const LetterNumber num,
        const bool isBkptExtension,
        Range &parentRange,
        const int subIntervalNum
    )
    {
        //        assert( false );
        static Range r;
        r = Range( word, pos, num, isBkptExtension );
        return r;
    }

    void createOutputFile( const int subsetThreadNum, const int i, const int j, const int cycle, const string &outputDirectory );
    std::ofstream outFile_;
};

typedef void ( IntervalHandlerBase::*IntervalHandler_FoundCallbackPtr ) (
    const int pileNum,
    const LetterCount &countsSoFar,
    const LetterCount &countsThisRange,
    const char *bwtSubstring,
    Range &thisRange,
    AlphabetFlag &propagateInterval,
    const int cycle
);

#endif
