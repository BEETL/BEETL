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

#ifndef EXTENDER_INTERVAL_HANDLER_HH
#define EXTENDER_INTERVAL_HANDLER_HH

#include "EndPosFile.hh"
#include "IntervalHandlerBase.hh"
#include "RangeStore.hh"


struct ExtenderIntervalHandler : public IntervalHandlerBase
{
    ExtenderIntervalHandler( EndPosFile &endPosFile );

    virtual ~ExtenderIntervalHandler() {}

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

    virtual void foundInAOnly
    ( const int pileNum,
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
        static Range r;
        r = Range( word, pos, num, isBkptExtension, true, parentRange.userData_ );
        return r;
    }

private:
    EndPosFile &endPosFile_;
};

#endif // EXTENDER_INTERVAL_HANDLER_HH
