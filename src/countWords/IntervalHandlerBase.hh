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

#ifndef INCLUDED_INTERVALHANDLER_BASE_HH
#define INCLUDED_INTERVALHANDLER_BASE_HH

#include "Alphabet.hh"
#include "LetterCount.hh"
#include "RangeStore.hh"
#include "Types.hh"


typedef bool AlphabetFlag[alphabetSize];

//
// IntervalHandler
//
// Idea here is that different algorithms can be implemented by defining
// new subclasses of IntervalHandler


struct IntervalHandlerBase
{
    virtual ~IntervalHandlerBase() {}

    virtual void foundInBoth
    ( const int pileNum,
      const LetterCount &countsThisRangeA, const LetterCount &countsThisRangeB,
      const Range &thisRangeA, const Range &thisRangeB,
      AlphabetFlag &propagateIntervalA, AlphabetFlag &propagateIntervalB )=0;

    virtual void foundInAOnly
    ( const int pileNum,
      const LetterCount &countsSoFarA,
      const LetterCount &countsThisRangeA,
      const Range &thisRangeA,
      AlphabetFlag &propagateIntervalA )=0;

    virtual void foundInBOnly
    ( const int pileNum,
      const LetterCount &countsSoFarB,
      const LetterCount &countsThisRangeB,
      const Range &thisRangeB,
      AlphabetFlag &propagateIntervalB )=0;
};
#endif
