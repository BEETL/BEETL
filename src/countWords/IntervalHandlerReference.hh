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

#ifndef INCLUDED_INTERVALHANDLER_REFERENCE_HH
#define INCLUDED_INTERVALHANDLER_REFERENCE_HH

#include "Config.hh"
#include "IntervalHandlerBase.hh"


//#define PROPAGATE_PREFIX 1

//
// IntervalHandler
//
// Idea here is that different algorithms can be implemented by defining
// new subclasses of IntervalHandler

struct IntervalHandlerReference : public IntervalHandlerBase
{
    IntervalHandlerReference( unsigned int minOcc ) : minOcc_( minOcc ) {}
    virtual ~IntervalHandlerReference() {}
    virtual void foundInBoth
    ( const int pileNum,
      const LetterCount &countsThisRangeA, const LetterCount &countsThisRangeB,
      const Range &thisRangeA, const Range &thisRangeB,
      AlphabetFlag &propagateIntervalA, AlphabetFlag &propagateIntervalB );

    virtual void foundInAOnly
    ( const int pileNum,
      const LetterCount &countsSoFarA,
      const LetterCount &countsThisRangeA,
      const Range &thisRangeA,
      AlphabetFlag &propagateIntervalA );

    virtual void foundInBOnly
    ( const int pileNum,
      const LetterCount &countsSoFarB,
      const LetterCount &countsThisRangeB,
      const Range &thisRangeB,
      AlphabetFlag &propagateIntervalB );

    const LetterCountType minOcc_;
};

#endif
