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

#ifndef INCLUDED_ERRORINTERVALHANDLER_HH
#define INCLUDED_ERRORINTERVALHANDLER_HH

using namespace std;

#include <iostream>
#include "Config.hh"
#include "errors/HiTECErrorLocation.hh"
#include "countWords/IntervalHandlerBase.hh"
#include <sstream>

//#define PROPAGATE_PREFIX 1

struct ErrorIntervalHandler : public IntervalHandlerBase
{
    ErrorIntervalHandler(
        int witnessLength,
        int backtrackerIterationNum,
        vector<HiTECErrorLocation> &inErrorsGathered
    ):
        witnessLength_( witnessLength ),
        iterationNumber_( backtrackerIterationNum ),
        errorsGathered_( inErrorsGathered )
    {}

    virtual ~ErrorIntervalHandler() {}
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

    virtual void foundInBOnly
    ( const int pileNum,
      const LetterCount &countsSoFarB,
      const LetterCount &countsThisRangeB,
      const Range &thisRangeB,
      AlphabetFlag &propagateIntervalB );

    int witnessLength_;
    int iterationNumber_;
    vector<HiTECErrorLocation> &errorsGathered_;
};

#endif
