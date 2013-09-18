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

#ifndef INCLUDED_HITECERROR_HH
#define INCLUDED_HITECERROR_HH

using namespace std;

#include <iostream>
#include "Config.hh"
#include "errors/ErrorTypes.hh"
#include "countWords/RangeStore.hh"

struct HiTECError : public Range
{
    HiTECError(
        const LetterCountType pos,
        const LetterCountType num,
        char correctTo,
        int witnessLength,
        vector<ReadErrorLocation> &inErrorsGathered
    ):
        witnessLength_( witnessLength ),
        errorsGathered_( inErrorsGathered )
    {}

    virtual ~HiTECError() {}
    void process(
        const int pileNum,
        const LetterCount &countsSoFar,
        const LetterCount &countsThisRange,
        AlphabetFlag &propagateInterval
    );
    int witnessLength_;
    vector<ReadErrorLocation> &errorsGathered_;
};

#endif
