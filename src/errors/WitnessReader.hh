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
 **l
 ** Citation: Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
 ** Lightweight BWT Construction for Very Large String Collections.
 ** Proceedings of CPM 2011, pp.219-231
 **
 **/

#ifndef INCLUDED_WITNESSREADER_HH
#define INCLUDED_WITNESSREADER_HH

#include "BwtReader.hh"
#include "IntervalHandlerBase.hh"
#include "LetterCount.hh"
#include "RangeStore.hh"
#include "Types.hh"
#include "libzoo/util/Logger.hh"

#include <string>

using namespace std;


class WitnessReader
{
public:
    WitnessReader(
        const string &lcpFileName,
        const string &bwtFileName,
        int witnessLength,
        int minimumSupport,
        bool rleBWT
    );
    virtual ~WitnessReader();
    int currentWitnessCount() const;
    LetterCount TotalCountSoFar();
    int currentWitnessBlockStart() const;
    LetterCount currentWitnessSupport();
    bool nextWitnessBlock( LetterCount &lc );
    void test();
private:
    FILE *pFile_;
    BwtReaderBase *bwtReader_;
    int lcpBuf_[ReadBufferSize];
    int filledTo_;
    int at_;
    int lastBlockEnd_;
    int filePos_;
    int witnessLength_;
    int minimumSupport_;
    int lastLcpBlockSupport_;
    LetterCount totalCountSoFar_;
    void refill_();
    bool nextCandidateLcpBlock_();
};

#endif
