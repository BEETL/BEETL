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
 **
 **
 **/

#ifndef SEQUENCE_EXTRACTOR_HH
#define SEQUENCE_EXTRACTOR_HH

#include "Types.hh"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;


class SequenceExtractor
{
public:
    SequenceExtractor();
    void init( const string &seqNumFilename );
    bool doWeExtractNextSequence();

private:
    bool isActive_;
    vector<SequenceNumber> seqNums_;
    SequenceNumber index_;
    SequenceNumber currentSeqNum_;
};

#endif // SEQUENCE_EXTRACTOR_HH
