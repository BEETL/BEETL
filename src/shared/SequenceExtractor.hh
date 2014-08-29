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

using std::vector;


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
