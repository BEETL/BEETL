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

#include "SequenceExtractor.hh"

#include <fstream>

using namespace std;


SequenceExtractor::SequenceExtractor() : isActive_( false ) {}

void SequenceExtractor::init( const string &seqNumFilename )
{
    seqNums_.clear();
    index_ = 0;
    currentSeqNum_ = -1;

    ifstream seqNumFile( seqNumFilename );
    string line;
    while ( getline( seqNumFile, line ) )
    {
        istringstream iss( line );
        SequenceNumber seqNum = ( SequenceNumber ) - 1;
        iss >> seqNum;
        if ( seqNum != ( SequenceNumber ) - 1 )
            seqNums_.push_back( seqNum );
    }

    std::sort( seqNums_.begin(), seqNums_.end() );
    isActive_ = true;
}

bool SequenceExtractor::doWeExtractNextSequence()
{
    if ( !isActive_ ) return true;

    ++currentSeqNum_;
    if ( index_ >= seqNums_.size() )
        return false;
    if ( seqNums_[index_] != currentSeqNum_ )
        return false;

    // Handle multiple occurrences of the same seqnum in the input file
    while ( seqNums_[++index_] == currentSeqNum_ ) {}

    return true;
}
