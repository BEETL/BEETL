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

#include "WitnessReader.hh"

using namespace std;

WitnessReader::WitnessReader(
    const string &lcpFileName,
    const string &bwtFileName,
    int witnessLength,
    int minimumSupport,
    bool rleBWT
)
{
    for ( int i = 0; i < ReadBufferSize; i++ )
        lcpBuf_[i] = 0;
    at_ = 0;
    filePos_ = 0;
    lastLcpBlockSupport_ = 0;
    filledTo_ = 0;
    pFile_ = fopen( lcpFileName.c_str(), "rb" );
    if ( rleBWT )
        bwtReader_ = new BwtReaderRunLength( bwtFileName );
    else
        bwtReader_ = new BwtReaderASCII( bwtFileName );
    witnessLength_ = witnessLength;
    minimumSupport_ = minimumSupport;
    totalCountSoFar_.clear();
    refill_();
}
WitnessReader::~WitnessReader()
{
    fclose( pFile_ );
    delete bwtReader_;
}
LetterCount WitnessReader::TotalCountSoFar()
{
    return totalCountSoFar_;
}
int WitnessReader::currentWitnessCount()
{
    return lastLcpBlockSupport_;
}
int WitnessReader::currentWitnessBlockStart()
{
    return filePos_ - filledTo_ + at_ - lastLcpBlockSupport_;
}
bool WitnessReader::nextWitnessBlock( LetterCount &lc )
{
    while ( nextCandidateLcpBlock_() )
    {
        //catch the bwt file up with the lcp...
        bwtReader_->readAndCount( totalCountSoFar_, currentWitnessBlockStart() - lastBlockEnd_ );
        //get the actual individual letter counts we're interested in...
        lc.clear();
        bwtReader_->readAndCount( lc, currentWitnessCount() );

        totalCountSoFar_ += lc;

        int totalSupport = currentWitnessCount() - lc.count_[
                               whichPile[( int )'$']
                           ];

        if ( totalSupport > minimumSupport_ )
        {
            lc.count_[
                whichPile[( int )'$']
            ] = 0;
            return true;
        }
    }
    bwtReader_->readAndCount( totalCountSoFar_ );
    return false;
}

void WitnessReader::refill_()
{
    filledTo_ = fread( lcpBuf_, sizeof( int ), ReadBufferSize, pFile_ );
    filePos_ += filledTo_;
    at_ = 0;
}

bool WitnessReader::nextCandidateLcpBlock_()
{
    lastBlockEnd_ = currentWitnessCount() + currentWitnessBlockStart();
    lastLcpBlockSupport_ = 1;
    while ( filledTo_ > 0 )
    {
        while ( at_ < filledTo_ )
        {
            if ( lcpBuf_[at_] >= witnessLength_ )
                lastLcpBlockSupport_++;
            else if ( lastLcpBlockSupport_ > minimumSupport_ )
                return true;
            else
                lastLcpBlockSupport_ = 1;
            at_++;
        }
        refill_();
    }
    return false;
}

void WitnessReader::test()
{
    char bwtChars[filledTo_];
    ( *bwtReader_ )( bwtChars, filledTo_ );
    for ( int i = 0; i < filledTo_; i++ )
        cout << bwtChars[i] << "        " << lcpBuf_[i] << endl;
}
