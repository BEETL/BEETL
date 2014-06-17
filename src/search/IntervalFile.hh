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

#ifndef INTERVAL_FILE_HH
#define INTERVAL_FILE_HH

#include "Types.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using std::string;
using std::vector;


struct IntervalRecord
{
    IntervalRecord(): kmer( "" ), position( 0 ), count( 0 ) {}
    IntervalRecord(
        string inKmer,
        LetterNumber inPosition,
        LetterNumber inCount
    ):
        kmer( inKmer ),
        position( inPosition ),
        count( inCount )
    {}

    friend std::ostream &operator<<( std::ostream &os, const IntervalRecord &obj );

    static bool bwtPositionCompare( const IntervalRecord &lhs, const IntervalRecord &rhs )
    {
        return lhs.position < rhs.position;
    }

    string kmer;
    LetterNumber position;
    LetterNumber count;
    vector<LetterNumber> dollarSignPositions;
};

class IntervalWriter
{
public:
    IntervalWriter( std::ostream &file );
    void write( const IntervalRecord &ir ) const;
    void writeV2( const IntervalRecord &ir ) const;
private:
    std::ostream &file_;
};

class IntervalReader
{
public:
    IntervalReader( std::istream &file ): file_( file ) {}
    bool read( IntervalRecord &ir );
    vector<IntervalRecord> readFullFileAsVector();
private:
    std::istream &file_;
};

#endif
