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
 **/

#include "ReadBuffer.hh"

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>

using namespace std;


bool ReadBufferBase::getNext
( SequenceNumber &seqNum, LetterNumber &seqPtr )
{
    if ( thisEntry_ == maxEntry_ )
    {
        if ( lastIter_ == true ) return false;

        int numRead( read( fdSeq_, seqBufBase_, blockSize_ ) );

#ifdef DEBUG
        std::cout << "RB: read " << numRead << " bytes, asked for " << blockSize_ << std::endl;
#endif

        if ( numRead == 0 ) return false;
        else if ( numRead < blockSize_ )
        {
            // ensure whole number of seqs are read, else a problem
            if ( ( numRead % readSize_ ) != 0 )
            {
                cerr << "Not all sequences could be read. Aborting." << endl;
                exit( EXIT_FAILURE );
            }
            lastIter_ = true;
            maxEntry_ = numRead / readSize_;
        }

#ifdef TRACK_SEQUENCE_NUMBER
        assert( read( fdNum_, seqNum_, sizeof( SequenceNumber )*maxEntry_ ) == sizeof( SequenceNumber )*maxEntry_ );
#ifdef DEBUG
        std::cout << "RB: read " << numRead << " seqnums" << endl;
#endif
#endif
        assert( read( fdPtr_, seqPtr_, sizeof( LetterNumber )*maxEntry_ ) == ( int )sizeof( LetterNumber )*maxEntry_ );

#ifdef DEBUG
        std::cout << "RB: read " << numRead << " seqptrs" << std::endl;
#endif
        thisEntry_ = 0;
    } // ~if
    seqBuf_ = seqBufBase_;
    seqBuf_ += ( readSize_ * thisEntry_ );
    //    seqBuf=seqBuf_;
#ifdef TRACK_SEQUENCE_NUMBER
    seqNum = seqNum_[thisEntry_];
#endif

    seqPtr = seqPtr_[thisEntry_];
    thisEntry_++;
    return true;
} // ~getNext
