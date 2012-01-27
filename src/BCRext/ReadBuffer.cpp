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

#include "ReadBuffer.hh"

#include <iostream>
#include <fstream>
#include <string>
//#include <cassert>
#include <cstring>


bool ReadBufferBase::getNext
( SequenceNumberType& seqNum, LetterCountType& seqPtr )    
{
    if (thisEntry_==maxEntry_)
    {
      if (lastIter_==true) return false;

      int numRead(read( fdSeq_, seqBufBase_, blockSize_));

#ifdef DEBUG
      std::cout << "RB: read " << numRead << " bytes, asked for " << blockSize_<< std::endl;
#endif
      
      if (numRead==0) return false;
      else if (numRead<blockSize_) 
      {
        // ensure whole number of seqs are read, else a problem
        if ((numRead % readSize_) != 0) {
            cerr << "Not all sequences could be read. Aborting." << endl;
            exit(-1);
        }
        lastIter_=true; 
	maxEntry_=numRead/readSize_; 
      }

#ifdef TRACK_SEQUENCE_NUMBER
      assert(read( fdNum_, seqNum_, sizeof(SequenceNumberType)*maxEntry_)==sizeof(SequenceNumberType)*maxEntry_);
#ifdef DEBUG
      std::cout << "RB: read " << numRead << " seqnums" <<endl;
#endif
#endif
      assert(read( fdPtr_, seqPtr_, sizeof(LetterCountType)*maxEntry_)==(int)sizeof(LetterCountType)*maxEntry_);

#ifdef DEBUG
      std::cout << "RB: read " << numRead << " seqptrs" <<std::endl;
#endif
      thisEntry_=0;
    } // ~if
    seqBuf_=seqBufBase_;
    seqBuf_+=(readSize_*thisEntry_);  
    //    seqBuf=seqBuf_;
#ifdef TRACK_SEQUENCE_NUMBER
    seqNum=seqNum_[thisEntry_];  
#endif

    seqPtr=seqPtr_[thisEntry_];  
    thisEntry_++;
    return true;
  } // ~getNext
