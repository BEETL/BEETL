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



#ifndef INCLUDED_BWTREADER_HH
#define INCLUDED_BWTREADER_HH

#include <string>
#include <cassert>
#include <cstdio>

#include "Alphabet.hh"
#include "Config.hh"
#include "Types.hh"

using namespace std;


class BwtWriterBase;
class LetterCount;



class BwtReaderBase
{
public:
  BwtReaderBase( const std::string& fileName );
  virtual ~BwtReaderBase();

  virtual unsigned int readAndCount( LetterCount& c, const LetterCountType numChars )=0; 
  virtual uint readAndCount( LetterCount& c );
  virtual uint readAndSend( BwtWriterBase& writer, const int numChars )=0;
  virtual uint readAndSend( BwtWriterBase& writer );
  virtual int operator()( char* p, int numChars )=0;
  virtual void rewindFile(void)=0;
  virtual LetterCountType tellg( void ) const=0;

protected:
  FILE* pFile_;

  char buf_[ReadBufferSize];
}; // ~class BwtReaderBase

class BwtReaderASCII : public BwtReaderBase
{
public:
  BwtReaderASCII( const std::string& fileName ) : 
    BwtReaderBase(fileName),
    currentPos_(0),
    lastChar_(notInAlphabet),
    runLength_(0)
  {
  }

  virtual ~BwtReaderASCII() {}

  virtual unsigned int readAndCount( LetterCount& c, const LetterCountType numChars ); 

  virtual uint readAndSend( BwtWriterBase& writer, const int numChars );

  virtual int operator()( char* p, int numChars );

  virtual void rewindFile(void);

  virtual LetterCountType tellg( void ) const; 

  protected:
  LetterCountType currentPos_;
  uchar lastChar_;
  uint runLength_;
}; // ~class BwtReaderASCII

class BwtReaderRunLength : public BwtReaderBase
{
public:
  BwtReaderRunLength( const std::string& fileName );

  virtual ~BwtReaderRunLength() {}

  virtual uint readAndCount( LetterCount& c, const LetterCountType numChars ); 

  virtual uint readAndSend( BwtWriterBase& writer, const int numChars );

  virtual int operator()( char* p, int numChars );

  virtual void rewindFile(void);

  virtual LetterCountType tellg( void ) const; 


  bool getRun(void);

protected:
  uint lengths_[256];
  uchar codes_[256];
  uchar buf_[ReadBufferSize];
  uint runLength_;
  uchar* pBuf_;
  uchar* pBufMax_;
  uchar lastChar_;
  bool finished_;
  LetterCountType currentPos_;

}; // class ~BwtReaderRunLength

// new input module to support Huffman encoded input
// migrated & adapted from compression.cpp in Tony's Misc CVS tree
// Tobias, 28/11/11

class BwtReaderHuffman : public BwtReaderBase
{
public:
  BwtReaderHuffman( const std::string& fileName );

  virtual uint readAndCount( LetterCount& c, const LetterCountType numChars ); 

  virtual uint readAndSend( BwtWriterBase& writer, const int numChars );

  virtual int operator()( char* p, int numChars );

  virtual void rewindFile(void);

  virtual LetterCountType tellg( void ) const; 

  uint getNum( int& i ); 

  bool getRun(void);

protected:
    BitBuffer soFar_;
    BitBuffer toAdd_;
    
    uint runLength_; // current run length
    uchar lastChar_; // last char read until now
    int bitsUsed_; // number of bits currently used
    int numInts_; // # of ints in the file
    bool finished_; // self explainatory
    bool nearlyFinished_; // self explainatory    
    TokenTable tokenTable_; // holds shortcuts for runlength decoding
    int intCounter_; // how many ints already have been processed from BWT file
    int numSymbols_; // current number of character already decoded and waiting for output
    int maxSymbols_; //max number of character already decoded and waiting for output
    int queueCounter_; // position in the queue of decoded sumbols
    uchar symBuf[huffmanBufferSize]; // extracted character from compressed BWT
    uint runBuf[huffmanBufferSize]; // runlength of character at same pos in upper array  
    LetterCountType currentPos_; // position in the file
    bool firstRun_; // self explainatory    


}; // class ~BwtReaderHuffman


#endif
