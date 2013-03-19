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

#include "Alphabet.hh"
#include "Config.hh"
#include "LetterCount.hh"
#include "Types.hh"

#include <cassert>
#include <cstdio>
#include <string>
#include <vector>


class BwtWriterBase;

class BwtReaderBase
{
public:
    BwtReaderBase( const string &fileName );
    virtual ~BwtReaderBase();

    virtual unsigned int readAndCount( LetterCount &c, const LetterCountType numChars ) = 0;
    uint readAndCount( LetterCount &c );
    virtual uint readAndSend( BwtWriterBase &writer, const int numChars ) = 0;
    virtual uint readAndSend( BwtWriterBase &writer );
    virtual int operator()( char *p, int numChars ) = 0;
    virtual void rewindFile( void ) = 0;
    virtual LetterCountType tellg( void ) const = 0;

protected:
    FILE *pFile_;

    char buf_[ReadBufferSize];
}; // ~class BwtReaderBase

class BwtReaderASCII : public BwtReaderBase
{
public:
    BwtReaderASCII( const string &fileName ) :
        BwtReaderBase( fileName ),
        currentPos_( 0 ),
        lastChar_( notInAlphabet ),
        runLength_( 0 )
    {
    }

    virtual ~BwtReaderASCII() {}

    virtual unsigned int readAndCount( LetterCount &c, const LetterCountType numChars );

    virtual uint readAndSend( BwtWriterBase &writer, const int numChars );

    virtual int operator()( char *p, int numChars );

    virtual void rewindFile( void );

    virtual LetterCountType tellg( void ) const;

protected:
    LetterCountType currentPos_;
    uchar lastChar_;
    uint runLength_;
}; // ~class BwtReaderASCII


class BwtReaderRunLength : public BwtReaderBase
{
public:
    BwtReaderRunLength( const string &fileName );

    virtual ~BwtReaderRunLength() {}

    virtual uint readAndCount( LetterCount &c, const LetterCountType numChars );

    virtual uint readAndSend( BwtWriterBase &writer, const int numChars );

    virtual int operator()( char *p, int numChars );

    virtual void rewindFile( void );

    virtual LetterCountType tellg( void ) const;

    bool getRun( void );

protected:
    uint lengths_[256];
    uchar codes_[256];
    uchar buf_[ReadBufferSize];
    uint runLength_;
    uchar *pBuf_;
    uchar *pBufMax_;
    uchar lastChar_;
    bool finished_;
    LetterCountType currentPos_;

}; // class ~BwtReaderRunLength

class BwtReaderRunLengthIndex : public BwtReaderRunLength
{
public:
    BwtReaderRunLengthIndex( const string &fileName );

    virtual ~BwtReaderRunLengthIndex() {}

    virtual uint readAndCount( LetterCount &c, const LetterCountType numChars );

    virtual uint readAndSend( BwtWriterBase &writer, const int numChars )
    {
        assert( 1 == 0 );
    }

    virtual int operator()( char *p, int numChars )
    {
        assert( 1 == 0 );
    }

    virtual void rewindFile( void );

    //  virtual LetterCountType tellg( void ) const;

    void buildIndex( FILE *pFile, const int indexBinSize );

    void initIndex( const LetterCount &current );

    //  bool getRun(void);

protected:


    //  string fileName_;
    string indexFileName_;

    LetterCount temp_;
    LetterCount currentIndex_; // count at last index at or before current pos
    LetterCount current_; // index right now
    LetterCount next_; // counts for next index point
    LetterCountType currentIndexPos_; // position in BWT at last index point at or before current pos
    LetterCountType currentFilePos_; // position in file at last index point at or before current pos
    LetterCountType nextPos_; // position in BWT for next index point
    LetterCountType nextFilePos_;// position in file for next index point
    bool isNextIndex_; // is there a next index point?

    FILE *pIndexFile_;

}; // class ~BwtReaderRunLengthIndex

class BwtReaderIncrementalRunLength : public BwtReaderBase
{
public:
    BwtReaderIncrementalRunLength( const string &fileName );

    virtual ~BwtReaderIncrementalRunLength() {}

    virtual uint readAndCount( LetterCount &c, const LetterCountType numChars );

    virtual uint readAndSend( BwtWriterBase &writer, const int numChars );

    virtual int operator()( char *p, int numChars );

    virtual void rewindFile( void );

    virtual LetterCountType tellg( void ) const;


    bool getRun( void );
    void defragment();

protected:
    uint lengths_[256];
    uchar codes_[256];
    uchar buf_[ReadBufferSize];
    uint runLength_;
    uchar *pBuf_;
    uchar *pBufMax_;
    uchar lastChar_;
    uchar lastMetadata_;
    bool finished_;
    LetterCountType currentPos_;
    size_t posInRamFile_;

    int fileNum_;
    vector<size_t> posInRamFiles_;
    vector<int> stackedFileNums_;

}; // class ~BwtReaderIncrementalRunLength

// new input module to support Huffman encoded input
// migrated & adapted from compression.cpp in Tony's Misc CVS tree
// Tobias, 28/11/11

class BwtReaderHuffman : public BwtReaderBase
{
public:
    BwtReaderHuffman( const string &fileName );

    virtual uint readAndCount( LetterCount &c, const LetterCountType numChars );

    virtual uint readAndSend( BwtWriterBase &writer, const int numChars );

    virtual int operator()( char *p, int numChars );

    virtual void rewindFile( void );

    virtual LetterCountType tellg( void ) const;

    uint getNum( int &i );

    bool getRun( void );

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
