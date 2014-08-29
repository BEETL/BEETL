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

//#define DONT_USE_MMAP


class BwtWriterBase;

class BwtReaderBase
{
public:
    BwtReaderBase( const string &filename );
    BwtReaderBase( const BwtReaderBase &obj );
    virtual ~BwtReaderBase();
    virtual BwtReaderBase *clone() const = 0;

    virtual LetterNumber readAndCount( LetterCount &c, const LetterNumber numChars ) = 0;
    LetterNumber readAndCount( LetterCount &c );
    virtual LetterNumber readAndSend( BwtWriterBase &writer, const LetterNumber numChars ) = 0;
    virtual LetterNumber readAndSend( BwtWriterBase &writer );
    virtual LetterNumber operator()( char *p, LetterNumber numChars ) = 0;
    virtual void rewindFile( void ) = 0;
    virtual LetterNumber tellg( void ) const = 0;
    virtual int seek( const LetterNumber posInFile, const LetterNumber baseNumber ) = 0;

    const string filename_;
protected:
    FILE *pFile_;
    vector<uchar> buf_;
}; // ~class BwtReaderBase

class BwtReaderASCII : public BwtReaderBase
{
public:
    BwtReaderASCII( const string &filename ) :
        BwtReaderBase( filename ),
        currentPos_( 0 ),
        lastChar_( notInAlphabet ),
        runLength_( 0 )
    {
    }

    BwtReaderASCII( const BwtReaderASCII &obj ) :
        BwtReaderBase( obj ),
        currentPos_( obj.currentPos_ ),
        lastChar_( obj.lastChar_ ),
        runLength_( obj.runLength_ )
    {
    }

    virtual ~BwtReaderASCII() {}
    virtual BwtReaderASCII *clone() const
    {
        return new BwtReaderASCII( *this );
    };

    virtual LetterNumber readAndCount( LetterCount &c, const LetterNumber numChars );

    virtual LetterNumber readAndSend( BwtWriterBase &writer, const LetterNumber numChars );

    virtual LetterNumber operator()( char *p, LetterNumber numChars );

    virtual void rewindFile( void );
    virtual LetterNumber tellg( void ) const;
    virtual int seek( const LetterNumber posInFile, const LetterNumber baseNumber );

protected:
    LetterNumber currentPos_;
    uchar lastChar_;
    uint runLength_;
}; // ~class BwtReaderASCII


class BwtReaderRunLengthBase : public BwtReaderBase
{
public:
    BwtReaderRunLengthBase( const string &filename );
    BwtReaderRunLengthBase( const BwtReaderRunLengthBase &obj );

    virtual ~BwtReaderRunLengthBase() {}

    virtual LetterNumber readAndCount( LetterCount &c, const LetterNumber numChars );

    virtual LetterNumber readAndSend( BwtWriterBase &writer, const LetterNumber numChars );

    virtual LetterNumber operator()( char *p, LetterNumber numChars );

    virtual void rewindFile( void );
    virtual LetterNumber tellg( void ) const;
    virtual int seek( const LetterNumber posInFile, const LetterNumber baseNumber );

    virtual bool getRun( void ) = 0;

protected:
    vector<uint> lengths_;
    vector<uchar> codes_;
    uchar *pBuf_;
    uchar *pBufMax_;
    bool finished_;

public: // exposed for buildIndex. TODO: make private again
    uchar lastChar_;
    uint runLength_;
    LetterNumber currentPos_;
    LetterNumber currentPosInFile_;
};

class BwtReaderRunLength : public BwtReaderRunLengthBase
{
public:
    BwtReaderRunLength( const string &filename );
    BwtReaderRunLength( const BwtReaderRunLength &obj );

    virtual ~BwtReaderRunLength() {}
    virtual BwtReaderRunLength *clone() const
    {
        return new BwtReaderRunLength( *this );
    };

    virtual bool getRun( void );
};

const vector<char> rleV3Header = { 'B', 'W', 'T', 13, 10, 26, 3, 0 };
class BwtReaderRunLengthV3 : public BwtReaderRunLengthBase
{
public:
    BwtReaderRunLengthV3( const string &filename );
    BwtReaderRunLengthV3( const BwtReaderRunLengthV3 &obj );

    virtual ~BwtReaderRunLengthV3() {}
    virtual BwtReaderRunLengthV3 *clone() const
    {
        return new BwtReaderRunLengthV3( *this );
    };

    virtual bool getRun( void );
    virtual void rewindFile( void );
    virtual LetterNumber tellg( void ) const;
    virtual int seek( const LetterNumber posInFile, const LetterNumber baseNumber );

protected:
    vector<uchar> symbolForRunLength1ForPile_;
    vector<LetterNumber> maxEncodedRunLengthForPile_;
    uchar firstContinuationSymbol_;
    LetterNumber maxEncodedRunLengthMultiplierForContinuationSymbol_;
    long firstDataByteInFile_;

    void prefetchNextByte();
    int prefetchedByte_;
};


class BwtReaderIncrementalRunLength : public BwtReaderBase
{
public:
    BwtReaderIncrementalRunLength( const string &filename );

    virtual ~BwtReaderIncrementalRunLength() {}
    virtual BwtReaderIncrementalRunLength *clone() const
    {
        return new BwtReaderIncrementalRunLength( *this );
    };

    virtual LetterNumber readAndCount( LetterCount &c, const LetterNumber numChars );

    virtual LetterNumber readAndSend( BwtWriterBase &writer, const LetterNumber numChars );

    virtual LetterNumber operator()( char *p, LetterNumber numChars );

    virtual void rewindFile( void );
    virtual LetterNumber tellg( void ) const;
    virtual int seek( const LetterNumber posInFile, const LetterNumber baseNumber ) { assert(false && "todo"); }

    bool getRun( void );
    void defragment();

protected:
    uint lengths_[256];
    uchar codes_[256];
    uint runLength_;
    uchar *pBuf_;
    uchar *pBufMax_;
    uchar lastChar_;
    uchar lastMetadata_;
    bool finished_;
    LetterNumber currentPos_;
    size_t posInRamFile_;

    int fileNum_;
    vector<size_t> posInRamFiles_;
    vector<int> stackedFileNums_;

}; // class ~BwtReaderIncrementalRunLength


#ifdef ACTIVATE_HUFFMAN

// new input module to support Huffman encoded input
// migrated & adapted from compression.cpp in Tony's Misc CVS tree
// Tobias, 28/11/11

class BwtReaderHuffman : public BwtReaderBase
{
public:
    BwtReaderHuffman( const string &filename );
    virtual BwtReaderHuffman *clone() const
    {
        assert( false && "todo" );
        return new BwtReaderHuffman( *this );
    };

    virtual LetterNumber readAndCount( LetterCount &c, const LetterNumber numChars );

    virtual LetterNumber readAndSend( BwtWriterBase &writer, const LetterNumber numChars );

    virtual LetterNumber operator()( char *p, LetterNumber numChars );

    virtual void rewindFile( void );
    virtual LetterNumber tellg( void ) const;
    virtual int seek( const LetterNumber posInFile, const LetterNumber baseNumber ) { assert(false && "todo"); }

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
    LetterNumber currentPos_; // position in the file
    bool firstRun_; // self explainatory


}; // class ~BwtReaderHuffman

#endif //ifdef ACTIVATE_HUFFMAN


class BwtReaderRunLengthRam : public BwtReaderBase
{
public:
    BwtReaderRunLengthRam( const string &filename );
    BwtReaderRunLengthRam( const BwtReaderRunLengthRam & );

    virtual ~BwtReaderRunLengthRam();
    virtual BwtReaderRunLengthRam *clone() const
    {
        return new BwtReaderRunLengthRam( *this );
    };

    virtual LetterNumber readAndCount( LetterCount &c, const LetterNumber numChars );

    virtual LetterNumber readAndSend( BwtWriterBase &writer, const LetterNumber numChars );

    virtual LetterNumber operator()( char *p, LetterNumber numChars );

    virtual void rewindFile( void );
    virtual LetterNumber tellg( void ) const;
    virtual int seek( const LetterNumber posInFile, const LetterNumber baseNumber ) { assert(false && "todo"); }

    bool getRun( void );

protected:
    uint lengths_[256];
    uchar codes_[256];
    uint runLength_;
    uchar lastChar_;
    LetterNumber currentPos_;

    char *fullFileBuf_;
    LetterNumber posInFullFileBuf_;
    LetterNumber sizeOfFullFileBuf_;

#ifndef DONT_USE_MMAP
    size_t mmapLength_;
#endif

private:
    bool isClonedObject_;
}; // class ~BwtReaderRunLengthRam


BwtReaderBase* instantiateBwtPileReader( const string &pileFilename, const string &useShm = "", const bool keepBwtInRam = false, const bool forceNotUseIndexClass = false );
vector <BwtReaderBase *> instantiateBwtPileReaders( const string &bwtPrefix, const string &useShm = "" );

#endif
