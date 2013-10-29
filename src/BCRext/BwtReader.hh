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

//#define DONT_USE_MMAP


class BwtWriterBase;

class BwtReaderBase
{
public:
    BwtReaderBase( const string &fileName );
    virtual ~BwtReaderBase();
    virtual BwtReaderBase *clone() const = 0;

    virtual LetterNumber readAndCount( LetterCount &c, const LetterNumber numChars ) = 0;
    LetterNumber readAndCount( LetterCount &c );
    virtual LetterNumber readAndSend( BwtWriterBase &writer, const LetterNumber numChars ) = 0;
    virtual LetterNumber readAndSend( BwtWriterBase &writer );
    virtual LetterNumber operator()( char *p, LetterNumber numChars ) = 0;
    virtual void rewindFile( void ) = 0;
    virtual LetterNumber tellg( void ) const = 0;

protected:
    FILE *pFile_;
    const string fileName_;

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

    BwtReaderASCII( const BwtReaderASCII &obj ) :
        BwtReaderBase( obj.fileName_ ),
        currentPos_( 0 ),
        lastChar_( notInAlphabet ),
        runLength_( 0 )
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

protected:
    LetterNumber currentPos_;
    uchar lastChar_;
    uint runLength_;
}; // ~class BwtReaderASCII


class BwtReaderRunLength : public BwtReaderBase
{
public:
    BwtReaderRunLength( const string &fileName );
    BwtReaderRunLength( const BwtReaderRunLength &obj );

    virtual ~BwtReaderRunLength() {}
    virtual BwtReaderRunLength *clone() const
    {
        return new BwtReaderRunLength( *this );
    };

    virtual LetterNumber readAndCount( LetterCount &c, const LetterNumber numChars );

    virtual LetterNumber readAndSend( BwtWriterBase &writer, const LetterNumber numChars );

    virtual LetterNumber operator()( char *p, LetterNumber numChars );

    virtual void rewindFile( void );

    virtual LetterNumber tellg( void ) const;

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
    LetterNumber currentPos_;

}; // class ~BwtReaderRunLength

class BwtReaderRunLengthIndex : public BwtReaderRunLength
{
public:
    BwtReaderRunLengthIndex( const string &fileName );
    //    BwtReaderRunLengthIndex( const BwtReaderRunLengthIndex & );

    BwtReaderRunLengthIndex( const BwtReaderRunLengthIndex &obj ) :
        BwtReaderRunLength( obj ),
        indexFileName_( obj.indexFileName_ ),
        pIndexFile_( obj.pIndexFile_ ),
        indexPosBwt_( obj.indexPosBwt_ ),
        indexPosFile_( obj.indexPosFile_ ),
        indexCount_( obj.indexCount_ ),
        indexNext_( obj.indexNext_ )
    {
        assert( pIndexFile_ == NULL ); // If it's not NULL, we may try to fclose it multiple times
    }


    virtual ~BwtReaderRunLengthIndex() {}
    virtual BwtReaderRunLengthIndex *clone() const
    {
        return new BwtReaderRunLengthIndex( *this );
    };

    virtual LetterNumber readAndCount( LetterCount &c, const LetterNumber numChars );

    virtual LetterNumber readAndSend( BwtWriterBase &writer, const LetterNumber numChars )
    {
        assert( 1 == 0 );
    }

    virtual LetterNumber operator()( char *p, LetterNumber numChars )
    {
        return BwtReaderRunLength::operator()( p, numChars );
        //        assert( 1 == 0 );
    }

    virtual void rewindFile( void );

    //  virtual LetterNumber tellg( void ) const;

    void buildIndex( FILE *pFile, const int indexBinSize );

    void initIndex( void );

    //  bool getRun(void);
protected:


    string indexFileName_;


    FILE *pIndexFile_;

    vector<LetterNumber> indexPosBwt_;
    vector<LetterNumber> indexPosFile_;
    vector<LetterCountCompact> indexCount_;
    uint32_t indexNext_;

}; // class ~BwtReaderRunLengthIndex

class BwtReaderIncrementalRunLength : public BwtReaderBase
{
public:
    BwtReaderIncrementalRunLength( const string &fileName );

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
    LetterNumber currentPos_;
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



class BwtReaderRunLengthRam : public BwtReaderBase
{
public:
    BwtReaderRunLengthRam( const string &fileName );
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

#endif
