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

#ifndef INCLUDED_BWTWRITER_HH
#define INCLUDED_BWTWRITER_HH

#include "Alphabet.hh"
#include "Config.hh"
#include "LetterCount.hh"
#include "Tools.hh"
#include "Types.hh"

#include <cstdio>
#include <map>
#include <string>


struct BwtWriterBase
{
    virtual ~BwtWriterBase() {};

    virtual void operator()( const char *p, LetterNumber numChars ) = 0;

    virtual void sendRun( char c, LetterNumber runLength ) = 0;
    virtual void sendRunOfPreExistingData( char c, LetterNumber runLength, int fileNum, size_t posInRamFile, LetterNumber remainingRunLength )
    {
        sendRun( c, runLength );
    };
    virtual char getLastChar()
    {
        assert( false && "virtual method needs implementing" );
    }
    virtual bool isIncremental()
    {
        return false;
    }
    virtual void flush()
    {
        assert( false && "virtual method needs implementing" );
    }
}; // ~BwtWriterBase

struct BwtWriterFile : public BwtWriterBase
{
    BwtWriterFile( const string &fileName );

    virtual ~BwtWriterFile();

    virtual void flush();


    FILE *pFile_;
}; // ~BwtWriterBase



struct BwtWriterASCII : public BwtWriterFile
{
    BwtWriterASCII( const string &fileName );

    virtual ~BwtWriterASCII();

    virtual void operator()( const char *p, LetterNumber numChars );

    virtual void sendRun( char c, LetterNumber runLength );
    virtual char getLastChar();

    char lastChar_;
};


struct BwtWriterRunLengthBase : public BwtWriterFile
{
    BwtWriterRunLengthBase( const string &fileName, const int baseFieldWidthInBits ):
        BwtWriterFile( fileName ),
        runLength_( 0 ), pBuf_( buf_ ), pBufMax_( buf_ + ReadBufferSize ), lastChar_( notInAlphabet )
#ifdef REPORT_COMPRESSION_RATIO
        , charsReceived_( 0 ), bytesWritten_( 0 )
#endif
        , baseFieldWidthInBits_( baseFieldWidthInBits )
        , lengthFieldWidthInBits_( 8 - baseFieldWidthInBits )
        , lengthFieldMask_( ~( ( uchar )0 ) << baseFieldWidthInBits )
    {
        assert( baseFieldWidthInBits_ >= 3 );
        assert( baseFieldWidthInBits_ <= 7 );
    }

    virtual ~BwtWriterRunLengthBase();

    virtual void operator()( const char *p, LetterNumber numChars );

    void sendChar( char c );

    virtual void encodeRun( char c, LetterNumber runLength );

    virtual void sendRun( char c, LetterNumber runLength );
    virtual char getLastChar();
    virtual void flush();

    LetterNumber runLength_;
    uchar buf_[ReadBufferSize];
    uchar *pBuf_;
    const uchar *pBufMax_;
    uchar lastChar_;
#ifdef REPORT_COMPRESSION_RATIO
    LetterNumber charsReceived_;
    LetterNumber bytesWritten_;
#endif
    const int baseFieldWidthInBits_;
    const int lengthFieldWidthInBits_;
    const uchar lengthFieldMask_;

protected:
    virtual void flushBuffer();

    //#define GENERATE_RLE_HISTOGRAM
#ifdef GENERATE_RLE_HISTOGRAM
    std::map< std::pair< uchar, LetterNumber >, LetterNumber > histogram_;
#endif
}; // ~BwtWriterRunLengthBase


struct BwtWriterRunLength : public BwtWriterRunLengthBase
{
    BwtWriterRunLength( const string &fileName ):
        BwtWriterRunLengthBase( fileName, 4 )
    {}

};
typedef BwtWriterRunLength BwtWriterRunLength_4_4;


struct BwtWriterRunLength_5_3 : public BwtWriterRunLengthBase
{
    BwtWriterRunLength_5_3( const string &fileName ):
        BwtWriterRunLengthBase( fileName, 3 )
    {}
};


// New Bwt encoding.
// First version has its encoding calculated offline, so not really dynamic
// A,C,G,T: 1-63
// N: 1
// $: 1, 2, 3
struct BwtWriterRunLengthDynamic : public BwtWriterRunLengthBase
{
    BwtWriterRunLengthDynamic( const string &fileName );

    virtual void encodeRun( char c, LetterNumber runLength );

protected:
    vector<uchar> symbolForRunLength1ForPile_;
    vector<LetterNumber> maxEncodedRunLengthForPile_;
};


struct BwtWriterIncrementalRunLength : public BwtWriterFile
{
    BwtWriterIncrementalRunLength( const string &fileName );

    virtual ~BwtWriterIncrementalRunLength();

    virtual void operator()( const char *p, LetterNumber numChars );

    void terminateLastInsertion();
    void sendChar( unsigned char c, unsigned char metadata );

    void encodeRun( char c, LetterNumber runLength );

    virtual void sendRun( char c, LetterNumber runLength );
    virtual void sendRunOfPreExistingData( char c, LetterNumber runLength, int fileNum, size_t posInRamFile, LetterNumber remainingRunLength );
    virtual bool isIncremental()
    {
        return true;
    }

    LetterNumber runLength_;
    uchar buf_[ReadBufferSize];
    uchar *pBuf_;
    const uchar *pBufMax_;
    uchar lastChar_;
#ifdef REPORT_COMPRESSION_RATIO
    LetterNumber charsReceived_;
    LetterNumber bytesWritten_;
#endif

private:
    uint fileNum_;
    uint fileNumInReader_;
    size_t filePosInReader_;
    LetterNumber remainingRunLengthInReader_;
    bool lastFileReturnNeeded_;
    unsigned char onHoldUntilNextReturn_letter_;
    unsigned char onHoldUntilNextReturn_runLength_;
    unsigned char onHoldUntilNextReturn_metadata_;
}; // ~BwtWriterIncrementalRunLength

// new output module to support huffman encoded output
// migrated & adapted from compression.cpp in Tony's Misc CVS tree
// Tobias, 28/11/11

struct BwtWriterHuffman : public BwtWriterFile
{

    BwtWriterHuffman( const string &fileName ) :
        BwtWriterFile( fileName ),
        bitsUsed_( 0 ), lastChar_( notInAlphabet ), runLength_( 0 ), huffmanBufferPos( 0 )
#ifdef REPORT_COMPRESSION_RATIO
        , charsReceived_( 0 ), bytesWritten_( 0 )
#endif
    {
        soFar_.ull = 0;
        numBuf_.ull = 0;

        for ( int i = 0; i < huffmanWriterBufferSize; i++ )
        {
            symBuf[i] = 0;
        }
    }
    void sendToken( unsigned long long code, LetterNumber length );

    virtual void operator()( const char *p, LetterNumber numChars );

    virtual void sendRun( char c, LetterNumber runLength );

    void sendNum( LetterNumber runLength );

    void emptyBuffer( void );

    void processBuffer( int itemsToPrint );

    virtual ~BwtWriterHuffman(); // ~BwtWriterHuffman()

    BitBuffer soFar_;
    BitBuffer toAdd_;
    BitBuffer numBuf_;
    int bitsUsed_;
    char lastChar_;
    LetterNumber runLength_;
    uint huffmanBufferPos;
    uchar symBuf[huffmanWriterBufferSize];

#ifdef REPORT_COMPRESSION_RATIO
    LetterNumber charsReceived_;
    LetterNumber bytesWritten_;
#endif
}; // ~BwtWriterHuffman

struct BwtWriterImplicit : public BwtWriterBase
{
    BwtWriterImplicit( BwtWriterBase *pWriter ) :
        BwtWriterBase(),
        firstSAP_( 0 ),
        pWriter_( pWriter ),
        lastChar_( notInAlphabet ), lastRun_( 0 ), inSAP_( false ) {}

    virtual ~BwtWriterImplicit();

    virtual void operator()( const char *p, LetterNumber numChars );

    virtual void sendRun( char c, LetterNumber runLength );

    void flushSAP( void );

    LetterCount countSAP_;
    int firstSAP_;

    BwtWriterBase *pWriter_;

    char lastChar_;
    LetterNumber lastRun_;


    bool inSAP_;
    //  FILE* pFile_;
}; // ~BwtWriterImplicit : public BwtWriterBase




#endif

