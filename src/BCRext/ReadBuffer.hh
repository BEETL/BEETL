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

#ifndef INCLUDED_READBUFFER_HH
#define INCLUDED_READBUFFER_HH

#include "Alphabet.hh"
#include "Config.hh"
#include "Types.hh"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>


struct ReadBufferBase
{
    ReadBufferBase( int readSize, int fdSeq, int fdNum, int fdPtr ) :
        thisEntry_( ReadBufferSize ),
        maxEntry_( ReadBufferSize ),
        readSize_( readSize ),
        blockSize_( readSize_*ReadBufferSize ),
        fdSeq_( fdSeq ),
        fdNum_( fdNum ),
        fdPtr_( fdPtr ),
        lastIter_( false )
    {
#ifdef DEBUG
        std::cout << "Creating abstract read buffer: "
                  << readSize_ << " bytes per base, "
                  << blockSize_ << " bytes per block" << std::endl;
#endif
        seqBufBase_ = new char[ blockSize_ ];
        seqBuf_ = seqBufBase_;
    }
    virtual ~ReadBufferBase()
    {
        delete [] seqBufBase_;
    }


    bool getNext(
        //char*& seqBuf,
        SequenceNumber &seqNum, LetterNumber &seqPtr );

    virtual void sendTo( FILE *pFile )
    {
        sendTo( pFile, readSize_ );
        //    assert(fwrite(seqBuf_, sizeof(char), readSize_, pFile)==readSize_);
    } // ~sendTo
    void sendTo( FILE *pFile, int readSize )
    {
        size_t bytesWritten = fwrite( seqBuf_, sizeof( char ), readSize, pFile );
        if ( bytesWritten != ( size_t )readSize )
        {
            std::cerr << "Unable to write " << ReadBufferSize
                      << " chars. Aborting." << std::endl;
            exit( -1 );
        }
    }


    virtual int operator[]( const int i ) = 0;
    virtual void convertFromASCII( void ) = 0;


    int thisEntry_;
    int maxEntry_;
    const int readSize_;
    const int blockSize_;
    const int fdSeq_;
    const int fdNum_;
    const int fdPtr_;

    //  FILE* inSeq_;
    // FILE* inNum_;
    // FILE* inPtr_;

    bool lastIter_;

#ifdef TRACK_SEQUENCE_NUMBER
    SequenceNumber seqNum_[ReadBufferSize];
#endif
    LetterNumber seqPtr_[ReadBufferSize];
    // start of storage buffer
    char *seqBufBase_;
    // start of storage buffer for last entry read
    char *seqBuf_;

};



struct ReadBufferASCII : public ReadBufferBase
{
    ReadBufferASCII( int seqSize, int fdSeq, int fdNum, int fdPtr ) :
        ReadBufferBase( 1 + seqSize, fdSeq, fdNum, fdPtr )
    {
#ifdef DEBUG
        std::cout << "Creating ASCII read buffer: " << seqSize << " symbols per read"
                  << std::endl;
#endif
    } // ~ctor
    virtual int operator[]( const int i )
    {
        return whichPile[( int )seqBuf_[i]];
    } // ~operator[]
    virtual void convertFromASCII( void )
    {} // ~convertFromASCII( void )
}; // ~struct ReadBufferASCII

struct ReadBuffer4Bits : public ReadBufferBase
{
    ReadBuffer4Bits( int seqSize, int fdSeq, int fdNum, int fdPtr ) :
        ReadBufferBase( convertBasesToBytes( seqSize ), fdSeq, fdNum, fdPtr ),
        seqSize_( seqSize )
    {
#ifdef DEBUG
        std::cout << "Creating 4-bits-per-base read buffer: "
                  << seqSize << " symbols per read" << std::endl;
#endif
    } // ~ctor
    virtual int operator[]( const int i )
    {
        char idx( seqBuf_[i >> 1] );
        if ( ( i % 2 ) != 0 ) idx >>= 4;
        idx &= ( char )0xF;
#ifdef DEBUG
        std::cout << i << " " << ( int )idx << std::endl;
#endif
        return ( int )idx;
    } // ~operator[]
    virtual void convertFromASCII( void )
    {
#ifdef DEBUG
        std::cout << seqBuf_;
#endif
        for ( int i( 0 ), j( 0 ); i < seqSize_; i++ )
        {
#ifdef DEBUG
            std::cout << seqBuf_[i];
#endif
            char code = ( char )whichPile[( int )seqBuf_[i]];
            if ( ( i % 2 ) == 0 )
                seqBuf_[j] = code;
            else
            {
                code <<= 4;
                seqBuf_[j++] |= code;
#ifdef DEBUG
                std::cout << std::endl << ( int )seqBuf_[j - 1] << std::endl;
#endif

            } // ~else
        }

    } // ~convertFromASCII


    int convertBasesToBytes( int numBases ) const
    {
        return ( ( numBases >> 1 ) + ( 1 * ( ( numBases & 0x1 ) != 0 ) ) );
    } // ~convertBasesToBytes( int numBases ) const

    const int seqSize_;
};

struct ReadBufferPrefix: public ReadBuffer4Bits
{
    ReadBufferPrefix
    ( int seqSize, int cycleNum, int fdSeq, int fdNum, int fdPtr ) :
        ReadBuffer4Bits( getPrefixBases( seqSize, cycleNum ), fdSeq, fdNum, fdPtr ),
        cycleNum_( cycleNum )
    {
#ifdef DEBUG
        std::cout << "Creating prefix read buffer for cycle " << cycleNum
                  << " of " << seqSize << std::endl;
#endif
    }

    int getPrefixBases( int seqSize, int cycleNum )
    {
        return ( seqSize - cycleNum + 2 );
    }

    virtual void sendTo( FILE *pFile )
    {
#ifdef DEBUG
        std::cout << "Outputting " << convertBasesToBytes( getPrefixBases( seqSize_, cycleNum_ + 1 ) ) << " bytes" << std::endl;
#endif

        ReadBufferBase::sendTo
        ( pFile, convertBasesToBytes( seqSize_ - 1 ) );
    } // ~sendTo


    const int cycleNum_;

};


#endif
