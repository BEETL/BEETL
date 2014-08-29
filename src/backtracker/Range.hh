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

#ifndef INCLUDED_RANGE_HH
#define INCLUDED_RANGE_HH

#include "Config.hh"
#include "Alphabet.hh"
#include "LetterCount.hh"
#include "Tools.hh"
#include "libzoo/util/TemporaryFilesManager.hh"

#include <iostream>
#include <string>


#define ENCODE_POSITIONS_AS_OFFSETS


// The following helper methods should be moved to a stream class





//#define WRITE_COMPRESSED_NUM_TEST_VERSION 1
#ifdef WRITE_COMPRESSED_NUM_TEST_VERSION

//static const LetterNumber tooBigMask( ~( ( ((LetterNumber)1) << 60 ) -1 ) );
//static const LetterNumber fitsIn4BytesMask( ~( ( ((LetterNumber)1) << 32 ) -1 ) );
//static const LetterNumber fitsIn2BytesMask( ~( ( ((LetterNumber)1) << 16 ) -1 ) );

// These are shifted by 4 bits to account for fact that num is shifted to
// accommodate extraByteCount
static const LetterNumber fitsIn4BytesMask( ~( ( ( ( LetterNumber )1 ) << 28 ) - 1 ) );
static const LetterNumber fitsIn2BytesMask( ~( ( ( ( LetterNumber )1 ) << 12 ) - 1 ) );


inline void writeCompressedNum( TemporaryFile *pFile, LetterNumber num )
{
    // Test for num being too big removed - can do this by logical AND with tooBigMask


    LetterNumber extraByteCount( sizeof( LetterNumber ) - 1 - 4 * ( ( num & fitsIn4BytesMask ) == 0 ) - 2 * ( ( num & fitsIn2BytesMask ) == 0 ) );
    // above is equivalent to the three lines below
    //    LetterNumber extraByteCount(sizeof(LetterNumber)-1); // starts at 7
    //    extraByteCount-=4*((num&fitsIn4BytesMask)==0); // becomes 3
    //    extraByteCount-=2*((num&fitsIn2BytesMask)==0); // becomes 1


    num = ( ( num << 4 ) | extraByteCount++ );

    fwrite( &num, extraByteCount, 1, pFile );

    //    if ( fwrite( &num, extraByteCount, 1, pFile ) != 1 )
    //    {
    //        cerr << "Could not write " << extraByteCount
    //             << "+1 chars to file. Aborting." << endl;
    //        exit( -1 );
    //    }
}


#else
void writeCompressedNum( TemporaryFile *pFile, LetterNumber num );
#endif

void writeBytes( TemporaryFile *pFile, bool val, int byteCount );
bool readCompressedNum( TemporaryFile *pFile, LetterNumber &num );
void readBytes( TemporaryFile *pFile, bool &val, const int byteCount );


class RangeState;

class Range
{
public:
    Range(
        const string &word,
        const LetterNumber pos,
        const LetterNumber num,
        const bool isBkptExtension = false,
        const bool hasUserData = false,
        void *userData = NULL
    ) :
        word_( word ),
        pos_( pos ),
        num_( num ),
        isBkptExtension_( isBkptExtension ),
        hasUserData_( hasUserData ),
        userData_( userData )
    {}

    Range(
        const bool hasUserData = false
    ) :
        pos_( 0 ),
        num_( 0 ),
        isBkptExtension_( false ),
        hasUserData_( hasUserData ),
        userData_( NULL )
    {}

    virtual ~Range() {}

    virtual void clear()
    {
        word_.clear();
        pos_ = 0;
        num_ = 0;
        isBkptExtension_ = 0;
        userData_ = NULL;
    }

    virtual bool writeTo( TemporaryFile *pFile, RangeState &currentState ) const;
    virtual bool readFrom( TemporaryFile *pFile, RangeState &currentState );
    virtual void prettyPrint( std::ostream &os ) const;

    string word_;
    LetterNumber pos_;
    LetterNumber num_;
    bool isBkptExtension_;
    bool  hasUserData_;
    void *userData_;
};

inline bool compareRangeByPos( const Range &r1, const Range &r2 )
{
    return r1.pos_ < r2.pos_;
}
inline bool compareRangeByPosInPair( const std::pair<Range, AlphabetSymbol> &rp1, const std::pair<Range, AlphabetSymbol> &rp2 )
{
    return rp1.first.pos_ < rp2.first.pos_;
}


//
// RangeState: this is a proxy class that sits in RangeStoreExternal
// and handles the optional storage of compressed sequence
// when PROPAGATE_SEQUENCE is set
//
struct RangeState
{
    RangeState( const bool propagateSequence );
    void clear( void );

    TemporaryFile *pFile_;
#ifdef ENCODE_POSITIONS_AS_OFFSETS
    LetterNumber lastProcessedPos_;
#endif

    RangeState &operator<<( const Range & );
    RangeState &operator>>( Range & );
    bool good();

private:
    bool propagateSequence_;
    string wordLast_;
    void addSeq( const string &seq );
    void getSeq( string &word );
};


#endif
