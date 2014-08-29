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

#include "Range.hh"

#include "libzoo/util/Logger.hh"

#include <cstring>

using namespace std;


bool Range::writeTo( TemporaryFile *pFile, RangeState &currentState ) const
{
    //    return true; // for debugging
#ifndef ENCODE_POSITIONS_AS_OFFSETS

    writeCompressedNum( pFile, ( pos_ << 1 ) | ( ( pos_ >> 63 ) & 1 ) ); // Move the last bit flag to bit 0 for better downstream compression
    writeCompressedNum( pFile, num_ );
    writeBytes( pFile, isBkptExtension_, 1 );

    if ( hasUserData_ )
    {
        assert( sizeof( LetterNumber ) == sizeof( void * ) );
        writeCompressedNum( pFile, reinterpret_cast<LetterNumber>( userData_ ) );
    }

    return true;

#else

    LetterNumber &lastProcessedPosOut = currentState.lastProcessedPos_;

    // Decode pos with its bit63 flag, and re-encode pos offset and flags together
    //  + move the flag bits to bit 0 for better downstream compression
    LetterNumber posWithoutBit63 = pos_ & 0x7FFFFFFFFFFFFFFFull;
    bool flag1 = pos_ >> 63;
    bool flag2 = isBkptExtension_;
    assert( posWithoutBit63 >= lastProcessedPosOut );
    LetterNumber offsetAndFlags = ( ( posWithoutBit63 - lastProcessedPosOut ) << 2 ) | ( flag1 << 1 ) | flag2;

    writeCompressedNum( pFile, offsetAndFlags );
    writeCompressedNum( pFile, num_ );

    if ( hasUserData_ )
    {
        assert( sizeof( LetterNumber ) == sizeof( void * ) );
        writeCompressedNum( pFile, reinterpret_cast<LetterNumber>( userData_ ) );
    }

    lastProcessedPosOut = posWithoutBit63 + num_;

    return true;

#endif
}

bool Range::readFrom( TemporaryFile *pFile, RangeState &currentState )
{
    if ( pFile == NULL || feof( pFile ) )
        return false;

#ifndef ENCODE_POSITIONS_AS_OFFSETS

    if ( readCompressedNum( pFile, pos_ ) )
    {
        pos_ = ( ( pos_ & 1 ) << 63 ) | ( pos_ >> 1 );
    }
    else
    {
        this->clear();
        return false;
    }

    if ( !readCompressedNum( pFile, num_ ) )
    {
        cerr << "getNum did not return true. Aborting." << endl;
        exit( -1 );
    }

    readBytes( pFile, isBkptExtension_, 1 );

    if ( hasUserData_ )
    {
        LetterNumber userData;
        readCompressedNum( pFile, userData );
        assert( sizeof( LetterNumber ) == sizeof( void * ) );
        userData_ = reinterpret_cast<void *>( userData );
    }

    return true;

#else

    LetterNumber offsetAndFlags;
    if ( readCompressedNum( pFile, offsetAndFlags ) == false )
    {
        this->clear();
        return false;
    }

    bool flag2 = offsetAndFlags & 1;
    bool flag1 = ( offsetAndFlags >> 1 ) & 1;
    LetterNumber posWithoutBit63 = ( offsetAndFlags >> 2 ) + currentState.lastProcessedPos_;

    pos_ = ( ( ( LetterNumber )flag1 ) << 63 ) | posWithoutBit63;
    isBkptExtension_ = flag2;

    if ( !readCompressedNum( pFile, num_ ) )
    {
        cerr << "getNum did not return true. Aborting." << endl;
        exit( -1 );
    }

    if ( hasUserData_ )
    {
        LetterNumber userData;
        readCompressedNum( pFile, userData );
        assert( sizeof( LetterNumber ) == sizeof( void * ) );
        userData_ = reinterpret_cast<void *>( userData );
    }

    currentState.lastProcessedPos_ = posWithoutBit63 + num_;

    /*
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "got range: " << fileStemIn_
                    << " " << thisRange.word_
                    << " " << thisRange.pos_
                    << " " << thisRange.num_ << " " << thisRange.isBkptExtension_ << endl;
    */

    return true;

#endif
}

void Range::prettyPrint( ostream &os ) const
{
    const LetterNumber posWithoutBit63 = pos_ & 0x7FFFFFFFFFFFFFFFull;
    const bool flag1 = pos_ >> 63;
    const bool flag2 = isBkptExtension_;

    os << "{pos=" << posWithoutBit63 << ", num=" << num_ << ", flag1=" << flag1 << ", flag2=" << flag2;
    if ( hasUserData_ )
        os << ", userData=" << userData_;
    os << "}";
}

RangeState::RangeState( const bool propagateSequence ) : pFile_( NULL ), propagateSequence_( propagateSequence )
{
    clear();
}

void RangeState::clear( void )
{
    if ( pFile_ != NULL ) fclose( pFile_ );
    pFile_ = NULL;

#ifdef ENCODE_POSITIONS_AS_OFFSETS
    lastProcessedPos_ = 0;
#endif

    if ( propagateSequence_ )
    {
        wordLast_ = "x";
        wordLast_[0] = notInAlphabet;
        //        wordLast_[1] = '\0';
    }
} // ~clear


void RangeState::addSeq( const string &seq )
{
#define COMPRESS_SEQ
#ifdef COMPRESS_SEQ
    if ( wordLast_[0] != notInAlphabet )
    {
        const char *pSeq( seq.c_str() );
        assert( seq.size() < 65536 );
        uint16_t pSeqLen = seq.size();
        const char *pLast( wordLast_.c_str() );
        while ( *pSeq == *pLast && pSeqLen > 0 )
        {
            ++pSeq;
            --pSeqLen;
            ++pLast;
        }
        Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "FF " << wordLast_ << " " << seq << " " << pSeq << endl;
        assert( pSeqLen < 256 );
        fwrite( &pSeqLen, sizeof( uint16_t ), 1, pFile_ );
        fwrite( pSeq, pSeqLen, 1, pFile_ );

        strcpy( ( char * )pLast, pSeq );
    }
    else
    {
#endif //ifdef COMPRESS_SEQ
        //        fprintf( pFile_, "%s\n", seq.c_str() );
        assert( seq.size() < 65536 );
        uint16_t seqLen = seq.size();
        fwrite( &seqLen, sizeof( uint16_t ), 1, pFile_ );
        fwrite( seq.c_str(), seqLen, 1, pFile_ );

#ifdef COMPRESS_SEQ
        wordLast_ = seq;
    }
#endif //ifdef COMPRESS_SEQ
}

void RangeState::getSeq( string &word )
{
#ifdef COMPRESS_SEQ
    if ( wordLast_[0] != notInAlphabet )
    {
        Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "FF " << word;
        //        word = wordLast_;
        unsigned int wordLen = wordLast_.size(); //strlen( wordLast_ ); //word.size();

        uint16_t seqLen;
        fread( &seqLen, sizeof( uint16_t ), 1, pFile_ );
        assert( seqLen < 256 );
        assert( seqLen <= wordLen );
        if ( seqLen > 0 )
        {
            if ( fread( ( char * )( wordLast_.c_str() ) + ( wordLen - seqLen ), seqLen, 1, pFile_ ) != 1 )
            {
                cerr << "Could not get valid data from file. Aborting." << endl;
                exit( -1 );
            }
        }
        //        wordLast_[ seqLen ] = 0;

        //unsigned int wordLastLen = seqLen; //strlen( wordLast_ );
        //wordLast_[ wordLastLen - 1] = '\0';
        //--wordLastLen;
        //        for ( unsigned int i( 0 ); i < wordLastLen; i++ )
        //            word[wordLen - wordLastLen + i] = wordLast_[i];
        word = wordLast_;
        Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << " -> " << word << endl;
    }
    else
#endif //ifdef COMPRESS_SEQ
    {
        uint16_t seqLen;
        fread( &seqLen, sizeof( uint16_t ), 1, pFile_ );
        assert( seqLen < 256 );
        wordLast_.resize( seqLen );
        if ( fread( ( char * )( wordLast_.c_str() ), seqLen, 1, pFile_ ) != 1 )
        {
            cerr << "Could not get valid data from file. Aborting." << endl;
            exit( -1 );
        }
        wordLast_[ seqLen ] = 0;

        word = wordLast_;
    }
}


RangeState &RangeState::operator<<( const Range &r )
{
    r.writeTo( pFile_, *this );
    if ( propagateSequence_ )
    {
        addSeq( r.word_ );
    }
    return *this;
}


RangeState &RangeState::operator>>( Range &r )
{
    if ( r.readFrom( pFile_, *this ) )
    {
        if ( propagateSequence_ )
        {
            getSeq( r.word_ );
        }
    }
    /*
        else
        {
            // Range object was cleared during readFrom
            // Error state is detectable by calling good()
        }
    */
    return *this;
}

bool RangeState::good()
{
    return ( this->pFile_ != NULL && !feof( this->pFile_ ) );
}


// Helper functions, to be moved somewhere else

#ifndef WRITE_COMPRESSED_NUM_TEST_VERSION
void writeCompressedNum( TemporaryFile *pFile, LetterNumber num )
{
    Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "AN: send " << num << endl;
    if ( ( num >> 60 ) != 0 )
    {
        cerr << "Overflow in RangeState::writeCompressedNum. Aborting." << endl;
        exit( -1 );
    }
    num <<= 4;
    LetterNumber num2 = num >> 8;
    unsigned char extraByteCount = 0;
    while ( num2 != 0 )
    {
        ++extraByteCount;
        num2 >>= 8;
    }
    if ( extraByteCount > 15 )
    {
        cerr << "Overflow(2) in RangeState::writeCompressedNum. Aborting." << endl;
        exit( -1 );
    }
    num |= extraByteCount;

    if ( fwrite( &num, 1 + extraByteCount, 1, pFile ) != 1 )
    {
        cerr << "Could not write " << extraByteCount
             << "+1 chars to file. Aborting." << endl;
        exit( -1 );
    }
}
#endif

void writeBytes( TemporaryFile *pFile, const bool val, const int byteCount )
{
    assert( byteCount == 1 && "TODO: implement more than 1 byte" );
    Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "AN: send " << val << endl;
    if ( fwrite( &val, 1, 1, pFile ) != 1 )
    {
        cerr << "Could not write 1 char to file. Aborting." << endl;
        exit( -1 );
    }
}

bool readCompressedNum( TemporaryFile *pFile, LetterNumber &num )
{
    num = 0;
    if ( fread( &num, 1, 1, pFile ) != 1 )
    {
        if ( !feof( pFile ) )
        {
            perror( "Error reading interval file: " );
        }
        return false;
    }
    Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "AN: got " << num << endl;
    int count = num & 0xF;
    if ( count > 0 && fread( ( ( char * )&num ) + 1, count, 1, pFile ) != 1 )
    {
        cerr << "Could not read " << count << " chars from file " << fileno( pFile ) << " pos " << ftell( pFile ) << " . Aborting." << endl;
        exit( -1 );
    }
    num >>= 4;

    Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "AN: decoded " << num << endl;
    return true;
}

void readBytes( TemporaryFile *pFile, bool &val, const int byteCount )
{
    assert( byteCount == 1 && "TODO: implement more than 1 byte" );
    val = 0;
    if ( fread( &val, 1, 1, pFile ) != 1 )
    {
        cerr << "Could not read 1 char to file. Aborting." << endl;
        exit( -1 );
    }
}
