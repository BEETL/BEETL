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

#include "RangeStore.hh"

#include "libzoo/util/Logger.hh"

#include <cstring>
#include <inttypes.h>
#include <unistd.h>

using namespace std;


bool RangeStore::isSubsetValid( const string &subset, const int cycle, const int pileNum, const int portionNum, const string &seq )
{
    switch ( subset.size() )
    {
        case 0:
            return true;

        case 1:
            if ( cycle == 1 && subset[subset.size() - cycle] != alphabet[portionNum] )
                return false;
            return true;

        default:
            if ( cycle < ( int )subset.size() && cycle >= 1 )
            {
                if ( subset[subset.size() - cycle - 1] != alphabet[pileNum] || subset[subset.size() - cycle] != alphabet[portionNum] )
                    return false;
            }
    }
    return true;
}

//
// RangeStoreRAM - hold BWT intervals in guess where
//

RangeStoreRAM::RangeStoreRAM() : pThis( &r1 ), pNext( &r2 ) {}

RangeStoreRAM::~RangeStoreRAM() {}

void RangeStoreRAM::swap( )
{
    pTemp = pNext;
    pNext = pThis;
    pThis = pTemp;
} // ~swap


void RangeStoreRAM::setPortion( int pileNum, int portionNum )
{
    i_ = ( *pThis )[pileNum][portionNum].begin();
    end_ = ( *pThis )[pileNum][portionNum].end();
} // ~setPortion

bool RangeStoreRAM::getRange( Range &thisRange )
{
    if ( i_ != end_ )
    {
        thisRange = *( i_++ );
        return true;
    }
    else return false;
}  // ~getRange


void RangeStoreRAM::addRange( const int pileNum, const int portionNum, const string &seq,
                              const LetterNumber pos, const LetterNumber num, const bool flags, const string &subset, const int cycle )
{
    assert( subset.empty() && "todo" );
    ( *pNext )[pileNum][portionNum].push_back( Range( seq, pos, num ) );
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "add range: " << alphabet[pileNum] << " " << alphabet[portionNum]
            << " " << seq << " " << pos << " " << num << endl;
} // ~addRange

// clear range store ready for next iter
void RangeStoreRAM::clear( void )
{
    for ( int i( 1 ); i < alphabetSize; ++i )
    {
        for ( int j( 1 ); j < alphabetSize; ++j )
        {
            ( *pThis )[i][j].clear();
        } // ~for j
    } // ~for i
} // ~clear

//
// RangeState: this is a proxy class that sits in RangeStoreExternal
// and handles the conversion between compressed representations on file
// and the structs the code expects
//

RangeState::RangeState() : pFile_( NULL )
{
    clear();
}

void RangeState::clear( void )
{
    //    wordLast_.clear();
    // next line makes sure string comparison fails at first char
    //  wordLast_+=notInAlphabet;
    wordLast_[0] = notInAlphabet;
    wordLast_[1] = '\0';
    //    posLast_ = 0;
    lastProcessedPos_ = 0;
    if ( pFile_ != NULL ) fclose( pFile_ );
    pFile_ = NULL;
} // ~clear

void RangeState::addSeq( const string &seq )
{
    //#define COMPRESS_SEQ
#ifdef COMPRESS_SEQ
    if ( wordLast_[0] != notInAlphabet )
    {
        const char *pSeq( seq.c_str() );
        assert( seq.size() < 65536 );
        uint16_t pSeqLen = seq.size();
        char *pLast( wordLast_ );
        while ( *pSeq == *pLast && pSeqLen > 0 )
        {
            ++pSeq;
            --pSeqLen;
            ++pLast;
        }
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "FF " << wordLast_ << " " << seq << " " << pSeq << endl;
        assert( pSeqLen < 256 );
        fwrite( &pSeqLen, sizeof( uint16_t ), 1, pFile_ );
        fwrite( pSeq, pSeqLen, 1, pFile_ );

        strcpy( pLast, pSeq );
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
        strcpy( wordLast_, seq.c_str() );
    }
#endif //ifdef COMPRESS_SEQ
}

void RangeState::getSeq( string &word )
{
#ifdef COMPRESS_SEQ
    if ( wordLast_[0] != notInAlphabet )
    {
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "FF " << word;
        word = wordLast_;
        unsigned int wordLen = strlen( wordLast_ ); //word.size();

        uint16_t seqLen;
        fread( &seqLen, sizeof( uint16_t ), 1, pFile_ );
        assert( seqLen < 256 );
        assert( seqLen <= wordLen );
        if ( seqLen == 0 || fread( wordLast_ + ( wordLen - seqLen ), seqLen, 1, pFile_ ) != 1 )
        {
            cerr << "Could not get valid data from file. Aborting." << endl;
            exit( -1 );
        }
        //        wordLast_[ seqLen ] = 0;

        unsigned int wordLastLen = seqLen; //strlen( wordLast_ );
        //wordLast_[ wordLastLen - 1] = '\0';
        //--wordLastLen;
        //        for ( unsigned int i( 0 ); i < wordLastLen; i++ )
        //            word[wordLen - wordLastLen + i] = wordLast_[i];
        word = wordLast_;
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << " -> " << word << endl;
    }
    else
    {
#endif //ifdef COMPRESS_SEQ
        uint16_t seqLen;
        fread( &seqLen, sizeof( uint16_t ), 1, pFile_ );
        assert( seqLen < 256 );
        if ( fread( wordLast_, seqLen, 1, pFile_ ) != 1 )
        {
            cerr << "Could not get valid data from file. Aborting." << endl;
            exit( -1 );
        }
        wordLast_[ seqLen ] = 0;

        word = wordLast_;
#ifdef COMPRESS_SEQ
    }
#endif //ifdef COMPRESS_SEQ
}

void RangeState::addNum( LetterNumber num )
{
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "AN: send " << num << endl;
    if ( ( num >> 60 ) != 0 )
    {
        cerr << "Overflow in RangeState::addNum. Aborting." << endl;
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
        cerr << "Overflow(2) in RangeState::addNum. Aborting." << endl;
        exit( -1 );
    }
    num |= extraByteCount;

    if ( fwrite( &num, 1 + extraByteCount, 1, pFile_ ) != 1 )
    {
        cerr << "Could not write " << ( int )extraByteCount
             << "+1 chars to file. Aborting." << endl;
        exit( -1 );
    }
}

bool RangeState::getNum( LetterNumber &num )
{
    num = 0;
    if ( fread( &num, 1, 1, pFile_ ) != 1 )
        return false;
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "AN: got " << num << endl;
    int count = num & 0xF;
    assert( count < ( int )sizeof( LetterNumber ) );
    if ( count > 0 && fread( ( ( char * )&num ) + 1, count, 1, pFile_ ) != 1 )
    {
        cerr << "Could not read " << count << " chars from file " << fileno( pFile_ ) << " pos " << ftell( pFile_ ) << " . Aborting." << endl;
        cerr << "Press any key to exit" << endl;
        for ( ;; )
        {
            sleep( 1 );
            if ( getchar() == '\n' )
                break;
        }
        exit( -1 );
    }
    num >>= 4;

    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "AN: decoded " << num << endl;
    return true;
}

void RangeState::addFlag( bool flag )
{
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "AN: send " << flag << endl;
    if ( fwrite( &flag, 1, 1, pFile_ ) != 1 )
    {
        cerr << "Could not write 1 char to file. Aborting." << endl;
        exit( -1 );
    }
}

bool RangeState::getFlag( bool &flag )
{
    flag = 0;
    return ( fread( &flag, 1, 1, pFile_ ) == 1 );
}


//
// RangeStoreExternal
//

RangeStoreExternal::RangeStoreExternal( const string fileStemIn,
                                        const string fileStemOut ) :
    fileStemIn_( fileStemIn ), fileStemOut_( fileStemOut ), stateIn_()
{
    string fileName;
    for ( int i( 1 ); i < alphabetSize; ++i )
    {
        for ( int j( 1 ); j < alphabetSize; ++j )
        {
            stateOut_[i][j].clear();
            // stateOut_[i][j].pFile_=NULL;
            // (*pThis)[i][j].clear();
            getFileName( fileStemIn_, i, j, fileName );
            if ( TemporaryRamFile::remove( fileName.c_str() ) == 0 )
            {
                Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Removed " << fileName << endl;
            }
            getFileName( fileStemOut_, i, j, fileName );
            if ( TemporaryRamFile::remove( fileName.c_str() ) == 0 )
            {
                Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Removed " << fileName << endl;
            }
        } // ~for j
    } // ~for i
} // ~ctor

RangeStoreExternal::~RangeStoreExternal() {}

void RangeStoreExternal::swap( void )
{
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "swap " << fileStemIn_ << " " << fileStemOut_ << endl;
    string temp = fileStemIn_;
    fileStemIn_ = fileStemOut_;
    fileStemOut_ = temp;
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << fileStemIn_ << " " << fileStemOut_ << endl;
} // ~swap

void RangeStoreExternal::setPortion( int pileNum, int portionNum )
{
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "set portion " << alphabet[pileNum] << alphabet[portionNum] << endl;
    stateIn_.clear();
    //    if (stateIn_.pFile_!=NULL) fclose(stateIn_.pFile_);
    string fileName;
    getFileName( fileStemIn_, pileNum, portionNum, fileName );
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Made input file name " << fileName << endl;
    stateIn_.pFile_ = TemporaryRamFile::fopen( fileName.c_str(), "rb" );
    if ( stateIn_.pFile_ == NULL )
    {
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Warning: no file " << fileName
                << " found, presuming no ranges of interest in this region"
                << endl;
    }
    stateIn_.lastProcessedPos_ = 0;
}

bool RangeStoreExternal::getRange( RangeState &stateFileIn, Range &thisRange )
{
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "get range: " << fileStemIn_ << endl;
    LetterNumber offsetAndFlags;
    if ( stateFileIn.pFile_ == NULL || feof( stateFileIn.pFile_ ) )
        return false;
    else if ( stateFileIn.getNum( offsetAndFlags ) == false )
    {
        return false;
    }
    else
    {
        bool flag2 = offsetAndFlags & 1;
        bool flag1 = ( offsetAndFlags >> 1 ) & 1;
        LetterNumber posWithoutBit63 = ( offsetAndFlags >> 2 ) + stateFileIn.lastProcessedPos_;

        thisRange.pos_ = ( ( ( LetterNumber )flag1 ) << 63 ) | posWithoutBit63;
        thisRange.isBkptExtension_ = flag2;

        if ( !stateFileIn.getNum( thisRange.num_ ) )
        {
            cerr << "getNum did not return true. Aborting." << endl;
        }

        stateFileIn.lastProcessedPos_ = posWithoutBit63 + thisRange.num_;

#ifdef PROPAGATE_PREFIX
        stateFileIn.getSeq( thisRange.word_ );
#endif
        //      assert(fgets(buf_,256,stateFileIn.pFile_)!=NULL);
        //  buf_[strlen(buf_)-1]='\0';
        //  thisRange.word_=buf_;

        //stateFileIn.getFlag( thisRange.isBkptExtension_ );


        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "got range: " << fileStemIn_
#ifdef PROPAGATE_PREFIX
                << " " << thisRange.word_
#endif
                << " " << thisRange.pos_
                << " " << thisRange.num_ << " " << thisRange.isBkptExtension_ << endl;
        return true;
    }
} // ~getRange

bool RangeStoreExternal::getRange( Range &thisRange )
{
    return getRange( stateIn_, thisRange );
}

bool RangeStoreExternal::isRangeKnown( const int pileNum, const int portionNum, const string &seq,
                                       const LetterNumber pos, const LetterNumber num, const bool flags, const string &subset, const int cycle )
{
    if ( !isSubsetValid( subset, cycle, pileNum, portionNum, seq ) )
        return false;

    // Checks whether the interval already existed at the previous cycle
    if ( stateInForComparison_[pileNum][portionNum].pFile_ == NULL )
    {
        string fileName;
        getFileName( fileStemIn_, pileNum, portionNum, fileName );
        stateInForComparison_[pileNum][portionNum].pFile_ = TemporaryRamFile::fopen( fileName.c_str(), "rb" );
        stateInForComparison_[pileNum][portionNum].lastProcessedPos_ = 0;

        lastRangeReadForComparison_[pileNum][portionNum].pos_ = 0;
        lastRangeReadForComparison_[pileNum][portionNum].num_ = 0;
    }
    if ( stateInForComparison_[pileNum][portionNum].pFile_ )
    {
        while ( lastRangeReadForComparison_[pileNum][portionNum].pos_ < pos )
        {
            if ( !getRange( stateInForComparison_[pileNum][portionNum], lastRangeReadForComparison_[pileNum][portionNum] ) )
                lastRangeReadForComparison_[pileNum][portionNum].pos_ = maxLetterNumber;
        }
        if ( pos == lastRangeReadForComparison_[pileNum][portionNum].pos_ &&
             num == lastRangeReadForComparison_[pileNum][portionNum].num_ )
        {
            Logger_if( LOG_SHOW_IF_VERBOSE )
            {
                Logger::out() << "Range detected as already processed: " << seq << " " << pos << " " << num << endl;
            }
            return false;
        }
    }
    return true;
}

void RangeStoreExternal::addRange( const int pileNum, const int portionNum, const string &seq,
                                   const LetterNumber pos, const LetterNumber num, const bool flags, const string &subset, const int cycle )
{
    if ( !isSubsetValid( subset, cycle, pileNum, portionNum, seq ) )
        return;

    Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
    {
        Logger::out() << "set range: " << fileStemOut_ << " " << alphabet[pileNum] << " " << alphabet[portionNum] << " " << seq << " " << pos << " " << num << endl;
    }

    if ( stateOut_[pileNum][portionNum].pFile_ == NULL )
    {
        string fileName;
        getFileName( fileStemOut_, pileNum, portionNum, fileName );
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Made output file name " << fileName << endl;

        stateOut_[pileNum][portionNum].pFile_ = TemporaryRamFile::fopen( fileName.c_str(), "wb", static_cast<uint64_t>( TemporaryFilesManager::get().ramLimitMB_ * 1024 * ( 1024 / 64 ) * 0.5 ) ); // Reserves half of the available RAM for temporary files
        stateOut_[pileNum][portionNum].lastProcessedPos_ = 0;
    }

    LetterNumber &lastProcessedPosOut = stateOut_[pileNum][portionNum].lastProcessedPos_;

    // Decode pos with its bit63 flag, and re-encode pos offset and flags together
    //  + move the flag bits to bit 0 for better downstream compression
    LetterNumber posWithoutBit63 = pos & 0x7FFFFFFFFFFFFFFFull;
    bool flag1 = pos >> 63;
    bool flag2 = flags;
    assert( posWithoutBit63 >= lastProcessedPosOut );
    LetterNumber offsetAndFlags = ( ( posWithoutBit63 - lastProcessedPosOut ) << 2 ) | ( flag1 << 1 ) | flag2;

    stateOut_[pileNum][portionNum].addNum( offsetAndFlags );
    stateOut_[pileNum][portionNum].addNum( num );

    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "set range: " << fileStemOut_ << " " << alphabet[pileNum] << " " << alphabet[portionNum]
            << " seq=" << seq << " pos=" << posWithoutBit63 << " length=" << num << " flag1=" << flag1 << " flag2=" << flag2 << " lastProcessedPosOut=" << lastProcessedPosOut << endl;

    lastProcessedPosOut = posWithoutBit63 + num;

#ifdef PROPAGATE_PREFIX
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "add seq " << seq.c_str() << " " << strlen( seq.c_str() ) << endl;
    //    fputs(seq.c_str(),fileOut_[pileNum][portionNum]);
    //    fprintf(stateOut_[pileNum][portionNum].pFile_,"%s\n",seq.c_str());
    stateOut_[pileNum][portionNum].addSeq( seq );
#endif

    //    stateOut_[pileNum][portionNum].addFlag( flags );
}

void RangeStoreExternal::clear( void )
{

    string fileName;

    for ( int i( 1 ); i < alphabetSize; ++i )
    {
        for ( int j( 1 ); j < alphabetSize; ++j )
        {
            // if (stateOut_[i][j].pFile_!=NULL)
            //  fclose( stateOut_[i][j].pFile_ );
            // stateOut_[i][j].pFile_=NULL;
            stateOut_[i][j].clear();
            stateInForComparison_[i][j].clear();

            getFileName( fileStemIn_, i, j, fileName );
            if ( TemporaryRamFile::remove( fileName.c_str() ) != 0 )
            {
                Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Could not remove " << fileName << endl;
            }

            // (*pThis)[i][j].clear();
        } // ~for j
    } // ~for i
    //    stateIn_.lastProcessedPos_ = 0;
} // ~clear

void RangeStoreExternal::getFileName( const string &stem, const int pile, const int portion,
                                      string &fileName )
{
    fileName = stem;
    fileName += '-';
    fileName += '0';
    fileName += ( char )( 48 + pile );
    fileName += '-';
    fileName += '0';
    fileName += ( char )( 48 + portion );
    if ( pile > 9 || portion > 9 )
    {
        cerr << "Alphabet seems to be larger than 9 chars. Aborting." << endl;
        exit( -1 );
    }

}

