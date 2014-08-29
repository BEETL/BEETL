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

#include "RangeStore.hh"

#include "libzoo/util/Logger.hh"

#include <algorithm>
#include <cstring>
#include <inttypes.h>
#include <sstream>
#include <unistd.h>

using namespace std;

//#define USE_RAM_FILES_FOR_READ_MODE
//#define USE_RAM_FILES_FOR_WRITE_MODE
//#define DONT_DELETE_PREVIOUS_CYCLE_FILES

#ifdef USE_RAM_FILES_FOR_READ_MODE
# define LOCAL_DEF__TEMPORARY_FILE__READ_MODE TemporaryRamFile
#else
# define LOCAL_DEF__TEMPORARY_FILE__READ_MODE TemporaryFile
#endif

#ifdef USE_RAM_FILES_FOR_WRITE_MODE
# define LOCAL_DEF__TEMPORARY_FILE__WRITE_MODE TemporaryRamFile
#else
# define LOCAL_DEF__TEMPORARY_FILE__WRITE_MODE TemporaryFile
#endif

bool RangeStore::isSubsetValid( const string &subset, const int cycle, const int pileNum, const int portionNum )
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
// RangeStoreExternal
//

RangeStoreExternal::RangeStoreExternal( const bool propagateSequence, const string fileStem )
    : fileStem_( fileStem )
    , stateIn_( propagateSequence )
    , stateOut_( alphabetSize, vector< RangeState >( alphabetSize, RangeState( propagateSequence ) ) )
    , stateInForComparison_( alphabetSize, vector< RangeState >( alphabetSize, RangeState( propagateSequence ) ) )
{
    setCycleNum( 0 );

    string fileName;
    for ( int i( 0 ); i < alphabetSize; ++i )
    {
        for ( int j( 0 ); j < alphabetSize; ++j )
        {
            stateOut_[i][j].clear();
            getFileName( fileStemIn_, i, j, fileName );
            if ( LOCAL_DEF__TEMPORARY_FILE__WRITE_MODE::remove( fileName.c_str() ) == 0 )
            {
                Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Removed " << fileName << endl;
            }
            getFileName( fileStemOut_, i, j, fileName );
            if ( LOCAL_DEF__TEMPORARY_FILE__WRITE_MODE::remove( fileName.c_str() ) == 0 )
            {
                Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Removed " << fileName << endl;
            }
        } // ~for j
    } // ~for i
} // ~ctor

RangeStoreExternal::~RangeStoreExternal()
{
}

/*
void RangeStoreExternal::swap( void )
{
    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "swap" << endl;
    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << fileStemIn_ << " " << fileStemOut_ << endl;
    string temp = fileStemIn_;
    fileStemIn_ = fileStemOut_;
    fileStemOut_ = temp;
    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << fileStemIn_ << " " << fileStemOut_ << endl;
} // ~swap
*/
void RangeStoreExternal::setCycleNum( const int cycleNum )
{
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "RangeStoreExternal: setting cycle num " << cycleNum << endl;
    //    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << fileStemIn_ << " " << fileStemOut_ << endl;
    //    const string fileStem_ = "compareIntervals";
    {
        ostringstream oss;
        oss << fileStem_ << "_cycle" << ( cycleNum - 1 );
        fileStemIn_ = oss.str();
    }
    {
        ostringstream oss;
        oss << fileStem_ << "_cycle" << cycleNum;
        fileStemOut_ = oss.str();
    }
    //    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << fileStemIn_ << " " << fileStemOut_ << endl;
} // ~swap

void RangeStoreExternal::setPortion( int pileNum, int portionNum )
{
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "set portion " << alphabet[pileNum] << alphabet[portionNum] << endl;
    stateIn_.clear();
    //    if (stateIn_.pFile_!=NULL) fclose(stateIn_.pFile_);
    string fileName;
    getFileName( fileStemIn_, pileNum, portionNum, fileName );
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Opening input portion file " << fileName << endl;
    stateIn_.pFile_ = LOCAL_DEF__TEMPORARY_FILE__READ_MODE::fopen( fileName.c_str(), "rb" );
    if ( stateIn_.pFile_ == NULL )
    {
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Warning: no file " << fileName
                << " found, presuming no ranges of interest in this region"
                << endl;
    }
}

void RangeStoreExternal::deleteInputPortion( int i, int j )
{
    stateInForComparison_[i][j].clear();

    string fileName;
    getFileName( fileStemIn_, i, j, fileName );
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Deleting input portion file " << fileName << endl;
#ifndef DONT_DELETE_PREVIOUS_CYCLE_FILES
    if ( LOCAL_DEF__TEMPORARY_FILE__READ_MODE::remove( fileName.c_str() ) != 0 )
    {
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Could not remove " << fileName << endl;
    }
#endif
    stateIn_.clear();
}

bool RangeStoreExternal::getRange( RangeState &stateFileIn, Range &thisRange )
{
    bool success;
    thisRange.clear();
    /* seems more optimised by removing this
        if ( !stateFileIn.good() )
        {
            success = false;
        }
        else
    */
    {
        stateFileIn >> thisRange;
        success = stateFileIn.good();
    }

    Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
    {
        Logger::out() << "get range: " << fileStemIn_ << " => ";
        if ( success )
        {
            thisRange.prettyPrint(  Logger::out() );
            Logger::out() << endl;
            if ( thisRange.pos_ == 0 ) usleep( 1000000 );
        }
        else
        {
            Logger::out() << "EOF" << endl;
        }
    }

    return success;
} // ~getRange

bool RangeStoreExternal::getRange( Range &thisRange )
{
    return getRange( stateIn_, thisRange );
}

bool RangeStoreExternal::isRangeKnown( const Range &r, const int pileNum, const int portionNum, const string &subset, const int cycle )
{
    if ( !isSubsetValid( subset, cycle, pileNum, portionNum ) )
        return true;

    // Checks whether the interval already existed at the previous cycle
    if ( stateInForComparison_[pileNum][portionNum].pFile_ == NULL )
    {
        string fileName;
        getFileName( fileStemIn_, pileNum, portionNum, fileName );
        stateInForComparison_[pileNum][portionNum].pFile_ = LOCAL_DEF__TEMPORARY_FILE__READ_MODE::fopen( fileName.c_str(), "rb" );
#ifdef ENCODE_POSITIONS_AS_OFFSETS
        stateInForComparison_[pileNum][portionNum].lastProcessedPos_ = 0;
#endif

        lastRangeReadForComparison_[pileNum][portionNum].pos_ = 0;
        lastRangeReadForComparison_[pileNum][portionNum].num_ = 0;
    }
    if ( stateInForComparison_[pileNum][portionNum].pFile_ )
    {
        while ( lastRangeReadForComparison_[pileNum][portionNum].pos_ < r.pos_ )
        {
            if ( !getRange( stateInForComparison_[pileNum][portionNum], lastRangeReadForComparison_[pileNum][portionNum] ) )
                lastRangeReadForComparison_[pileNum][portionNum].pos_ = maxLetterNumber;
        }
        if ( r.pos_ == lastRangeReadForComparison_[pileNum][portionNum].pos_ &&
             r.num_ == lastRangeReadForComparison_[pileNum][portionNum].num_ )
        {
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                Logger::out() << "Range detected as already processed: " << r.pos_ << " " << r.num_ << endl;
            }
            return true;
        }
    }
    return false;
}

void RangeStoreExternal::addRange( const Range &r, const int pileNum, const int portionNum, const string &subset, const int cycle )
{
    if ( !isSubsetValid( subset, cycle, pileNum, portionNum ) )
        return;

    Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
    {
        Logger::out() << "set range: " << fileStemOut_ << " " << alphabet[pileNum] << " " << alphabet[portionNum]
                      << " " << r.word_
                      << " ";
        r.prettyPrint( Logger::out() );
        Logger::out() << endl;
    }

    if ( stateOut_[pileNum][portionNum].pFile_ == NULL )
    {
        string fileName;
        getFileName( fileStemOut_, pileNum, portionNum, fileName );
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Made output file name " << fileName << endl;

#ifdef USE_RAM_FILES_FOR_WRITE_MODE
        stateOut_[pileNum][portionNum].pFile_ = LOCAL_DEF__TEMPORARY_FILE__WRITE_MODE::fopen( fileName.c_str(), "wb", static_cast<uint64_t>( TemporaryFilesManager::get().ramLimitMB_ * 1024 * ( 1024 / 64 ) * 0.5 ) ); // Reserves half of the available RAM for temporary files
#else
        stateOut_[pileNum][portionNum].pFile_ = LOCAL_DEF__TEMPORARY_FILE__WRITE_MODE::fopen( fileName.c_str(), "wb" );
#endif

        //#ifdef PROPAGATE_SEQUENCE
        stateOut_[pileNum][portionNum].lastProcessedPos_ = 0;
        //#endif
    }

    stateOut_[pileNum][portionNum] << r;
}


void RangeStoreExternal::addOutOfOrderRange( const Range &r, const int pileNum, const AlphabetSymbol portionNum, const string &subset, const int cycle )
{
    assert( pileNum == 0 && "only used for reordering $ pile after using end-pos file permutation" );
    assert( subset.empty() && "todo: implement with subset" );
    outOfOrderRangesForPile0_.push_back( make_pair( r, portionNum ) );
}


void RangeStoreExternal::clear( bool doDeleteFiles )
{
    // Flush out-of-order ranges if necessary
    if ( !outOfOrderRangesForPile0_.empty() )
    {
        std::sort( outOfOrderRangesForPile0_.begin(), outOfOrderRangesForPile0_.end(), compareRangeByPosInPair );
        for ( pair< Range, AlphabetSymbol > &rp : outOfOrderRangesForPile0_ )
        {
            Range &r = rp.first;
            //            AlphabetSymbol portionNum = rp.second;
            addRange( r, 0, 0/*portionNum*/, "", 0 ); //subset_, cycle_ ); // subset not implemented, cycle only used with subset
        }
        outOfOrderRangesForPile0_.clear();
    }

    // Clean up
    for ( int i( 0 ); i < alphabetSize; ++i )
    {
        for ( int j( 0 ); j < alphabetSize; ++j )
        {
            stateOut_[i][j].clear();
            stateInForComparison_[i][j].clear();

            if ( doDeleteFiles )
            {
#ifndef DONT_DELETE_PREVIOUS_CYCLE_FILES
                string fileName;
                getFileName( fileStemIn_, i, j, fileName );
                if ( LOCAL_DEF__TEMPORARY_FILE__READ_MODE::remove( fileName.c_str() ) != 0 )
                {
                    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Could not remove " << fileName << endl;
                }
#endif
            }
        }
    }
}

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
