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

RangeStoreExternal::~RangeStoreExternal()
{
}

void RangeStoreExternal::swap( void )
{
    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "swap" << endl;
    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << fileStemIn_ << " " << fileStemOut_ << endl;
    string temp = fileStemIn_;
    fileStemIn_ = fileStemOut_;
    fileStemOut_ = temp;
    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << fileStemIn_ << " " << fileStemOut_ << endl;
} // ~swap

void RangeStoreExternal::setPortion( int pileNum, int portionNum )
{
    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "set portion " << alphabet[pileNum] << alphabet[portionNum] << endl;
    stateIn_.clear();
    //    if (stateIn_.pFile_!=NULL) fclose(stateIn_.pFile_);
    string fileName;
    getFileName( fileStemIn_, pileNum, portionNum, fileName );
    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "Made input file name " << fileName << endl;
    stateIn_.pFile_ = TemporaryRamFile::fopen( fileName.c_str(), "rb" );
    if ( stateIn_.pFile_ == NULL )
    {
        Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "Warning: no file " << fileName
                                                << " found, presuming no ranges of interest in this region"
                                                << endl;
    }
}

bool RangeStoreExternal::getRange( RangeState &stateFileIn, Range &thisRange )
{
    thisRange.clear();
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "get range: " << fileStemIn_ << endl;
    if ( !stateFileIn.good() )
        return false;
    else
    {
        stateFileIn >> thisRange;
        return stateFileIn.good();
    }
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
        stateInForComparison_[pileNum][portionNum].pFile_ = TemporaryRamFile::fopen( fileName.c_str(), "rb" );
        stateInForComparison_[pileNum][portionNum].lastProcessedPos_ = 0;

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
#ifdef PROPAGATE_PREFIX
                      << " " << r.word_
#endif
                      << " " << r.pos_ << " " << r.num_ << endl;
    }

    if ( stateOut_[pileNum][portionNum].pFile_ == NULL )
    {
        string fileName;
        getFileName( fileStemOut_, pileNum, portionNum, fileName );
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Made output file name " << fileName << endl;

        stateOut_[pileNum][portionNum].pFile_ = TemporaryRamFile::fopen( fileName.c_str(), "wb", static_cast<uint64_t>( TemporaryFilesManager::get().ramLimitMB_ * 1024 * ( 1024 / 64 ) * 0.5 ) ); // Reserves half of the available RAM for temporary files

#ifdef PROPAGATE_PREFIX
        stateOut_[pileNum][portionNum].lastProcessedPos_ = 0;
#endif
    }

    stateOut_[pileNum][portionNum] << r;
}


void RangeStoreExternal::clear( bool doDeleteFiles )
{

    string fileName;

    for ( int i( 1 ); i < alphabetSize; ++i )
    {
        for ( int j( 1 ); j < alphabetSize; ++j )
        {
            stateOut_[i][j].clear();
            stateInForComparison_[i][j].clear();

            if ( doDeleteFiles )
            {
                getFileName( fileStemIn_, i, j, fileName );
                if ( TemporaryRamFile::remove( fileName.c_str() ) != 0 )
                {
                    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Could not remove " << fileName << endl;
                }
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
