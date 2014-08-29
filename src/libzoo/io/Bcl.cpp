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

#include "Bcl.hh"

#include "libzoo/cli/Common.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <unistd.h>

#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>

#ifdef _OPENMP
#include <omp.h>
#endif //ifdef _OPENMP

//Removed due to LGPL license: #include "../redist/gzstream.hh"

using namespace std;


BclRunFolder::BclRunFolder( const string &runFolder, const string &laneFormat, const string &tileFormat )
    : runFolder_( runFolder )
    , cycleCount_( 0 )
    , laneAndTileIndex_( -1 )
{
#ifdef _OPENMP
    omp_set_nested( 1 );
#endif //ifdef _OPENMP

    generateLanesAndTilesList( laneFormat, tileFormat );
}

void BclRunFolder::generateLanesAndTilesList( const string &laneFormat, const string &tileFormat )
{
    vector<string> lanes;
    vector<string> tiles;

    // Lanes discovery
    if ( !laneFormat.empty() )
    {
        lanes.push_back( laneFormat );
    }
    else
    {
        for ( unsigned laneCount = 0; laneCount < 100; ++laneCount )
        {
            stringstream laneName;
            laneName << "L0" << ( laneCount % 10 ) << ( laneCount / 10 );
            string filename = runFolder_ + "/Data/Intensities/BaseCalls/" + laneName.str();
            if ( doesFileExist( filename ) )
            {
                clog << "Discovered lane " << filename << endl;
                lanes.push_back( laneName.str() );
            }
        }
        if ( lanes.empty() )
        {
            cerr << "Error: No lane found in " << runFolder_ << "/Data/Intensities/BaseCalls/" << endl;
        }
    }

    // Tiles discovery
    if ( !tileFormat.empty() )
    {
        tiles.push_back( tileFormat );
        lanesAndTiles_.push_back( make_pair( lanes[0], tileFormat ) );
    }
    else
    {
        for ( unsigned laneNum = 0; laneNum < lanes.size(); ++laneNum )
        {
            DIR *dp;
            dp = opendir( ( runFolder_ + "/Data/Intensities/BaseCalls/" + lanes[laneNum] + "/C1.1/" ).c_str() );

            if ( dp )
            {
                struct dirent *ep;
                while ( ( ep = readdir( dp ) ) != 0 )
                {
                    string filename( ep->d_name );
                    if ( endsWith( filename, ".stats" ) )
                    {
                        string tileName = filename.substr( 0, filename.length() - 6 );
                        clog << "Discovered tile " << tileName << endl;
                        tiles.push_back( string( ep->d_name ) );
                        lanesAndTiles_.push_back( make_pair( lanes[laneNum], tileName ) );
                    }
                }
                closedir( dp );
            }
            if ( tiles.empty() )
            {
                cerr << "Error: No tile found in " << runFolder_ << "/Data/Intensities/BaseCalls/" << endl;
            }
        }
    }

    clog << "Found " << lanesAndTiles_.size() << " tiles" << endl;
}

unsigned int BclRunFolder::getCycleCount()
{
    if ( cycleCount_ == 0 )
    {
        bclFiles_.clear();
        // Count how many cycles are available
        for ( cycleCount_ = 1; ; ++cycleCount_ )
        {
            stringstream filenameBase;
            if ( !lanesAndTiles_.empty() )
            {
                const string &lane = lanesAndTiles_[0].first;
                const string &tile = lanesAndTiles_[0].second;
                filenameBase << runFolder_ << "/Data/Intensities/BaseCalls/" << lane << "/C" << cycleCount_ << ".1/" << tile << ".bcl";
            }
            else
                filenameBase << runFolder_ << "/Data/Intensities/BaseCalls/L001/C" << cycleCount_ << ".1/s_1_1101.bcl";
            string filename1 = filenameBase.str();
            string filename2 = filenameBase.str() + ".gz";
            if ( doesFileExist( filename1 ) )
            {
            }
            else if ( doesFileExist( filename2 ) )
            {
            }
            else
            {
                --cycleCount_;
                break;
            }
        }
        clog << "BclRunFolder: Found " << cycleCount_ << " cycles" << endl;;
    }
    return cycleCount_;
}

bool BclRunFolder::getNextLaneAndTileNames( string &nextLane, string &nextTile )
{
    ++laneAndTileIndex_;
    if ( laneAndTileIndex_ < ( int )lanesAndTiles_.size() )
    {
        nextLane = lanesAndTiles_[laneAndTileIndex_].first;
        nextTile = lanesAndTiles_[laneAndTileIndex_].second;
        return true;
    }
    return false;
}

void BclRunFolder::initReader( const string &lane, const string &tile, const unsigned int firstCycle, const unsigned int lastCycle )
{
    lane_ = lane;
    tile_ = tile;
    firstCycle_ = firstCycle;
    lastCycle_ = lastCycle;

    // Close previous files if needed
    {
        filterFile_.close();
        for ( unsigned i = 0; i < bclFiles_.size(); ++i )
        {
            if ( bclFiles_[i] )
            {
                delete bclFiles_[i];
                bclFiles_[i] = 0;
            }
        }
    }

    // Open filter file
    {
        stringstream filenameBase;
        filenameBase << runFolder_ << "/Data/Intensities/BaseCalls/" << lane_ << "/" << tile_ << ".filter";
        string filename = filenameBase.str();
        cerr << "BclRunFolder: Opening " << filename << endl;
        filterFile_.open( filename.c_str(), ios_base::binary );

        // Check version
        unsigned int version;
        filterFile_.read( reinterpret_cast<char *>( &version ), 4 );
        filterFile_.read( reinterpret_cast<char *>( &version ), 4 );
        if ( version != 3 )
        {
            cerr << "Error: We only support filter files version 3. " << filename << " reports version " << version << endl;
            exit ( -1 );
        }
    }

    // Extract length
    filterFile_.read( reinterpret_cast<char *>( &readCount_ ), 4 );

    // Open BCL files
    bclFiles_.resize( lastCycle_ - firstCycle_ + 1 );
    #pragma omp parallel for
    for ( unsigned i = firstCycle_; i <= lastCycle_; ++i )
    {
        stringstream filenameBase;
        filenameBase << runFolder_ << "/Data/Intensities/BaseCalls/" << lane_ << "/C" << i << ".1/" << tile_ << ".bcl";

        string filename = filenameBase.str();
        istream *newFile = new ifstream( filename.c_str(), ios_base::binary );
/*
        // Trying .bcl.gz files, removed due to igzstream's LGPL license
        if ( !newFile->good() )
        {
            filename += ".gz";
            newFile = new igzstream( filename.c_str() );
        }
*/
        if ( !newFile->good() )
            #pragma omp critical (IO)
        {
            cerr << "Error opening " << filename << ". Aborting." << endl;
            exit( -1 );
        }

        bclFiles_[i - firstCycle_] = newFile;

        // Check length
        unsigned int readCountCheck;
        newFile->read( reinterpret_cast<char *>( &readCountCheck ), 4 );
        if ( readCountCheck != readCount_ )
            #pragma omp critical (IO)
        {
            cerr << "Error: " << filename << " reports " << readCountCheck << " entries instead of " << readCount_ << "." << endl;
            exit( -1 );
        }
    }

    currentReadNum_ = 0;
    lastReportedReadNum_ = 0;
    clog << "BclRunFolder: " << readCount_ << " reads detected." << endl;
}

bool BclRunFolder::getRead( vector< unsigned int > &bclValues, bool *passFilter )
// passFiltered == NULL => skip all reads where passFilter bit (from .filter file) == 0
//                      else report passFilter bit
{
    char passFilterVal;
    bool acceptableRead;
    bclValues.clear();
    do
    {
        if ( firstCycle_ == 0 )
        {
            return false;
        }

        // Read filter file
        filterFile_.get( passFilterVal );
        acceptableRead = ( passFilterVal || passFilter );

        // Read BCL files
        bool error = false;
        bclValues.resize( lastCycle_ - firstCycle_ + 1 );
        #pragma omp parallel for num_threads(6)
        for ( unsigned int i = firstCycle_; i <= lastCycle_; ++i )
        {
            unsigned char c;
            if ( bclFiles_[i - firstCycle_]->get( ( char & )c ) )
            {
                if ( acceptableRead )
                {
                    bclValues[i - firstCycle_] = c;
                }
            }
            else
            {
                #pragma omp critical (BCL_ERROR)
                error = true;
            }
        }

        if ( error )
        {
            return false;
        }
        ++currentReadNum_;
    }
    while ( !acceptableRead );

    if ( passFilter )
    {
        *passFilter = ( passFilterVal != '\0' );
    }

    return true;
}

void BclRunFolder::reportProgress( float step )
{
    float lastReport = ( float )lastReportedReadNum_ / readCount_;
    float potentialReport = ( float )currentReadNum_ / readCount_;
    if ( floorf( lastReport / step ) != floorf( potentialReport / step ) )
    {
        lastReportedReadNum_ = currentReadNum_;
        clog << "Progress: " << floorf( potentialReport * 100 ) << "% complete" << endl;
    }
}
