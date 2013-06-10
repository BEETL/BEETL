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

#include "TemporaryFilesManager.hh"

#include "Logger.hh"
#include "Types.hh"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;


void TemporaryFilesManager::setTempPath( const string &path, const bool createUniqueSubDirectory )
{
    assert( tempPath_.empty() );

    if ( !createUniqueSubDirectory )
    {
        tempPath_ = path;
        tempPathWasCreated_ = false;
    }
    else
    {
        // Temporary change to the destination directory in order to use mkdtemp function
        char oldPath[10000];
        if ( getcwd( oldPath, 10000 ) == NULL )
        {
            cerr << "Error getting the current directory path" << endl;
            exit( -1 );
        }

        tempPathParent_ = ( path.empty() ? "" : ( path + "/" ) ) + "BEETL-Temp";
        mkdir( tempPathParent_.c_str(), 0700 ); // If this fails, it will be detected and reported below
        chdir( tempPathParent_.c_str() );

        vector<char> fullTempPath( 6, 'X' );
        fullTempPath.push_back( '\0' );
        const char *ret = mkdtemp( &fullTempPath[0] );
        if ( ret == NULL )
        {
            cerr << "Error creating temporary subdirectory in \"" << path << "\"'s subdir:" << tempPathParent_ << endl;
            exit( -1 );
        }
        tempPath_ = tempPathParent_ + "/" + string( &fullTempPath[0] );
        tempPathWasCreated_ = true;

        chdir( oldPath );
    }

    // Todo: Check that the directory exists

    Logger::out( LOG_SHOW_IF_VERBOSE ) << "Temporary files go in " << tempPath_ << endl;
}

void TemporaryFilesManager::addFilename( const string &filename )
{
    if ( std::find( filenames_.begin(), filenames_.end(), filename ) == filenames_.end() )
    {
        filenames_.push_back( filename );
    }
}

void TemporaryFilesManager::cleanupAllFiles()
{
    Logger::out( LOG_SHOW_IF_VERBOSE ) << "Removing " << filenames_.size() << " temporary files" << endl;

    for ( unsigned int i = 0; i < filenames_.size(); ++i )
    {
        const string &filename = filenames_[i];
        // cout << "Removing temporary file " << filename << endl;

        if ( remove( filename.c_str() ) != 0 )
            cerr << "TemporaryFilesManager: Error deleting file " << filename << endl;
    }
}

void TemporaryFilesManager::cleanup()
{
    cleanupAllFiles();

    if ( tempPathWasCreated_ )
    {
        if ( rmdir( tempPath_.c_str() ) != 0 )
        {
            cerr << "Warning: Could not delete temporary directory " << tempPath_ << endl;
        }
        rmdir( tempPathParent_.c_str() );
    }

    tempPath_.clear();
}


TemporaryFile::TemporaryFile( const char *filename, const char *mode )
{
    open( filename, mode );
}

void TemporaryFile::open( const char *filename, const char *mode )
{
    assert( filename[0] != '/' ); // No absolute filename allowed. I should probably test for no subdirectories as well.
    const string &tempPath = TemporaryFilesManager::get().tempPath_;
    string fullFilename;
    if ( tempPath.empty() )
        fullFilename = string( filename );
    else
        fullFilename = tempPath + "/" + string( filename );

    f_ = ::fopen( fullFilename.c_str(), mode );

    if ( mode[0] == 'w' && f_ )
        TemporaryFilesManager::get().addFilename( fullFilename );
}

static std::map<string, TemporaryRamFile * > existingRamFiles;

TemporaryRamFile::TemporaryRamFile( const char *filename, const char *mode, const uint64_t maxRAM  )
    : filename_( filename )
    , mode_( mode )
    , currentPos_( 0 )
{
    TemporaryRamFile *existingRamFile = existingRamFiles[filename_];

    switch ( mode[0] )
    {
        case 'w':
        {
            if ( existingRamFile )
            {
                delete existingRamFile;
            }

            const uint64_t reserveRAM = max<uint64_t>( maxRAM, 1 ); // minimum RAM = 1 byte
            buf_.reserve( reserveRAM );
            break;
        }

        case 'r':
        {
            if ( existingRamFile )
            {
                buf_.swap( existingRamFile->buf_ );
                delete existingRamFile;
            }
            break;
        }

        default:
            assert( false && "invalid file open mode" );
    }

    // We open the file anyway, to overwrite any existing one where needed
    open( filename_.c_str(), mode_.c_str() );

    existingRamFiles[filename_] = this;
}

TemporaryRamFile *TemporaryRamFile::fopen( const char *filename, const char *mode, const uint64_t maxRAM )
{
    TemporaryRamFile *result = new TemporaryRamFile( filename, mode, maxRAM );
    if ( result->buf_.capacity() || result->f_ )
        return result;
    else
    {
        existingRamFiles[string( filename )] = NULL;
        delete result;
        return NULL;
    }
}

size_t TemporaryRamFile::read( void *ptr, size_t size, size_t nmemb )
{
    size_t totalSize = size * nmemb;
    size_t availableInRam = ( buf_.size() > currentPos_ ) ? ( buf_.size() - currentPos_ ) : 0;
    size_t readFromRam = min( totalSize, availableInRam );
    size_t readFromDisk = totalSize - readFromRam;

    // Part to read to RAM
    if ( readFromRam )
    {
        memcpy( ptr, &buf_[currentPos_], readFromRam );
        currentPos_ += readFromRam;
    }

    // Part to read to disk
    if ( readFromDisk )
    {
        char *ptr2 = reinterpret_cast<char *>( ptr ); // Casting for pointer addition below
        if ( !f_ || fread( ptr2 + readFromRam, readFromDisk, 1, f_ ) == 0 )
            return readFromRam / size;
        currentPos_ += readFromDisk;
    }

    return size ? nmemb : 0;
}

size_t TemporaryRamFile::write( const void *ptr, size_t size, size_t nmemb )
{
    size_t totalSize = size * nmemb;
    size_t capacityLeft = buf_.capacity() - buf_.size();
    size_t writeToRam = min( totalSize, capacityLeft );
    size_t writeToDisk = totalSize - writeToRam;

    // Part to write to RAM
    if ( writeToRam )
    {
        size_t pos = buf_.size();
        buf_.resize( pos + writeToRam );
        memcpy( &buf_[pos], ptr, writeToRam );

        currentPos_ += writeToRam;
    }

    // Part to write to disk
    if ( writeToDisk )
    {
        const char *ptr2 = reinterpret_cast<const char *>( ptr ); // Casting for pointer addition below
        if ( fwrite( ptr2 + writeToRam, writeToDisk, 1, f_ ) == 0 )
            return writeToRam / size;

        currentPos_ += writeToDisk;
    }

    return size ? nmemb : 0;
}
