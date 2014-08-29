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

#include "TemporaryFilesManager.hh"

#include "libzoo/util/Logger.hh"

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


void TemporaryFilesManager_cleanupAtExit()
{
    TemporaryFilesManager::get().cleanup();
}

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
        mode_t directoryCreationMode = 0777;

        // When available, add username to temporary path name to make it user-independent and to rely a bit less on the extended permissions
        char *username = getlogin();
        if (username) {
            tempPathParent_ += string("-") + string(username);
            directoryCreationMode = 0700;
        }

        mkdir( tempPathParent_.c_str(), directoryCreationMode ); // If this fails, it will be detected and reported below
        chmod( tempPathParent_.c_str(), directoryCreationMode ); // Useful in case umask affected mkdir's permissions
        if ( chdir( tempPathParent_.c_str() ) != 0 )
        {
            cerr << "Warning: Cannot enter directory " << tempPathParent_ << endl;
        }

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

        if ( chdir( oldPath ) != 0 )
        {
            cerr << "Warning: Cannot restore path to " << oldPath << endl;
        }
        atexit( TemporaryFilesManager_cleanupAtExit );
    }

    // Todo: Check that the directory exists

    Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Temporary files go in " << tempPath_ << endl;
}

void TemporaryFilesManager::setRamLimit( const size_t ramLimitMB )
{
    ramLimitMB_ = ramLimitMB;
    Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Temporary RAM files limit set to " << ramLimitMB_ << " MB" << endl;
}

void TemporaryFilesManager::addFilename( const string &filename )
{
    #pragma omp critical (ACCESS_FILE_MANAGER_FILENAMES)
    if ( std::find( filenames_.begin(), filenames_.end(), filename ) == filenames_.end() )
    {
        filenames_.push_back( filename );
    }
}

void TemporaryFilesManager::cleanupAllFiles()
{
    #pragma omp critical (ACCESS_FILE_MANAGER_FILENAMES)
    {
        Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Removing " << filenames_.size() << " temporary files" << endl;

        for ( unsigned int i = 0; i < filenames_.size(); ++i )
        {
            const string &filename = filenames_[i];
            // cout << "Removing temporary file " << filename << endl;

            if ( remove( filename.c_str() ) != 0 )
            {
                Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "TemporaryFilesManager: Could not delete file " << filename << endl;
            }
        }

        filenames_.clear();
    }
}

void TemporaryFilesManager::cleanup()
{
    cleanupAllFiles();

    if ( tempPathWasCreated_ )
    {
        int attemptRemaining = 3;
        while ( attemptRemaining-- > 0 )
        {
            if ( rmdir( tempPath_.c_str() ) != 0 )
            {
                if ( attemptRemaining == 0 )
                {
                    cerr << "Warning: Could not delete temporary directory " << tempPath_ << endl;
                    perror( "Reason: " );
                }
                else
                {
                    Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Waiting for files to synchronise before deleting temp directory. " << attemptRemaining << " attempts remaining..." << endl;
                    ( void ) system( "sync" );
                    sleep( 1 );
                }
            }
            else
                break;
        }
        rmdir( tempPathParent_.c_str() );
    }

    tempPath_.clear();
    tempPathWasCreated_ = false;
}


TemporaryFile::TemporaryFile( const char *filename, const char *mode )
#ifdef BUFFERED_WRITE_TEST_VERSION
, p_( &buf_[0] )
, pBufMax_( p_ + TempFileBufSize )
#endif
{
    open( filename, mode );
}

void TemporaryFile::open( const char *filename, const char *mode )
{
    string fullFilename = getFullFilename( filename );
    f_ = ::fopen( fullFilename.c_str(), mode );

    if ( mode[0] == 'w' && f_ )
        TemporaryFilesManager::get().addFilename( fullFilename );
}

bool TemporaryFile::remove( const char *filename )
{
#ifdef DISABLE_WRITES_AND_REMOVES
    return true;
#endif
    string fullFilename = getFullFilename( filename );
    return ::remove( fullFilename.c_str() );
}

string TemporaryFile::getFullFilename( const string &filename )
{
    assert( filename[0] != '/' ); // No absolute filename allowed. I should probably test for no subdirectories as well.
    const string &tempPath = TemporaryFilesManager::get().tempPath_;
    string fullFilename;
    if ( tempPath.empty() )
        fullFilename = string( filename );
    else
        fullFilename = tempPath + "/" + string( filename );

    return fullFilename;
}


static std::map<string, TemporaryRamFile * > existingRamFiles;

TemporaryRamFile::~TemporaryRamFile()
{
    Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "TemporaryRamFile::~TemporaryRamFile " << filename_ << " , mode=" << mode_ << endl;
}

TemporaryRamFile::TemporaryRamFile( const char *filename, const char *mode, const uint64_t maxRAM  )
    : filename_( filename )
    , mode_( mode )
    , currentPos_( 0 )
{
    Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "TemporaryRamFile::TemporaryRamFile " << filename_ << " , mode=" << mode_ << endl;

    TemporaryRamFile *existingRamFile;
    #pragma omp critical (ACCESS_EXISTING_RAM_FILES)
    existingRamFile = existingRamFiles[filename_];

    switch ( mode[0] )
    {
        case 'w':
        {
            if ( existingRamFile )
            {
                delete existingRamFile;
            }

            const uint64_t reserveRAM = max<uint64_t>( maxRAM, 1 ); // minimum RAM = 1 byte
            buf_.reset( new NoInitCharVector );
            buf_->reserve( reserveRAM );
            #pragma omp critical (ACCESS_EXISTING_RAM_FILES)
            existingRamFiles[filename_] = this;
        }
        break;

        case 'r':
        {
            if ( existingRamFile )
            {
                buf_ = existingRamFile->buf_;
                //                delete existingRamFile;
            }
        }
        break;

        default:
            assert( false && "invalid file open mode" );
    }

    // We open the file anyway, to overwrite any existing one where needed
    open( filename_.c_str(), mode_.c_str() );
}

TemporaryRamFile *TemporaryRamFile::fopen( const char *filename, const char *mode, const uint64_t maxRAM )
{
    Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "TemporaryRamFile::fopen " << filename << " , mode=" << mode << ", maxRAM=" << maxRAM << endl;

    TemporaryRamFile *result = new TemporaryRamFile( filename, mode, maxRAM );
    if ( ( result->buf_ && result->buf_->capacity() ) || result->f_ )
        return result;
    else
    {
        #pragma omp critical (ACCESS_EXISTING_RAM_FILES)
        existingRamFiles[string( filename )] = NULL;
        delete result;
        return NULL;
    }
}

bool TemporaryRamFile::remove( const char *filename )
{
    Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "TemporaryRamFile::remove " << filename << endl;

    bool ret = false;
    #pragma omp critical (ACCESS_EXISTING_RAM_FILES)
    {
        std::map<string, TemporaryRamFile * >::iterator it = existingRamFiles.find( filename );
        if ( it != existingRamFiles.end() )
        {
            TemporaryRamFile *existingRamFile = it->second;
            delete existingRamFile;
            existingRamFiles.erase( it );

            // delete real file
            const string &tempPath = TemporaryFilesManager::get().tempPath_;
            string fullFilename;
            if ( tempPath.empty() )
                fullFilename = string( filename );
            else
                fullFilename = tempPath + "/" + string( filename );

            ::remove( fullFilename.c_str() );

            ret = true;
        }
    }
    return ret;
}

void TemporaryRamFile::close()
{
    TemporaryFile::close();
    if ( mode_[0] == 'r' )
    {
        buf_.reset();
        delete this;
    }
}

size_t TemporaryRamFile::read( void *ptr, size_t size, size_t nmemb )
{
    size_t totalSize = size * nmemb;
    size_t availableInRam = ( buf_->size() > currentPos_ ) ? ( buf_->size() - currentPos_ ) : 0;
    size_t readFromRam = min( totalSize, availableInRam );
    size_t readFromDisk = totalSize - readFromRam;

    // Part to read to RAM
    if ( readFromRam )
    {
        memcpy( ptr, &( buf_->data[currentPos_] ), readFromRam );
        currentPos_ += readFromRam;
        /*
                if (currentPos_ == 0)
                {
                    memcpy( localBuf_, buf_->data, 1024*1024 );
                }

                size_t startPos = currentPos_;
                size_t endPos = currentPos_ + readFromRam;
                size_t startPosInLocalBuf = startPos % (1024*1024);
                size_t endPosInLocalBuf = startPosInLocalBuf + readFromRam;

                if (endPosInLocalBuf < 1024*1024)
                {
                    memcpy( ptr, &( localBuf_[startPosInLocalBuf] ), readFromRam );
                }
                else
                {
                    memcpy( ptr, &( localBuf_[startPosInLocalBuf] ), 1024*1024-startPosInLocalBuf );
                    startPos += 1024*1024-startPosInLocalBuf;
                    ptr = ((char*)ptr)+1024*1024;
                    endPosInLocalBuf -= 1024*1024;

                    while (endPosInLocalBuf >= 1024*1024)
                    {
                        memcpy( ptr, buf_->data + startPos , 1024*1024 );
                        startPos += 1024*1024;
                        ptr = ((char*)ptr)+1024*1024;
                        endPosInLocalBuf -= 1024*1024;
                    }
                    memcpy( localBuf_, buf_->data + startPos , 1024*1024 );
                    memcpy( ptr, localBuf_, endPosInLocalBuf );
                }
                currentPos_ = endPos;
        */
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
    size_t capacityLeft = buf_->capacity() - buf_->size();
    size_t writeToRam = min( totalSize, capacityLeft );
    size_t writeToDisk = totalSize - writeToRam;

    // Part to write to RAM
    if ( writeToRam )
    {
        size_t pos = buf_->size();
        buf_->resize( pos + writeToRam );
        memcpy( ( char * ) & ( buf_->data[pos] ), ptr, writeToRam );

        //        const char *ptr2 = static_cast<const char *>( ptr );
        //        buf_->insert( buf_->end(), ptr2, ptr2 + writeToRam );

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

bool TemporaryRamFile::eof()
{
    size_t availableInRam = ( buf_->size() > currentPos_ ) ? ( buf_->size() - currentPos_ ) : 0;

    if ( availableInRam )
        return false;
    else
        return TemporaryFile::eof();
}
