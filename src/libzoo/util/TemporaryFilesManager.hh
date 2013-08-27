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

#ifndef TEMPORARY_FILES_MANAGER_HH
#define TEMPORARY_FILES_MANAGER_HH

#include <cstdio>
#include <inttypes.h>
#include <memory>
#include <string>
#include <vector>

using std::string;
using std::vector;
using std::shared_ptr;

class TemporaryFilesManager
{
private:
    TemporaryFilesManager() {}
    TemporaryFilesManager( TemporaryFilesManager const & ); // no impl to avoid copies of singleton
    void operator=( TemporaryFilesManager const & ); // no impl to avoid copies of singleton

public:
    static TemporaryFilesManager &get()
    {
        static TemporaryFilesManager singleton;
        return singleton;
    }

    void setTempPath( const string &path, const bool createUniqueSubDirectory = true );
    void setRamLimit( const size_t ramLimit );
    void addFilename( const string &filename );
    void cleanupAllFiles(); // Delete all existing temporary files
    void cleanup(); // Delete all existing temporary files and temp subdirectory

    string tempPath_;
    size_t ramLimitMB_;
private:
    string tempPathParent_;
    bool tempPathWasCreated_;
    vector<string> filenames_;
};


class TemporaryFile
{
public:
    TemporaryFile( ) : f_( NULL ) {}
    TemporaryFile( const char *filename, const char *mode );
    virtual ~TemporaryFile() {}

    static TemporaryFile *fopen( const char *filename, const char *mode )
    {
        TemporaryFile *result = new TemporaryFile( filename, mode );
        if ( result->f_ )
            return result;
        else
        {
            delete result;
            return NULL;
        }
    }

    void open( const char *filename, const char *mode );

    virtual size_t read( void *ptr, size_t size, size_t nmemb )
    {
        return ::fread( ptr, size, nmemb, f_ );
    }
    virtual size_t write( const void *ptr, size_t size, size_t nmemb )
    {
        return ::fwrite( ptr, size, nmemb, f_ );
    }
    virtual size_t tell()
    {
        return ::ftell( f_ );
    }
    void flush()
    {
        if ( f_ )
            ::fflush( f_ );
    }
    virtual void close()
    {
        if ( f_ )
            ::fclose( f_ );
    }
    int fileno()
    {
        if ( f_ )
            return ::fileno( f_ );
        else
            return -1;
    }
    virtual bool eof()
    {
        if ( f_ )
            return ::feof( f_ );
        else
            return true;
    }

    friend size_t fread( void *ptr, size_t size, size_t nmemb, TemporaryFile *stream )
    {
        return stream->read( ptr, size, nmemb );
    }

    friend size_t fwrite( const void *ptr, size_t size, size_t nmemb, TemporaryFile *stream )
    {
        return stream->write( ptr, size, nmemb );
    }

    friend size_t ftell( TemporaryFile *stream )
    {
        return stream->tell();
    }
    friend void fflush( TemporaryFile *stream )
    {
        stream->flush();
    }
    friend void fclose( TemporaryFile *stream )
    {
        stream->close();
    }
    friend int fileno( TemporaryFile *stream )
    {
        return stream->fileno();
    }
    friend bool feof( TemporaryFile *stream )
    {
        return stream->eof();
    }

    static bool remove( const char *filename )
    {
        return ::remove( filename );
    }

protected:
    FILE *f_;
};


// RAM file with disk fallback
class TemporaryRamFile : public TemporaryFile
{
public:
    TemporaryRamFile( ) : TemporaryFile() {}
    TemporaryRamFile( const char *filename, const char *mode, const uint64_t maxRAM = 0 );
    virtual ~TemporaryRamFile();

    static TemporaryRamFile *fopen( const char *filename, const char *mode, const uint64_t maxRAM = 0 );
    static bool remove( const char *filename );

    virtual void close();
    virtual size_t read( void *ptr, size_t size, size_t nmemb );
    virtual size_t write( const void *ptr, size_t size, size_t nmemb );

    virtual size_t tell()
    {
        return currentPos_;
    }

    virtual bool eof();

private:
    const string filename_;
    const string mode_;
    shared_ptr< vector< char > > buf_;
    size_t currentPos_;
};


#endif // TEMPORARY_FILES_HH
