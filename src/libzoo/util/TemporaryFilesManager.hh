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

#ifndef TEMPORARY_FILES_MANAGER_HH
#define TEMPORARY_FILES_MANAGER_HH

#include <cstdio>
#include <inttypes.h>
#include <memory>
#include <stdlib.h>
#include <string>
#include <vector>

using std::string;
using std::vector;
using std::shared_ptr;


//#define DISABLE_WRITES_AND_REMOVES


class TemporaryFilesManager
{
private:
    TemporaryFilesManager() : ramLimitMB_( 0 ), tempPathWasCreated_( false ) {}
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

static const int TempFileBufSize( 32768 );

class TemporaryFile
{
public:
    TemporaryFile( )
        : f_( NULL )
#ifdef BUFFERED_WRITE_TEST_VERSION
        , buf_( TempFileBufSize )
        , p_( &buf_[0] )
        , pBufMax_( p_ + TempFileBufSize )
#endif
    {}
    TemporaryFile( const char *filename, const char *mode );
    virtual ~TemporaryFile() {}

    static TemporaryFile *fopen( const char *filename, const char *mode )
    {
#ifdef DISABLE_WRITES_AND_REMOVES
        if ( mode[0] == 'w' )
            return new TemporaryFile( "tmp", mode );
#endif
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

#ifdef BUFFERED_WRITE_TEST_VERSION
    void flushBuffer( void )
    {
        ::fwrite( &buf_[0], 1, p_ - &buf_[0], f_ );
        p_ = &buf_[0];
    }

    virtual size_t write( const void *ptr, size_t size, size_t nmemb )
    {
#ifdef DISABLE_WRITES_AND_REMOVES
        return size ? nmemb : 0;
#endif
        char *pIn( ( char * )ptr );
        for ( size_t i( 0 ); i < size * nmemb; i++ )
        {
            if ( p_ == pBufMax_ ) flushBuffer();
            *p_++ = *pIn++;
        }


        return nmemb;
    }
#else
    virtual size_t write( const void *ptr, size_t size, size_t nmemb )
    {
#ifdef DISABLE_WRITES_AND_REMOVES
        return size ? nmemb : 0;
#endif
        return ::fwrite( ptr, size, nmemb, f_ );
    }
#endif


    virtual size_t tell()
    {
        // TBD - gives wrong answers for buffered write
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
        {
#ifdef BUFFERED_WRITE_TEST_VERSION
            flushBuffer();
#endif
            ::fclose( f_ );
        }
        delete this;
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

    static bool remove( const char *filename );

protected:
    static string getFullFilename( const string &filename );

    FILE *f_;
#ifdef BUFFERED_WRITE_TEST_VERSION
    vector<char> buf_;
    char *p_;
    char *pBufMax_;
#endif //ifdef BUFFERED_WRITE_TEST_VERSION
};


struct NoInitChar
{
    char value;
    NoInitChar() {} // do nothing, especially not an initialisation
};

class NoInitCharVector
{
public:
    NoInitCharVector()
        : data( NULL )
        , capacity_( 0 )
        , size_( 0 )
    {}

    ~NoInitCharVector()
    {
        free( data );
    }

    void reserve( size_t n )
    {
        data = ( char * )malloc( n );
        capacity_ = n;
    }
    void resize( size_t n )
    {
        size_ = n;
    }
    size_t capacity() const
    {
        return capacity_;
    }
    size_t size() const
    {
        return size_;
    }

    char *data;

private:
    size_t capacity_;
    size_t size_;
};

// RAM file with disk fallback
class TemporaryRamFile : public TemporaryFile
{
public:
    TemporaryRamFile( ) : TemporaryFile(), currentPos_( 0 ) {}
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
    shared_ptr< NoInitCharVector > buf_;
    size_t currentPos_;
    //    char localBuf_[1024*1024];
};


#endif // TEMPORARY_FILES_HH
