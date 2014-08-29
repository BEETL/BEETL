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

#ifndef FAST_OFSTREAM_HH
#define FAST_OFSTREAM_HH

#include <iostream>
#include <stdio.h>
#include <vector>

using namespace std;

#ifdef HAVE_BOOST
# include <boost/thread.hpp>
# include <boost/bind.hpp>


class FastOFStream
{
public:
    FastOFStream( )
        : fd( 0 )
        , buffer( NULL )
        , currentlyFlushingBuffer( NULL )
        , bufferPtr( NULL )
        , buffer2( NULL )
        , terminate_( false )
        , usedAsBaseClass_( true ) // <- This is the important one
        , bufferOffsetWithFile_( 0 )
    {}
    FastOFStream( const char *filename, const char *mode );
    virtual ~FastOFStream();
    virtual size_t fwrite( const void *ptr, size_t size, size_t count );
    virtual size_t ftell();
    virtual void fflush();
    virtual void fclose();

    //    void operator=( const FILE *f );

    friend size_t fwrite( const void *ptr, size_t size, size_t count, FastOFStream *stream )
    {
        return stream->fwrite( ptr, size, count );
    }

    friend size_t ftell( FastOFStream *stream )
    {
        return stream->ftell();
    }
    friend size_t fflush( FastOFStream *stream )
    {
        stream->fflush();
    }
    friend size_t fclose( FastOFStream *stream )
    {
        stream->fclose();
    }
    friend size_t fread( void *ptr, size_t size, size_t count, FastOFStream *stream )
    {
        if ( size )
            return read( stream->fd, ptr, size * count ) / size;
        else
            return 0;
    }
    friend int fileno( FastOFStream *stream ) const
    {
        return stream->fd;
    }

protected:
    string originalFilename_;
    int fd;
    char *buffer, *currentlyFlushingBuffer, *bufferPtr, *buffer2;
    bool terminate_;
#ifdef HAVE_BOOST
    boost::mutex mutex_;
    boost::condition_variable condition_;
    boost::thread flusherThread_;
#endif //HAVE_BOOST

    void FlushThreadFunc();

private:
    bool usedAsBaseClass_;
    size_t bufferOffsetWithFile_;
};


#else //HAVE_BOOST


class FastOFStream
{
public:
    FastOFStream( ) : f_( NULL ) {}
    FastOFStream( const char *filename, const char *mode )
    {
        f_ = fopen( filename, mode );
    }
    ~FastOFStream() {}
    size_t fwrite( const void *ptr, size_t size, size_t nmemb )
    {
        return ::fwrite( ptr, size, nmemb, f_ );
    }
    size_t ftell()
    {
        return ::ftell( f_ );
    }
    void fflush()
    {
        ::fflush( f_ );
    }
    void fclose()
    {
        ::fclose( f_ );
    }

    friend size_t fwrite( const void *ptr, size_t size, size_t nmemb, FastOFStream *stream )
    {
        return stream->fwrite( ptr, size, nmemb );
    }

    friend size_t ftell( FastOFStream *stream )
    {
        return stream->ftell();
    }
    friend size_t fflush( FastOFStream *stream )
    {
        stream->fflush();
    }
    friend size_t fclose( FastOFStream *stream )
    {
        stream->fclose();
    }

private:
    FILE *f_;
};


#endif //HAVE_BOOST


#endif //ifndef FAST_OFSTREAM_HH
