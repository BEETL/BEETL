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

#ifdef HAVE_BOOST

#include "FastOFStream.hh"

#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <cassert>


#define BUFSIZE (1024*1024)
//(128*1048576)

#define VERBOSE 0


FastOFStream::FastOFStream( const char *filename, const char *mode )
    : originalFilename_( filename )
    , currentlyFlushingBuffer( 0 )
    , terminate_( false )
    , flusherThread_( boost::bind( &FastOFStream::FlushThreadFunc, this ) )
    , usedAsBaseClass_( false )
    , bufferOffsetWithFile_( 0 )
{
    if ( VERBOSE ) cout << "FastOFStream: Opening " << filename << " in " << mode << " mode" << endl;
    //    assert ( mode[0] == 'w' && mode[1] == 0 );

    if ( mode[0] == 'w' ) // todo: remove read version
        fd = open( filename, O_WRONLY | O_CREAT | O_TRUNC, S_IRWXU );
    else if ( mode[0] == 'r' ) // todo: remove read version
        fd = open( filename, O_RDONLY );
    else
        assert( false && "Invalid mode to open file" );
    //    assert ( fd > 0 );
    if ( fd > 0 )
    {
        buffer = new char[BUFSIZE];
        buffer2 = new char[BUFSIZE];
        bufferPtr = buffer;
    }
    else
    {
        buffer = NULL;
        buffer2 = NULL;
    }
}

/*
void FastOFStream::operator=( const FILE* f )
{

}
*/

FastOFStream::~FastOFStream()
{
    if ( usedAsBaseClass_ )
        return;

    if ( VERBOSE ) cout << "FastOFStream destructor" << endl;

    fclose();

    /*
    cout << " press a key to continue... was file; " << originalFilename_ << endl;
    int x;
    cin >> x;
    */
}

void FastOFStream::FlushThreadFunc()
{
    while ( !terminate_ )
    {
        //            cout << "hi" << (void*)this << endl;
        {
            boost::mutex::scoped_lock lock( mutex_ );
            condition_.wait( lock );
            if ( currentlyFlushingBuffer )
            {
                //                    cout << "Flushing buffer!" << endl;
                write( fd, currentlyFlushingBuffer, BUFSIZE );

                assert ( buffer2 == 0 );
                buffer2 = currentlyFlushingBuffer;
                currentlyFlushingBuffer = 0;
            }
        }
    }
}

size_t FastOFStream::fwrite( const void *ptr0, size_t size, size_t count )
{
    const char *ptr = reinterpret_cast<const char *>( ptr0 );
    size_t n = size * count;
    assert ( n <= BUFSIZE );

    while ( n > 0 )
    {
        size_t nextN = 0;
        if ( bufferPtr + n > buffer + BUFSIZE )
        {
            size_t thisN = buffer + BUFSIZE - bufferPtr;
            nextN = n - thisN;
            n = thisN;
        }

        memcpy( bufferPtr, ptr, n );
        bufferPtr += n;
        assert ( bufferPtr <= buffer + BUFSIZE );

        // check if we reached the end of the buffer
        if ( bufferPtr == buffer + BUFSIZE )
        {
            // flush full buffer to disk

            // Wait if previous buffer is still being flushed
            while ( currentlyFlushingBuffer != 0 )
            {
                cout << "oh no i have to wait!" << endl;
                usleep( 10000 );
            }

            // swap buffers
            currentlyFlushingBuffer = buffer;
            buffer = buffer2;
            buffer2 = 0;
            bufferPtr = buffer;
            bufferOffsetWithFile_ += BUFSIZE;
            //            cout << "flushing" << (void*)this << endl;
            condition_.notify_all();
        }

        ptr += n;
        n = nextN;
    }
    return count;
}

size_t FastOFStream::ftell()
{
    size_t result = bufferOffsetWithFile_ + ( bufferPtr - buffer );
#ifdef DEBUG
    cout << "FastOFStream::ftell: bufferOffsetWithFile_=" << bufferOffsetWithFile_ << ", bufferPtr=" << ( void * )bufferPtr << ", buffer=" << ( void * )buffer << ", result=" << ( void * )result << endl;
#endif
    return result;
}

void FastOFStream::fflush()
{
    fsync( fd );
}

void FastOFStream::fclose()
{
    // Terminate flushing thread
    terminate_ = true;
    condition_.notify_all();
    if ( fd > 0 )
    {
        while ( currentlyFlushingBuffer != 0 )
        {
            if ( VERBOSE ) cout << "Waiting for flushing to complete!" << endl;
            usleep( 1000 );
        }

        // Flush the rest of the current buffer to disk
        if ( VERBOSE ) cout << "Flushing end of buffer" << endl;
        write( fd, buffer, bufferPtr - buffer );

        close( fd );
    }

    if ( VERBOSE ) cout << "Joining thread" << endl;
    if ( !flusherThread_.timed_join( boost::posix_time::milliseconds( 500 ) ) )
    {
        cerr << "Timed join expired" << endl;
        flusherThread_.interrupt();
    }

    if ( buffer )
    {
        if ( VERBOSE ) cout << "Freeing buffers" << endl;
        delete [] buffer;
        delete [] buffer2;
    }
    if ( VERBOSE ) cout << "All complete!" << endl;

    buffer = NULL;
    buffer2 = NULL;
    fd = 0;
}


#endif //ifdef HAVE_BOOST
