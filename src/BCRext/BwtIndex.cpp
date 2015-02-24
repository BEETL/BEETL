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

#include "BwtIndex.hh"

#include "BwtReader.hh"
#include "libzoo/util/Logger.hh"

#include <algorithm>
#include <unistd.h>
#include <sys/types.h>
#ifndef DONT_USE_MMAP
# include <fcntl.h>
# include <sys/mman.h>
# include <sys/stat.h>
# include <sys/types.h>
#endif

using namespace std;


template< class T >
BwtReaderIndex<T>::BwtReaderIndex( const string &filename, const string &optionalSharedMemoryPath ):
    T( filename ),
    indexFilename_( filename + ".idx" ),
    //    isNextIndex_( false ),
    pIndexFile_( NULL )
{
    //    current_.clear();
    initIndex( optionalSharedMemoryPath );
}


template< class T >
void BwtReaderIndex<T>::rewindFile( void )
{
    // rewind file and set all vars as per constructor
    //    current_.clear();
    indexNext_ = 0;
    //    initIndex();
    T::rewindFile();
} // ~rewindFile


template< class T >
LetterNumber BwtReaderIndex<T>::readAndCount( LetterCount &c, const LetterNumber numChars )
{
#ifdef DEBUG_RAC
    std::cout << "BR RLI readAndCount " << numChars << " chars " << endl;
    std::cout << "Before: " << currentPos_ << " " << ftell( T::pFile_ ) << " ";
    std::cout << c << endl;;
#endif

    LetterNumber charsLeft( numChars );
    uint32_t indexLast;

#ifdef DEBUG_RAC
    if ( indexNext_ != indexSize_ )
        assert( currentPos_ <= indexPosBwt_[indexNext_] );
#endif

    // gotcha: numChars can be set to maxLetterNumber so no expressions should
    // add to it - wraparound issues!

    // if indexLast==indexPosBwtSize we know we have gone past last index point
    // or that none are present at all
    if ( ( indexNext_ != indexSize_ )
         && ( numChars > ( indexPosBwt_[indexNext_] - T::currentPos_ ) ) )
    {
        // count interval spans at least one index point

        // how many index points does the count interval span?
        indexLast = indexNext_;
        while ( ( indexLast != indexSize_ )
                && ( numChars > ( indexPosBwt_[indexLast] - T::currentPos_ ) ) )
        {
            indexLast++;
        }
        indexLast--;

        if ( indexNext_ <= indexLast )
        {
            // more than one index point in count interval - can use index
            if ( ! ( T::currentPos_ == 0 && charsLeft >= indexPosBwt_[indexNext_] ) )
                charsLeft -= T::readAndCount( c, indexPosBwt_[indexNext_] - T::currentPos_ );
            else
            {
                charsLeft -= indexPosBwt_[0];
                c += indexCount_[0];
                if ( indexNext_ == indexLast )
                    T::seek( indexPosFile_[0], indexPosBwt_[0] );
            }
            //            assert(T::currentPos_==indexNext_);

            if ( indexNext_ != indexLast )
            {
                charsLeft -= ( indexPosBwt_[indexLast] - indexPosBwt_[indexNext_] );

                // update counts and also indexNext_
                while ( ++indexNext_ <= indexLast )
                {
                    c += indexCount_[indexNext_];
#ifdef DEBUG_RAC_VERBOSE
                    std::cout << indexNext_ << " " << indexPosBwt_[indexNext_] << " " << indexPosFile_[indexNext_] <<  " " << indexCount_[indexNext_] << endl;
#endif
                } //

                // skip to last index point and reset buffers
                T::seek( indexPosFile_[indexLast], indexPosBwt_[indexLast] );
            }
            else
            {
                assert( T::currentPos_ == indexPosBwt_[indexLast] );
                ++indexNext_;
            }
/*
            T::runLength_ = 0;
            T::pBuf_ = T::buf_ + ReadBufferSize;
            T::pBufMax_ = T::buf_ + ReadBufferSize;
*/
        } // if more than one index point
        // if we're in this clause we've gone past at least one index
        indexLast++;
        assert( indexLast <= indexSize_ );
    }
#ifdef DEBUG_RAC
    std::cout << "After (RLI) skip: " << T::currentPos_ << " " << ftell( T::pFile_ ) << " " << c << endl;
#endif
    // now read as normal until done
    charsLeft -= T::readAndCount( c, charsLeft );
    //    assert(T::currentPos_==desiredPos);
#ifdef DEBUG_RAC
    std::cout << "After (RLI) final read: " << T::currentPos_ << " " << ftell( T::pFile_ ) << " " << c << endl;
#endif
    return ( numChars - charsLeft );

}


template< class T >
void BwtReaderIndex<T>::initIndex( const string &optionalSharedMemoryPath )
{
    indexNext_ = 0;

    bool useSharedMemory = !optionalSharedMemoryPath.empty();
    string shmFilename1, shmFilename2, shmFilename3;
    if ( useSharedMemory )
    {
        string filenameWithoutSlash = T::filename_;
        std::replace( filenameWithoutSlash.begin(), filenameWithoutSlash.end(), '/', '_' );
        shmFilename1 = optionalSharedMemoryPath + "/BeetlIndexPosFile_" + filenameWithoutSlash;
        shmFilename2 = optionalSharedMemoryPath + "/BeetlIndexCount_" + filenameWithoutSlash;
        shmFilename3 = optionalSharedMemoryPath + "/BeetlIndexPosBwt_" + filenameWithoutSlash;
        if ( readWriteCheck( shmFilename1.c_str(), false, false ) )
        {
            // Load vectors from shared memory
            {
                cerr << "Info: Using mmap'ed index " << shmFilename1 << endl;
                int fd = open( shmFilename1.c_str(), O_RDONLY );
                assert( fd >= 0 );
                off_t fileSize = lseek( fd, 0, SEEK_END );
                lseek( fd, 0, SEEK_SET );

                char *mmappedFile = ( char * )mmap( NULL, fileSize, PROT_READ, MAP_SHARED /*| MAP_LOCKED | MAP_POPULATE*/, fd, 0 );
                if ( mmappedFile == ( void * ) - 1 )
                {
                    perror( "Error: Map failed" );
                    assert( false );
                }
                indexSize_ = *reinterpret_cast<uint32_t *>( mmappedFile );
                indexPosFile_ = reinterpret_cast<LetterNumber *>( mmappedFile + sizeof( indexSize_ ) );
                close( fd );
            }

            {
                int fd = open( shmFilename2.c_str(), O_RDONLY );
                assert( fd >= 0 );
                off_t fileSize = lseek( fd, 0, SEEK_END );
                lseek( fd, 0, SEEK_SET );

                char *mmappedFile = ( char * )mmap( NULL, fileSize, PROT_READ, MAP_SHARED /*| MAP_LOCKED | MAP_POPULATE*/, fd, 0 );
                if ( mmappedFile == ( void * ) - 1 )
                {
                    perror( "Error: Map failed" );
                    assert( false );
                }
                assert( indexSize_ == *reinterpret_cast<uint32_t *>( mmappedFile ) );
                indexCount_ = reinterpret_cast<LETTER_COUNT_CLASS *>( mmappedFile + sizeof( indexSize_ ) );
                close( fd );
            }
            {
                int fd = open( shmFilename3.c_str(), O_RDONLY );
                assert( fd >= 0 );
                off_t fileSize = lseek( fd, 0, SEEK_END );
                lseek( fd, 0, SEEK_SET );

                char *mmappedFile = ( char * )mmap( NULL, fileSize, PROT_READ, MAP_SHARED /*| MAP_LOCKED | MAP_POPULATE*/, fd, 0 );
                if ( mmappedFile == ( void * ) - 1 )
                {
                    perror( "Error: Map failed" );
                    assert( false );
                }
                assert( indexSize_ == *reinterpret_cast<uint32_t *>( mmappedFile ) );
                indexPosBwt_ = reinterpret_cast<LetterNumber *>( mmappedFile + sizeof( indexSize_ ) );
                close( fd );
            }
            return;
        }
    }

    LetterNumber currentPosBwt( 0 );
    uint8_t unusedAlphabetEntries( 0 );

    if ( pIndexFile_ != NULL ) fclose( pIndexFile_ );
    pIndexFile_ = fopen( indexFilename_.c_str(), "r" );
    if ( pIndexFile_ == NULL )
    {
//        Logger::error() << "Error opening index file " << indexFilename_;
//        exit( -1 );
    }
    else
    {
        // read file header
        bool isIndexV2 = false;
        uint8_t sizeOfAlphabet = 0;
        uint8_t sizeOfLetterNumber = 0;
        uint16_t sizeOfLetterCountCompact = 0;
        vector<char> buf( indexV1Header.size() );
        fread( buf.data(), indexV1Header.size(), 1, pIndexFile_ );
        if ( equal( buf.begin(), buf.end(), indexV1Header.begin() ) )
        {
            // index v1 detected
            fread( &sizeOfAlphabet, sizeof( uint8_t ), 1, pIndexFile_ );
            fread( &sizeOfLetterNumber, sizeof( uint8_t ), 1, pIndexFile_ );
            fread( &sizeOfLetterCountCompact, sizeof( uint16_t ), 1, pIndexFile_ );
        }
        else if ( equal( buf.begin(), buf.end(), indexV2Header.begin() ) )
        {
            // index v2 detected
            isIndexV2 = true;
            fread( &sizeOfAlphabet, sizeof( uint8_t ), 1, pIndexFile_ );
            fread( &sizeOfLetterNumber, sizeof( uint8_t ), 1, pIndexFile_ );
            sizeOfLetterCountCompact = sizeof( LetterCountCompact ); // unused in index v2
        }
        else
        {
            // default value from previous header-less format
            sizeOfAlphabet = 7;
            sizeOfLetterNumber = 8;
            sizeOfLetterCountCompact = 4*sizeOfAlphabet;

            rewind( pIndexFile_ );
        }
        if ( sizeOfAlphabet > alphabetSize )
        {
            Logger::error() << "WARNING: Index file " << indexFilename_ << " was built with alphabetSize == " << (int)sizeOfAlphabet << " whereas the current tools are using alphabetSize == " << alphabetSize << ".\n => You should rebuild the index files with beetl-index (or rebuild the tools using the same data widths (specified in Types.hh))." << endl;
            unusedAlphabetEntries = sizeOfAlphabet - alphabetSize;
        }
        else if ( sizeOfAlphabet < alphabetSize )
        {
            Logger::error() << "ERROR: Index file " << indexFilename_ << " was built with alphabetSize == " << (int)sizeOfAlphabet << " whereas the current tools are using alphabetSize == " << alphabetSize << ".\n => You should rebuild the index files with beetl-index (or rebuild the tools using the same data widths (specified in Types.hh))." << endl;
            exit( -1 );
        }
        if ( sizeOfLetterNumber != sizeof( LetterNumber ) )
        {
            Logger::error() << "ERROR: Index file " << indexFilename_ << " was built with sizeof(LetterNumber) == " << (int)sizeOfLetterNumber << " whereas the current tools are using sizeof(LetterNumber) == " << sizeof( LetterNumber ) << ".\n => You should rebuild the index files with beetl-index (or rebuild the tools using the same data widths (specified in Types.hh))." << endl;
            exit( -1 );
        }
        if ( sizeOfLetterCountCompact != sizeof( LetterCountCompact ) + 4 * unusedAlphabetEntries ) // allow 32 bits per unused entry to be automatically ignored
        {
            Logger::error() << "ERROR: Index file " << indexFilename_ << " was built with sizeof(LetterCountCompact) == " << sizeOfLetterCountCompact << " whereas the current tools are using sizeof(LetterCountCompact) == " << sizeof( LetterCountCompact ) << " + " << unusedAlphabetEntries << "unused alphabet entries.\n => You should rebuild the index files with beetl-index (or rebuild the tools using the same data widths (specified in Types.hh))." << endl;
            exit( -1 );
        }

        indexPosFile0_.push_back( 0 );
        while ( fread( &indexPosFile0_.back(), sizeof( LetterNumber ), 1, pIndexFile_ ) == 1 )
        {
            indexCount0_.push_back( LETTER_COUNT_CLASS() );
            if (!isIndexV2)
            {
                // In Index v1, counts were always stored using compact 32 bits values, which now need to be scaled to LETTER_COUNT_CLASS
                for (int i=0; i<alphabetSize; ++i)
                {
                    assert ( fread( &indexCount0_.back().count_[i], sizeof( uint32_t ), 1, pIndexFile_ ) == 1 );
                }
                uint32_t unusedEntry;
                for (int i=0; i<unusedAlphabetEntries; ++i)
                {
                    assert ( fread( &unusedEntry, sizeof( uint32_t ), 1, pIndexFile_ ) == 1 );
                }
            }
            else
            {
                for (int i=0; i<alphabetSize; ++i)
                {
                    int byteCount;
                    assert ( fread( &byteCount, 1, 1, pIndexFile_ ) == 1 );
                    if (byteCount)
                    {
#ifdef USE_COMPACT_STRUCTURES
                        if ( byteCount > sizeof(LetterNumberCompact) )
                        {
                            Logger::error() << "ERROR: Index file " << indexFilename_ << " contains large values. BEETL needs to be built without USE_COMPACT_STRUCTURES in BwtIndex.hh." << endl;
                            exit( -1 );
                        }
#endif
                        assert ( fread( &indexCount0_.back().count_[i], byteCount, 1, pIndexFile_ ) == 1 );
                    }
                }
            }
            for ( int i( 0 ); i < alphabetSize; i++ )
                currentPosBwt += indexCount0_.back().count_[i];
            indexPosBwt0_.push_back( currentPosBwt );
#ifdef DEBUG_RAC_VERBOSE
            cout << indexPosBwt0_.back() << " " << indexPosFile0_.back() << " " << indexCount0_.back() << endl;
#endif

            // skip unused alphabet entries, and check that they were indeed useless
            for (int i=0; i<unusedAlphabetEntries; ++i)
            {
                uint32_t unusedEntry;
                assert( fread( &unusedEntry, sizeof( uint32_t ), 1, pIndexFile_ ) == 1 );
                assert( unusedEntry == 0 && "Error: Trying to ignore an index entry, which contains a non-zero value" );
            }

            indexPosFile0_.push_back( 0 );
        } // ~while
        indexPosFile0_.pop_back();
        fclose( pIndexFile_ );
        pIndexFile_ = NULL;
    } // ~if
    indexSize_ = indexPosBwt0_.size();
    assert( indexSize_ == indexPosFile0_.size() );
    assert( indexSize_ == indexCount0_.size() );
    //  rewindFile();

    indexPosBwt_ = indexPosBwt0_.data();
    indexPosFile_ = indexPosFile0_.data();
    indexCount_ = indexCount0_.data();

    // Save vectors to shared memory
    if ( useSharedMemory && !indexPosBwt0_.empty() )
    {
        {
            ofstream os( shmFilename1 );
            if ( !os.good() )
            {
                cerr << "Error creating " << shmFilename1 << endl;
                exit( -1 );
            }
            os.write( reinterpret_cast<const char *>( &indexSize_ ), sizeof( indexSize_ ) );
            os.write( reinterpret_cast<const char *>( indexPosFile0_.data() ), indexSize_ * sizeof( indexPosFile0_[0] ) );
        }
        {
            ofstream os( shmFilename2 );
            os.write( reinterpret_cast<const char *>( &indexSize_ ), sizeof( indexSize_ ) );
            os.write( reinterpret_cast<const char *>( indexCount0_.data() ), indexSize_ * sizeof( indexCount0_[0] ) );
        }
        {
            ofstream os( shmFilename3 );
            os.write( reinterpret_cast<const char *>( &indexSize_ ), sizeof( indexSize_ ) );
            os.write( reinterpret_cast<const char *>( indexPosBwt0_.data() ), indexSize_ * sizeof( indexPosBwt0_[0] ) );
        }
    }
} // ~initIndex




// Index creation

void buildIndex( BwtReaderBase *reader0, FILE *pIndexFile, const int indexBinSize )
{
    BwtReaderRunLengthBase *reader = dynamic_cast< BwtReaderRunLengthBase* >( reader0 );
    const int runsPerChunk( indexBinSize );
    int runsThisChunk( 0 );
    LetterCount countsThisChunk;
    LetterNumber runsSoFar( 0 ), chunksSoFar( 0 );
    bool lastRun = false;

    if (reader == NULL)
    {
        Logger::out() << "Warning: cannot index file " << reader0->filename_ << endl;
        return;
    }
    reader->currentPos_ = 0;

    // Write file header
    assert( fwrite( indexV2Header.data(), indexV2Header.size(), 1, pIndexFile ) == 1 );
    uint8_t sizeOfAlphabet = alphabetSize;
    uint8_t sizeOfLetterNumber = sizeof( LetterNumber );
    fwrite( &sizeOfAlphabet, sizeof( uint8_t ), 1, pIndexFile );
    fwrite( &sizeOfLetterNumber, sizeof( uint8_t ), 1, pIndexFile );


    while ( !lastRun )
    {
        lastRun = !reader->getRun();
        if (!lastRun)
        {
            runsSoFar++;
            runsThisChunk++;

            countsThisChunk.count_[whichPile[reader->lastChar_]] += reader->runLength_;
            assert( countsThisChunk.count_[whichPile[reader->lastChar_]] >= reader->runLength_ && "Error: Overflow in buildIndex" );
            reader->currentPos_ += reader->runLength_;
        }
        if ( runsThisChunk == runsPerChunk || lastRun )
        {
#ifdef DEBUG_RAC
            cout << reader->currentPos_ << " " << runsSoFar << " " << countsThisChunk << endl;
#endif

            // don't bother writing this as can deduce by summing countsThisChunk
            //            assert
            //            ( fwrite( &reader->currentPos_, sizeof( LetterNumber ), 1, pIndexFile ) == 1 );

            LetterNumber posInFile = reader->tellg();

            assert
            ( fwrite( &posInFile, sizeof( LetterNumber ), 1, pIndexFile ) == 1 );

            // In index format v2, we write each LetterCount independently, encoding the number of bytes as first byte
            for (int i=0; i<alphabetSize; ++i)
            {
                LetterNumber val = countsThisChunk.count_[i];
                int bytesNeeded = 0;
                while (val >> (8*bytesNeeded))
                    ++bytesNeeded;
                assert( fwrite( &bytesNeeded, 1, 1, pIndexFile ) == 1 );
                if (bytesNeeded)
                    assert( fwrite( &val, bytesNeeded, 1, pIndexFile ) == 1 );
            }

            chunksSoFar++;
            runsThisChunk = 0;
            countsThisChunk.clear();
        }
    }
    cout << "buildIndex: read " << reader->currentPos_ << " bases compressed into " << runsSoFar << " runs" << " over " << reader->tellg() << " bytes." << endl;
    cout << "buildIndex: generated " << chunksSoFar << " index points." << endl;
} // ~buildIndex





// Explicit template instantiations
template class BwtReaderIndex<BwtReaderRunLength>;
template class BwtReaderIndex<BwtReaderRunLengthV3>;
