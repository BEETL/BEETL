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

#include "BwtReader.hh"

#include "BwtIndex.hh"
#include "BwtWriter.hh"
#include "CountWords.hh"
#include "LetterCount.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#ifndef DONT_USE_MMAP
# include <fcntl.h>
# include <sys/mman.h>
# include <sys/stat.h>
# include <sys/types.h>
#endif

using namespace std;


//#define DEBUG_RAC 1
//#define DEBUG_RAC_VERBOSE 1

//
// BwtReaderBase member function definitions
//

BwtReaderBase::BwtReaderBase( const string &filename ) :
    filename_( filename )
    , pFile_( fopen( filename.c_str(), "r" ) )
    , buf_( ReadBufferSize )
{
    if ( pFile_ == NULL )
    {
        cerr << "!! BwtReaderBase: failed to open file " << filename << endl;
        exit( EXIT_FAILURE );
    }
#ifdef DEBUG
    cout << "BwtReaderBase: opened file " << filename << " " << pFile_ << endl;
#endif
}

BwtReaderBase::BwtReaderBase( const BwtReaderBase &obj ):
    filename_( obj.filename_ )
    , pFile_( fopen( obj.filename_.c_str(), "r" ) )
    , buf_( obj.buf_ )
{
    if ( pFile_ == NULL )
    {
        cerr << "!! BwtReaderBase: failed to re-open file " << filename_ << endl;
        exit( EXIT_FAILURE );
    }
#ifdef DEBUG
    cout << "BwtReaderBase: re-opened file " << filename_ << " " << pFile_ << endl;
#endif

    // Restore file to same position as before
    fseek( pFile_, ftell( obj.pFile_ ), SEEK_SET );
} // ~ctor

BwtReaderBase::~BwtReaderBase()
{
    if ( pFile_ )
        fclose( pFile_ );
#ifdef DEBUG
    cout << "BwtReaderBase: closed file " << pFile_ << endl;
#endif
}



LetterNumber BwtReaderBase::readAndCount( LetterCount &c )
{
    // will call the relevant virtual function
    return readAndCount( c, maxLetterNumber );
}

LetterNumber BwtReaderBase::readAndSend( BwtWriterBase &writer )
{
    // will call the relevant virtual function
    return readAndSend( writer, maxLetterNumber );
}


//
// BwtReaderASCII member function definitions
//

void BwtReaderASCII::rewindFile( void )
{
    rewind( pFile_ );
    currentPos_ = 0;
} // ~rewindFile

LetterNumber BwtReaderASCII::tellg( void ) const
{
    return currentPos_;
} // ~tellg

int BwtReaderASCII::seek( const LetterNumber posInFile, const LetterNumber baseNumber )
{
    assert( posInFile == baseNumber );
    currentPos_ = baseNumber;
    return fseek( pFile_, posInFile, SEEK_SET );
}


LetterNumber BwtReaderASCII::readAndCount( LetterCount &c, const LetterNumber numChars )
{
#ifdef DEBUG
    std::cout << "BR ASCII readAndCount " << numChars << " chars " << endl;
#endif
    LetterNumber charsLeft( numChars ), charsToRead, charsRead;
    while ( charsLeft > 0 )
    {
        charsToRead = ( ( charsLeft > ReadBufferSize ) ? ReadBufferSize : charsLeft );
        charsRead = fread( buf_.data(), sizeof( char ), charsToRead, pFile_ );
#ifdef DEBUG
        std::cout << "Reading " << charsRead << " chars ";
#endif
        for ( LetterNumber i( 0 ); i < charsRead; i++ )
        {
#ifdef DEBUG
            std::cout << buf_[i];
#endif
            assert( whichPile[( int )buf_[i]] != nv && "Letter not in alphabet" );
            c.count_[whichPile[( int )buf_[i]]]++;
        }
#ifdef DEBUG
        std::cout << std::endl;
#endif
        charsLeft -= charsRead;
        if ( charsRead < charsToRead )
        {
            // did not get everything asked for! return num of chars actually found
            currentPos_ += ( numChars - charsLeft );
            return ( numChars - charsLeft );
        } // ~if
    } // ~while
    currentPos_ += numChars;
    return numChars;
} // ~int BwtReaderASCII::readAndCount( LetterCount& c, const LetterNumber numChars )

LetterNumber BwtReaderASCII::readAndSend( BwtWriterBase &writer, const LetterNumber numChars )
{
#ifdef DEBUG
    std::cout << "BR ASCII readAndSend " << numChars << " chars " << endl;
#endif

    LetterNumber totalRead = 0;
    // read readbufferzise bytes, if there are less bytes
    // ordered only fetch the last missing bytes
    LetterNumber readNextPass =
        ( ( numChars - totalRead ) < ReadBufferSize ) ? ( numChars - totalRead ) : ReadBufferSize;

    // std::cout << "Reading " << numChars << " chars " << endl;


    while ( totalRead < ( LetterNumber )numChars )
    {
        LetterNumber numRead = 0;

        // try to read buffersize byte from file
        numRead = fread( buf_.data(), sizeof ( uchar ), readNextPass, pFile_ );
        totalRead += numRead;
        if ( numRead == 0 ) break;

        readNextPass =
            ( ( numChars - totalRead ) < ReadBufferSize ) ? ( numChars - totalRead ) : ReadBufferSize;
        //std::cout << "next pass " << numRead << " chars " << endl;

        //    writer( buf_, numRead );
#define XXX 1
#ifdef XXX
        LetterNumber charsLeft = numRead;

        for ( LetterNumber counter = 0; counter < charsLeft; counter++ )
        {
            if ( buf_[counter] == lastChar_ )
            {
                runLength_++; // same char, increase runlength counter
            }
            else
            {
                // new char, print previous run
                writer.sendRun( lastChar_, runLength_ );
                // reset runlength to new char
                lastChar_ = buf_[counter];
                runLength_ = 1;
            }
        } // ~for
#endif

        currentPos_ += numRead;
    }
#ifdef XXX
    writer.sendRun( lastChar_, runLength_ ); // send out last run
    runLength_ = 0; // next call to his function will then again start a new run
#endif
    return totalRead;

} // ~LetterNumber BwtReaderASCII::readAndSend( BwtWriterBase& writer, const LetterNumber numChars )

LetterNumber BwtReaderASCII::operator()( char *p, LetterNumber numChars )
{
#ifdef DEBUG
    std::cout << "BR ASCII () " << numChars << " chars " << endl;
#endif
    //    std::cout << "want " << numChars << " chars" << std::endl;
    return fread( p, sizeof( char ), numChars, pFile_ );
} // ~operator()



//
// BwtReaderRunLengthBase member function definitions
//

BwtReaderRunLengthBase::BwtReaderRunLengthBase( const string &filename ):
    BwtReaderBase( filename ),
    lengths_( 256 ),
    codes_( 256 ),
    pBuf_( buf_.data() + ReadBufferSize ),
    pBufMax_( buf_.data() + ReadBufferSize ),
    finished_( false ),
    lastChar_( notInAlphabet ),
    runLength_( 0 ),
    currentPos_( 0 ),
    currentPosInFile_( 0 )
{
    for ( unsigned int i( 0 ); i < 256; i++ )
    {
        lengths_[i] = 1 + ( i >> 4 );
        uint j = ( i & 0xF );
        codes_[i] = ( j < alphabetSize ) ? alphabet[j] : notInAlphabet;
    } // ~for i
} // ~ctor

BwtReaderRunLengthBase::BwtReaderRunLengthBase( const BwtReaderRunLengthBase &obj ):
    BwtReaderBase( obj ),
    lengths_( obj.lengths_ ),
    codes_( obj.codes_ ),
    pBuf_( obj.pBuf_ ),
    pBufMax_( obj.pBufMax_ ),
    finished_( obj.finished_ ),
    lastChar_( obj.lastChar_ ),
    runLength_( obj.runLength_ ),
    currentPos_( obj.currentPos_ ),
    currentPosInFile_( obj.currentPosInFile_ )
{
} // ~ctor


void BwtReaderRunLengthBase::rewindFile( void )
{
    // rewind file and set all vars as per constructor
    rewind( pFile_ );
    runLength_ = 0;
    pBuf_ = buf_.data() + ReadBufferSize;
    pBufMax_ = buf_.data() + ReadBufferSize;
    lastChar_ = notInAlphabet;
    currentPos_ = 0;
    currentPosInFile_ = 0;
    finished_ = false;
} // ~rewindFile

LetterNumber BwtReaderRunLengthBase::tellg( void ) const
{
    return currentPosInFile_;
} // ~tellg

int BwtReaderRunLengthBase::seek( const LetterNumber posInFile, const LetterNumber baseNumber )
{
    currentPos_ = baseNumber;
    currentPosInFile_ = posInFile;

    runLength_ = 0;
    pBuf_ = buf_.data() + ReadBufferSize;
    pBufMax_ = buf_.data() + ReadBufferSize;
    lastChar_ = notInAlphabet;
    finished_ = false;

    return fseek( pFile_, posInFile, SEEK_SET );
}


LetterNumber BwtReaderRunLengthBase::readAndCount( LetterCount &c, const LetterNumber numChars )
{
#ifdef DEBUG_RAC
    std::cout << "BR RL readAndCount " << numChars << " chars " << endl;
    std::cout << "Before: " << currentPos_ << " " << ftell( pFile_ ) << " ";
    std::cout << c << endl;
#endif

    LetterNumber charsLeft( numChars );
    while ( charsLeft > runLength_ )
    {
        if (runLength_)
        {
            c.count_[whichPile[lastChar_]] += runLength_;
            charsLeft -= runLength_;
        }
#ifdef DEBUG_RAC
        std::cout << "R&C: " << currentPos_ << " " << ftell( pFile_ ) << " " << charsLeft << " " << runLength_ << " " << lastChar_ << " " << c << endl;
#endif
        if ( getRun() == false )
        {
            currentPos_ += ( numChars - charsLeft );
#ifdef DEBUG_RAC
            std::cout << "After (end): " << currentPos_ << " " << ftell( pFile_ ) << " " << c << endl;
#endif
            return ( numChars - charsLeft );
            // assert(1==0);


        } // ~if
    } // ~while

    c.count_[whichPile[lastChar_]] += charsLeft;
    runLength_ -= charsLeft;
    currentPos_ += numChars;
#ifdef DEBUG_RAC
    std::cout << "After (not at end): " << currentPos_ << " " << ftell( pFile_ ) << " " << c << endl;
#endif

    return numChars;
} // ~BwtReaderRunLengthBase::readAndCount( LetterCount& c, const LetterNumber numChars )


LetterNumber BwtReaderRunLengthBase::readAndSend( BwtWriterBase &writer, const LetterNumber numChars )
{
#ifdef DEBUG
    std::cout << "BR RL readAndSend " << numChars << " chars " << endl;
#endif
    LetterNumber charsLeft( numChars );
    while ( charsLeft > runLength_ )
    {
        //      int fred(whichPile[lastChar_]);
        writer.sendRun( lastChar_, runLength_ );
        //      c.count_[whichPile[lastChar_]]+=runLength_;
        charsLeft -= runLength_;
        if ( getRun() == false )
        {
            currentPos_ += ( numChars - charsLeft );
            return ( numChars - charsLeft );
            // assert(1==0);
        } // ~if
    } // ~while

    writer.sendRun( lastChar_, charsLeft );
    //    c.count_[whichPile[lastChar_]]+=charsLeft;
    runLength_ -= charsLeft;
    currentPos_ += numChars;
    return numChars;
} //~BwtReaderRunLengthBase::readAndSend(BwtWriterBase& writer, const LetterNumber numChars)

LetterNumber BwtReaderRunLengthBase::operator()( char *p, LetterNumber numChars )
{
#ifdef DEBUG
    std::cout << "BR RL () :  asked for " << numChars << " " << lastChar_ << " "
              << runLength_ << " " << pFile_ << std::endl;
#endif
    LetterNumber charsLeft( numChars );
    //    return fread( p, sizeof(char), numChars, pFile_ );
    while ( charsLeft > runLength_ )
    {
#ifdef DEBUG
        std::cout << "BR RL () :  setting " << lastChar_ << " "
                  << runLength_ << " " << pFile_ << std::endl;
#endif

        memset( p, lastChar_, runLength_ );
        p += runLength_;

        charsLeft -= runLength_;
        if ( getRun() == false )
        {
            // runLength_=0;
#ifdef DEBUG
            std::cout << "B read " << numChars - charsLeft << " out of "
                      << numChars << std::endl;
#endif
            currentPos_ += ( numChars - charsLeft );
            return ( numChars - charsLeft );
        } // ~if
    } // ~while
#ifdef DEBUG
    std::cout << "BR RL () :  last try - setting " << lastChar_ << " "
              << runLength_ << " " << pFile_ << std::endl;
#endif
    //     runLength_=lengths_[lastChar_];
    //  lastChar_=codes_[lastChar_];
#ifdef DEBUG
    std::cout << "BR RL () :  last try - setting " << lastChar_ << " "
              << runLength_ << " " << pFile_ << std::endl;
#endif

    memset( p, lastChar_, charsLeft );

    runLength_ -= charsLeft;
#ifdef DEBUG
    std::cout << "B delivered " << numChars << " " << charsLeft << " "
              << pFile_ << std::endl;
#endif
    currentPos_ += numChars;
    return numChars;
} // ~operator()



//
// BwtReaderRunLength member function definitions
//

BwtReaderRunLength::BwtReaderRunLength( const string &filename ):
    BwtReaderRunLengthBase( filename )
{
} // ~ctor

BwtReaderRunLength::BwtReaderRunLength( const BwtReaderRunLength &obj ):
    BwtReaderRunLengthBase( obj )
{
} // ~ctor

bool BwtReaderRunLength::getRun( void )
{
    if ( pBuf_ == pBufMax_ )
    {
        if ( finished_ )
        {
            runLength_ = 0;
            return false;
        }
        else
        {
            LetterNumber numRead( fread( buf_.data(), sizeof( uchar ), ReadBufferSize, pFile_ ) );
            if ( numRead == 0 )
            {
                runLength_ = 0;
                return false;
            }
            else if ( numRead < ReadBufferSize )
            {
                finished_ = true;
                pBufMax_ = buf_.data() + numRead;
            }
            pBuf_ = buf_.data();
        } // ~else
    } // ~if
    runLength_ = lengths_[*pBuf_];
    lastChar_ = codes_[*pBuf_];
#ifdef DEBUG
    cout << "Got run: " << runLength_ << " of " << lastChar_ << endl;
#endif
    ++pBuf_;
    ++currentPosInFile_;

    return true;

} // ~getRun



//
// BwtReaderRunLengthV3 member function definitions
//

BwtReaderRunLengthV3::BwtReaderRunLengthV3( const string &filename ):
    BwtReaderRunLengthBase( filename ),
    symbolForRunLength1ForPile_( 0 ),
    maxEncodedRunLengthForPile_( 0 ),
    firstContinuationSymbol_( 0 ),
    maxEncodedRunLengthMultiplierForContinuationSymbol_( 0 ),
    firstDataByteInFile_( 0 ),
    prefetchedByte_( 0 )
{
    // Check file header
    uchar c;
    c = fgetc( pFile_ );
    assert( c == 'B' );
    c = fgetc( pFile_ );
    assert( c == 'W' );
    c = fgetc( pFile_ );
    assert( c == 'T' );
    c = fgetc( pFile_ ); // \r\n sequence to check for invalid dos/unix format conversions
    assert( c == 13 );
    c = fgetc( pFile_ );
    assert( c == 10 );
    c = fgetc( pFile_ ); // Ctrl-Z, making some text viewers stop here and being non-ASCII to avoid confusion with text files
    assert( c == 26 );

    // Format version number, on 2 bytes to help identify endianness problems
    uint16_t formatVersion;
    assert( fread( &formatVersion, sizeof( formatVersion ), 1, pFile_ ) == 1 );
    assert( formatVersion == 3 );


    symbolForRunLength1ForPile_.resize( alphabetSize );
    maxEncodedRunLengthForPile_.resize( alphabetSize );
    firstContinuationSymbol_ = 0;
    maxEncodedRunLengthMultiplierForContinuationSymbol_ = -1;

    uint8_t base, rangeLength;
    uint16_t firstRunLength;
    uint16_t firstBytecode = 0;

    while (firstBytecode < 256 )
    {
        assert( fread( &base, sizeof( base ), 1, pFile_ ) == 1 );
        assert( fread( &rangeLength, sizeof( rangeLength ), 1, pFile_ ) == 1 );
        assert( fread( &firstRunLength, sizeof( firstRunLength ), 1, pFile_ ) == 1 );

        for (int i=0; i<rangeLength; ++i)
        {
            lengths_[ firstBytecode+i ] = firstRunLength+i;
            codes_[ firstBytecode+i ] = base;
        }

        if (base != '+')
        {
            assert( firstRunLength == 1 ); // current restriction
            symbolForRunLength1ForPile_[ whichPile[base] ] = firstBytecode;
            maxEncodedRunLengthForPile_[ whichPile[base] ] = firstRunLength + rangeLength - 1;
        }
        else
        {
            assert( firstRunLength == 0 ); // current restriction
            firstContinuationSymbol_ = firstBytecode;
            maxEncodedRunLengthMultiplierForContinuationSymbol_ = firstRunLength + rangeLength - 1;
        }

        firstBytecode += rangeLength;
    }
    assert (firstBytecode == 256 );

    currentPosInFile_ = firstDataByteInFile_ = ftell( pFile_ );
    prefetchNextByte();
} // ~ctor

BwtReaderRunLengthV3::BwtReaderRunLengthV3( const BwtReaderRunLengthV3 &obj ):
    BwtReaderRunLengthBase( obj ),
    symbolForRunLength1ForPile_( obj.symbolForRunLength1ForPile_ ),
    maxEncodedRunLengthForPile_( obj.maxEncodedRunLengthForPile_ ),
    firstContinuationSymbol_( obj.firstContinuationSymbol_ ),
    maxEncodedRunLengthMultiplierForContinuationSymbol_( obj.maxEncodedRunLengthMultiplierForContinuationSymbol_ ),
    firstDataByteInFile_( obj.firstDataByteInFile_ ),
    prefetchedByte_( obj.prefetchedByte_ )
{
} // ~ctor

bool BwtReaderRunLengthV3::getRun( void )
{
    if ( prefetchedByte_ == EOF )
    {
        lastChar_ = notInAlphabet;
        runLength_ = 0;
        return false;
    }

    LetterNumber currentContinuationMultiplier = 0;
    do
    {
        if (codes_[prefetchedByte_] != '+') // continuation code
        {
            lastChar_ = codes_[prefetchedByte_];
            runLength_ = lengths_[prefetchedByte_];
            currentContinuationMultiplier = maxEncodedRunLengthForPile_[ whichPile[lastChar_] ];
        }
        else
        {
            assert( currentContinuationMultiplier != 0 );
            runLength_ += lengths_[prefetchedByte_] * currentContinuationMultiplier;
            currentContinuationMultiplier *= maxEncodedRunLengthMultiplierForContinuationSymbol_ + 1;
        }

        ++currentPosInFile_;
        prefetchNextByte();
    }
    while (prefetchedByte_ != EOF && codes_[prefetchedByte_] == '+'); // continuation code

#ifdef DEBUG
    cout << "Got run: " << runLength_ << " of " << lastChar_ << endl;
#endif

    return true;
}

void BwtReaderRunLengthV3::prefetchNextByte()
{
    if (prefetchedByte_ == EOF)
        return;

    if ( pBuf_ == pBufMax_ )
    {
        LetterNumber numRead = fread( buf_.data(), sizeof( uchar ), ReadBufferSize, pFile_ );
        if ( numRead < ReadBufferSize )
        {
            pBufMax_ = buf_.data() + numRead;

            if ( numRead == 0 )
            {
                prefetchedByte_ = EOF;
                return;
            }
        }
        pBuf_ = buf_.data();
    }

    prefetchedByte_ = *pBuf_;
    ++pBuf_;
//    ++currentPosInFile_;
}

void BwtReaderRunLengthV3::rewindFile( void )
{
    seek( firstDataByteInFile_, 0 );
}

LetterNumber BwtReaderRunLengthV3::tellg( void ) const
{
    return currentPosInFile_;
}

int BwtReaderRunLengthV3::seek( const LetterNumber posInFile, const LetterNumber baseNumber )
{
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "seek filename=" << filename_ << " posInFile=" << posInFile << " baseNumber=" << baseNumber << endl;

    int ret = BwtReaderRunLengthBase::seek( posInFile, baseNumber );
    prefetchedByte_ = 0;
    prefetchNextByte();
    return ret;
}



// BwtReaderIncrementalRunLength member function definitions
//

extern vector< vector<unsigned char> > ramFiles; // declared in BwtWriter; todo: move those to another a new header file

BwtReaderIncrementalRunLength::BwtReaderIncrementalRunLength( const string &filename ):
    BwtReaderBase( filename ), runLength_( 0 ),
    pBuf_( buf_.data() + ReadBufferSize ), pBufMax_( buf_.data() + ReadBufferSize ),
    lastChar_( notInAlphabet ),
    lastMetadata_( 0 ),
    finished_( false ),
    currentPos_ ( 0 ),
    posInRamFile_ ( 0 )
{
    //    cout << "BwtReaderIncrementalRunLength: Opening " << filename << endl;
    if ( fread( &fileNum_, sizeof( fileNum_ ), 1, pFile_ ) != 1 )
    {
        fileNum_ = -1;
        finished_ = true;
    }
    //    cout << "  = file #" << fileNum_ << endl;
    fileNum_ %= 5;
    //    cout << "   => file #" << fileNum_ << endl;
    assert( ( int )ramFiles.size() > fileNum_ );

    for ( uint i( 0 ); i < 256; i++ )
    {
        lengths_[i] = 1 + ( i >> 4 );
        uint j = ( i & 0xF );
        codes_[i] = ( j < alphabetSize ) ? alphabet[j] : notInAlphabet;
    } // ~for i
} // ~ctor

void BwtReaderIncrementalRunLength::rewindFile( void )
{
    // rewind file and set all vars as per constructor
    rewind( pFile_ );
    runLength_ = 0;
    pBuf_ = buf_.data() + ReadBufferSize;
    pBufMax_ = buf_.data() + ReadBufferSize;
    lastChar_ = notInAlphabet;
    currentPos_ = 0;
    finished_ = false;
} // ~rewindFile

LetterNumber BwtReaderIncrementalRunLength::tellg( void ) const
{
    return currentPos_;
} // ~tellg

LetterNumber BwtReaderIncrementalRunLength::readAndCount( LetterCount &c, const LetterNumber numChars )
{
#ifdef DEBUG
    std::cout << "BR RL readAndCount " << numChars << " chars " << endl;
#endif
    LetterNumber charsLeft( numChars );
    while ( charsLeft > runLength_ )
    {
        // Below is not great design, at first call of this function it accesses an
        // out-of-range array element. Fortunately it always adds zero to it! :)
        c.count_[whichPile[lastChar_]] += runLength_;
        charsLeft -= runLength_;
        if ( getRun() == false )
        {
            currentPos_ += ( numChars - charsLeft );
            return ( numChars - charsLeft );
            //    assert(1==0);
        } // ~if
    } // ~while

    c.count_[whichPile[lastChar_]] += charsLeft;
    runLength_ -= charsLeft;
    currentPos_ += numChars;
    return numChars;
} // ~BwtReaderIncrementalRunLength::readAndCount( LetterCount& c, const LetterNumber numChars )

LetterNumber BwtReaderIncrementalRunLength::readAndSend( BwtWriterBase &writer, const LetterNumber numChars )
{
#ifdef DEBUG
    std::cout << "BR RL readAndSend " << numChars << " chars " << endl;
#endif
    bool isWriterIncremental = writer.isIncremental();
    LetterNumber charsLeft( numChars );
    while ( charsLeft > runLength_ )
    {
        if ( !isWriterIncremental )
            writer.sendRun( lastChar_, runLength_ );
        charsLeft -= runLength_;
        if ( getRun() == false )
        {
            currentPos_ += ( numChars - charsLeft );
            return ( numChars - charsLeft );
            //    assert(1==0);
        } // ~if
    } // ~while

    runLength_ -= charsLeft;
    currentPos_ += numChars;
    if ( !isWriterIncremental )
        writer.sendRun( lastChar_, charsLeft );
    else
        writer.sendRunOfPreExistingData( lastChar_, charsLeft, fileNum_, posInRamFile_, runLength_ );
    return numChars;
} //~BwtReaderIncrementalRunLength::readAndSend(BwtWriterBase& writer, const LetterNumber numChars)

LetterNumber BwtReaderIncrementalRunLength::operator()( char *p, LetterNumber numChars )
{
#ifdef DEBUG
    std::cout << "BR RL () :  asked for " << numChars << " " << lastChar_ << " "
              << runLength_ << " " << pFile_ << std::endl;
#endif
    LetterNumber charsLeft( numChars );
    //    return fread( p, sizeof(char), numChars, pFile_ );
    while ( charsLeft > runLength_ )
    {
#ifdef DEBUG
        std::cout << "BR RL () :  setting " << lastChar_ << " "
                  << runLength_ << " " << pFile_ << std::endl;
#endif

        memset( p, lastChar_, runLength_ );
        p += runLength_;

        charsLeft -= runLength_;
        if ( getRun() == false )
        {
            //    runLength_=0;
#ifdef DEBUG
            std::cout << "B read " << numChars - charsLeft << " out of "
                      << numChars << std::endl;
#endif
            currentPos_ += ( numChars - charsLeft );
            return ( numChars - charsLeft );
        } // ~if
    } // ~while
#ifdef DEBUG
    std::cout << "BR RL () :  last try - setting " << lastChar_ << " "
              << runLength_ << " " << pFile_ << std::endl;
#endif
    //     runLength_=lengths_[lastChar_];
    //  lastChar_=codes_[lastChar_];
#ifdef DEBUG
    std::cout << "BR RL () :  last try - setting " << lastChar_ << " "
              << runLength_ << " " << pFile_ << std::endl;
#endif

    memset( p, lastChar_, charsLeft );

    runLength_ -= charsLeft;
#ifdef DEBUG
    std::cout << "B delivered " << numChars << " " << charsLeft << " "
              << pFile_ << std::endl;
#endif
    currentPos_ += numChars;
    return numChars;
} // ~operator()

bool BwtReaderIncrementalRunLength::getRun( void )
{
again:
    if ( lastMetadata_ != 0 )
    {
        if ( ( lastMetadata_ & ~0x80 ) != 0 ) // if there is a subcall, whether or not the Return flag is also present
        {
            // push current file onto stack
            // we need a current position associated to each file
            int calledFileLevel = lastMetadata_ & ~0x80;
            assert( calledFileLevel < 128 );
            int oldFileNum = fileNum_;
            int newFileNum = 5 * calledFileLevel + ( fileNum_ % 5 );
            if ( ( int )posInRamFiles_.size() <= newFileNum )
            {
                posInRamFiles_.resize( newFileNum + 1 );
            }
            posInRamFiles_[oldFileNum] = posInRamFile_;
            posInRamFile_ = posInRamFiles_[newFileNum];
            stackedFileNums_.push_back( oldFileNum );
            fileNum_ = newFileNum;
        }
        else
        {
            // return bit (only) is set
            assert( lastMetadata_ == 0x80 );
            while ( ( lastMetadata_ & 0x80 ) != 0 && !stackedFileNums_.empty() )
            {
                // return from subfile
                assert( !stackedFileNums_.empty() );
                posInRamFiles_[fileNum_] = posInRamFile_;
                fileNum_ = stackedFileNums_.back();
                stackedFileNums_.pop_back();
                posInRamFile_ = posInRamFiles_[fileNum_];

                // retrieve the caller's lastMetadata to check its Return flag (the subcall may have been changed, but the return flag should persist)
                lastMetadata_ = ramFiles[fileNum_][posInRamFile_ - 1];
            }
        }
    }
    if ( finished_ || fileNum_ == -1 || posInRamFile_ >= ramFiles[fileNum_].size() )
    {
        finished_ = true;
        runLength_ = 0;
        return false;
    }
    const unsigned char c = ramFiles[fileNum_][posInRamFile_];
    lastMetadata_ = ramFiles[fileNum_][posInRamFile_ + 1];
    if ( posInRamFile_ + 2 >= ramFiles[fileNum_].size() )
    {
        // Add a Return bit if we reach the end of the file, for easier processing afterwards
        if ( fileNum_ > 5 )
        {
            lastMetadata_ |= 0x80;
        }
    }
#ifdef READ_DATA_FROM_FILES_FOR_DEBUGGING
    if ( pBuf_ == pBufMax_ )
    {
        if ( finished_ )
        {
            runLength_ = 0;
            return false;
        }
        else
        {
            LetterNumber numRead( fread( buf_.data(), sizeof( uchar ),
                                         ReadBufferSize, pFile_ ) );
            if ( numRead == 0 )
            {
                runLength_ = 0;
                return false;
            }
            else if ( numRead < ReadBufferSize )
            {
                //                finished_=true;
                pBufMax_ = buf_.data() + numRead;
            }
            pBuf_ = buf_.data();
        } // ~else
    } // ~if
    assert( c == *pBuf_ );
#endif //ifdef READ_DATA_FROM_FILES_FOR_DEBUGGING
    runLength_ = lengths_[( int )c];
    lastChar_ = codes_[( int )c];
#ifdef DEBUG
    cout << "Got run: " << runLength_ << " of " << lastChar_ << endl;
#endif
    pBuf_ += 2;
    posInRamFile_ += 2;

    if ( c == 0xFF )
        goto again;

    return true;

} // ~getRun

void BwtReaderIncrementalRunLength::defragment( void )
{
#define DEFRAGMENTATION_MAX_RUN_LENGTH 12
    vector<unsigned char> newRamFile;
    unsigned char prevLetter = 0;
    unsigned char prevRunLength = 0;
    assert ( runLength_ == 0 && "defragment shouldn't be called after any other operation" );

    while ( getRun() )
    {
        if ( prevRunLength == 0 )
        {
            prevLetter = whichPile[( int )lastChar_];
            prevRunLength = runLength_;//lastChar_ >> 4;
        }
        else
        {
            unsigned char newLetter =  whichPile[( int )lastChar_];
            unsigned char newRunLength = runLength_;
            if ( ( newLetter != prevLetter ) ) // || (prevRunLength >= DEFRAGMENTATION_MAX_RUN_LENGTH))
            {
                newRamFile.push_back( prevLetter | ( ( prevRunLength - 1 ) << 4 ) );
                newRamFile.push_back( 0 );
                prevLetter = newLetter;
                prevRunLength = newRunLength;
            }
            else
            {
                unsigned char totalRunLength = ( prevRunLength + newRunLength );
                if ( totalRunLength <= DEFRAGMENTATION_MAX_RUN_LENGTH )
                {
                    newRamFile.push_back( prevLetter | ( ( totalRunLength - 1 ) << 4 ) );
                    newRamFile.push_back( 0 );
                    prevRunLength = 0;
                }
                else
                {
                    prevRunLength = totalRunLength;
                    do
                    {
                        newRamFile.push_back( prevLetter | ( DEFRAGMENTATION_MAX_RUN_LENGTH - 1 ) << 4 );
                        newRamFile.push_back( 0 );
                        prevRunLength -= DEFRAGMENTATION_MAX_RUN_LENGTH;
                    }
                    while ( prevRunLength > DEFRAGMENTATION_MAX_RUN_LENGTH );
                }
            }
        }
    }

    if ( prevRunLength )
    {
        newRamFile.push_back( prevLetter | ( ( prevRunLength - 1 ) << 4 ) );
        newRamFile.push_back( 0 );
    }


    size_t sizeBefore = 0;
    for ( unsigned int i = fileNum_; i < ramFiles.size(); i += 5 )
    {
        vector<unsigned char> emptyVec;
        sizeBefore += ramFiles[i].size();
        //        ramFiles[i].clear();
        ramFiles[i].swap( emptyVec ); // deallocates vector memory
    }
    ramFiles[fileNum_].swap( newRamFile );
    size_t sizeAfter = ramFiles[fileNum_].size();
    Logger_if( LOG_SHOW_IF_VERBOSE )
    {
        Logger::out() << "defragment " << fileNum_ << " : size before= " << sizeBefore << " size after= " << sizeAfter << endl;
    }
}

//
// BwtReaderHuffman member function definitions
//

#ifdef ACTIVATE_HUFFMAN

BwtReaderHuffman::BwtReaderHuffman( const string &filename ):
    BwtReaderBase( filename ),
    runLength_( 0 ),
    lastChar_( notInAlphabet ),
    bitsUsed_( 0 ),
    finished_( false ),
    nearlyFinished_( false ),
    intCounter_( 0 ),
    numSymbols_( 0 ),
    maxSymbols_( 0 ),
    queueCounter_( 1 ), // needed for the first call of getRun()
    currentPos_( 0 )

{
    // just make sure everything is fine
    fseek( pFile_, 0, SEEK_END );
    long fileSize( ftell( pFile_ ) );
    // kill everything if file size is not a multiple of 4byte (~32 bit)
    assert( ( fileSize % sizeof( unsigned int ) ) == 0 ); // read with == 4 byte
    numInts_ = fileSize / sizeof( unsigned int ); // how many ints to read from file
    //cerr << filename << ": " << fileSize << " bytes/"<< numInts_ << " blocks" << endl;
    fseek( pFile_, 0, SEEK_SET );
    //init arrays

    for ( int i = 0; i < huffmanBufferSize; i++ )
    {
        symBuf[i] = 0;
        runBuf[i] = 0;
    }
    soFar_.ull = 0;
    toAdd_.ull = 0;

    // init the token lookup table, does not need to be cleared when the file
    // is rewind since its more or less static
    // TODO: hardcode the token table? stays the same for each program call

    unsigned int codeMask;
    for ( unsigned int i( 0 ); i < numTokens; i++ )
    {
        tokenTable_[i] = 0xFF;
        for ( unsigned int j( 0 ); j < numSingleCodes; j++ )
        {
            codeMask = ( 1 << singleCharLength[j] ) - 1; // (*2^3==) -1
            if ( ( i & codeMask ) == singleCharCode[j] )
            {
                assert ( tokenTable_[i] == 0xFF );
                tokenTable_[i] = ( j << 1 );
                //   cerr << "TT @ " << i << " is "<< itoa(tokenTable_[i],2) << endl;
            }
        } // ~for j
        for ( unsigned int j( 0 ); j < numDoubleCodes; j++ )
        {
            codeMask = ( 1 << doubleCharLength[j] ) - 1;
            if ( ( i & codeMask ) == doubleCharCode[j] )
            {
                assert ( tokenTable_[i] == 0xFF );
                tokenTable_[i] = ( ( j << 1 ) | 0x1 );
                //  cerr << "TT @ " << i << " is "<< itoa(tokenTable_[i],2) << endl;

            }
        } // ~for j
        //      assert (tokenTable_[i]!=0xFF); some tokens can have no prefix
        // that corresponds to a valid code
    } // ~for i

} // ~ctor

void BwtReaderHuffman::rewindFile( void )
{
    // rewind file and set all vars as per constructor
    rewind( pFile_ );
    runLength_ = 0;
    lastChar_ = notInAlphabet;
    currentPos_ = 0;
    finished_ = false;
    bitsUsed_ = 0;
    soFar_.ull = 0;
    numSymbols_ = 0;
    queueCounter_ = 0;
    maxSymbols_ = 0;
    intCounter_ = 0;
    firstRun_ = true; // fixes 3 Bit huffman error (BTL-17)
    nearlyFinished_ = false;

    for ( int i = 0; i < huffmanBufferSize; i++ )
    {
        symBuf[i] = 0;
        runBuf[i] = 0;
    }
} // ~rewindFile

LetterNumber BwtReaderHuffman::tellg( void ) const
{
    return currentPos_;
} // ~tellg

LetterNumber BwtReaderHuffman::readAndCount( LetterCount &c, const LetterNumber numChars )
{

    LetterNumber charsLeft( numChars );
    while ( charsLeft > runLength_ )
    {
        // Below is not great design, at first call of this function it accesses an
        // out-of-range array element. Fortunately it always adds zero to it! :)
        c.count_[whichPile[lastChar_]] += runLength_;
        charsLeft -= runLength_;
        if ( getRun() == false )
        {
            currentPos_ += ( numChars - charsLeft );
            return ( numChars - charsLeft );
        } // ~if
    } // ~while

    c.count_[whichPile[lastChar_]] += charsLeft;
    runLength_ -= charsLeft;
    currentPos_ += numChars;
    return numChars;
} // ~BwtReaderHuffman::readAndCount( LetterCount& c, const LetterNumber numChars )

LetterNumber BwtReaderHuffman::readAndSend( BwtWriterBase &writer, const LetterNumber numChars )
{
    if ( numChars == 0 )
    {
        return numChars;   // exit directy
    }

    LetterNumber charsLeft( numChars );
    while ( charsLeft > runLength_ )
    {

        writer.sendRun( lastChar_, runLength_ );
        charsLeft -= runLength_;

        if ( getRun() == false )
        {
            currentPos_ += ( numChars - charsLeft );
            return ( numChars - charsLeft );
        } // ~if
    } // ~while
    writer.sendRun( lastChar_, charsLeft );
    runLength_ -= charsLeft;
    currentPos_ += numChars;
    return numChars;
} //~BwtReaderHuffman::readAndSend(BwtWriterBase& writer, const LetterNumber numChars)


bool BwtReaderHuffman::getRun( void )
{
    numSymbols_ = -1; // needed for loops

    // no more data available AND or buffer is also empty -> finished here
    if ( finished_ && ( queueCounter_ > maxSymbols_ ) ) return false;


    // there is still data, we only have to fill out buffer
    if ( queueCounter_ > maxSymbols_ && !nearlyFinished_ )
    {

        unsigned int codeNum = 0;
        unsigned int codeSize = 0;
        unsigned int runLength = 0;
        unsigned int elementsRead = 0;

        toAdd_.ull = 0; // init ull AND both ui's with zeros

        elementsRead = fread( &toAdd_.ui, sizeof ( unsigned int ), 1, pFile_ );
        // try to read 32 bits in a row, if not return false
        if ( elementsRead == 1 )
        {
            toAdd_.ull <<= bitsUsed_; // left shift of how many bits used, start 0
            soFar_.ull |= toAdd_.ull; // glue both 32bit words
            bitsUsed_ += 32; // we have succesfully read 32 bit
            intCounter_++; // and we have used one int for that
        }


        if ( firstRun_ ) // first call, have to fill the integer buffer
        {
            toAdd_.ull = 0; // init
            elementsRead = fread( &toAdd_.ui, sizeof ( unsigned int ), 1, pFile_ );

            if ( elementsRead == 1 )
            {
                toAdd_.ull <<= bitsUsed_; // left shift of how many bits used, start 0
                soFar_.ull |= toAdd_.ull; // glue both 32bit words together
                bitsUsed_ += 32; // we have succesfully read 32 bit
                intCounter_++; // and we have used one int for that
            }
            firstRun_ = false;
        }

        while ( bitsUsed_ > 32 ) // as long as we have some more bits than 1 int uses
        {

            codeNum = tokenTable_[soFar_.ui & tokenMask]; // get codenum

            if ( ( codeNum & 0x1 ) == 0 ) // single code
            {
                codeNum >>= 1;
                // we have just read the stop sign
                if ( codeNum == finalCharCode )
                {
                    nearlyFinished_ = true;
                    break;
                }
                codeSize = singleCharLength[codeNum]; // how large is the code
                soFar_.ull >>= codeSize; // shift by these bits
                bitsUsed_ -= codeSize; // substract bits
                runLength = 1; // single run only
                numSymbols_++; // new symbol in buffer
                symBuf[numSymbols_] = alphabet[codeNum]; // add to buffer
                runBuf[numSymbols_] = runLength;          // ....
            }// ~if
            else   // double code
            {
                codeNum >>= 1;
                codeSize = doubleCharLength[codeNum];
                soFar_.ull >>= codeSize;
                bitsUsed_ -= codeSize;
                runLength = getNum( intCounter_ );
                numSymbols_++;
                symBuf[numSymbols_] = alphabet[codeNum];
                runBuf[numSymbols_] = runLength;
            } // ~else.

        } // ~while

        if ( intCounter_ == ( numInts_ ) && !nearlyFinished_ )
        {

            while ( 1 )
            {
                codeNum = tokenTable_[soFar_.ui & tokenMask];
                if ( ( codeNum & 0x1 ) == 0 )
                {
                    codeNum >>= 1;
                    if ( codeNum == finalCharCode )
                    {
                        nearlyFinished_ = true;
                        break;
                    }
                    codeSize = singleCharLength[codeNum];
                    soFar_.ull >>= codeSize;
                    bitsUsed_ -= codeSize;
                    assert( bitsUsed_ > 0 );
                    runLength = 1;
                    numSymbols_++;
                    symBuf[numSymbols_] = alphabet[codeNum];
                    runBuf[numSymbols_] = runLength;
                }// ~if
                else
                {
                    codeNum >>= 1;
                    codeSize = doubleCharLength[codeNum];
                    soFar_.ull >>= codeSize;
                    bitsUsed_ -= codeSize;
                    assert( bitsUsed_ > 0 );
                    int i( 0 );
                    runLength = getNum( i );
                    numSymbols_++;
                    symBuf[numSymbols_] = alphabet[codeNum];
                    runBuf[numSymbols_] = runLength;
                } // ~else
            } // ~while
        } // ~if
        //

        maxSymbols_ = numSymbols_;
        queueCounter_ = 0; // reset
        numSymbols_ = -1; // reset for next run
    }

    if ( nearlyFinished_ && queueCounter_ > maxSymbols_ ) // GET THIS RIGHT
    {
        finished_ = true; // now we are really finished
        return false;
    }
    else
    {
        runLength_ = runBuf[queueCounter_];
        lastChar_ = symBuf[queueCounter_];
        queueCounter_++;
        return true;
    }

} // ~getRun

unsigned int BwtReaderHuffman::getNum( int &i )
{
    unsigned int n( soFar_.ui & 0xF ); // only last 4 bit -> 4 bit encodig the number
    soFar_.ull >>= 4; // process shift
    bitsUsed_ -= 4; // set counter
    assert( bitsUsed_ > 0 );

    if ( n != 0xF ) // that would be excacly 15
    {
        n++;
        n++;
    }
    else
    {
        n = 0;
        int bitShift( 0 );
        bool carryOn;
        do
        {
            carryOn = ( ( soFar_.ui & 0x80 ) != 0 ); // test if 1000 0000 bit is set

            n |= ( ( soFar_.ui & 0x7F ) << bitShift ); // extract last 7 bits, shit by 0 in first iteration

            bitShift += 7; // next iter will shift by 7
            soFar_.ull >>= 8; // shift in the next 8 bits from left to right
            bitsUsed_ -= 8; // we used 8 bit so far
            assert( bitsUsed_ > 0 );


            if ( ( carryOn == true ) && ( bitsUsed_ < 8 ) ) // true if
            {
                toAdd_.ull = 0;
                i++;
                assert( fread( &toAdd_.ui, sizeof ( unsigned int ), 1, pFile_ ) == 1 );
                toAdd_.ull <<= bitsUsed_;
                soFar_.ull |= toAdd_.ull;
                bitsUsed_ += 32;
            } // ~if
        }// ~while
        while ( carryOn == true );
        n += 17;
    } // ~else

#ifdef DEBUG
    cout << "getNum " << n << endl;
#endif

    return n;

} // ~getNum

// just deprecated code to make everything compile
LetterNumber BwtReaderHuffman::operator()( char *p, LetterNumber numChars )
{
    assert( 1 == 0 );
    return -1;
} // ~operator()

#endif //ifdef ACTIVATE_HUFFMAN



//
// BwtReaderRunLengthRam member function definitions
//

BwtReaderRunLengthRam::BwtReaderRunLengthRam( const string &filename ):
    BwtReaderBase( filename ),
    runLength_( 0 ),
    lastChar_( notInAlphabet ),
    currentPos_( 0 ),
    isClonedObject_( false )
{
    // Find file size
    fseek( pFile_, 0, SEEK_END );
    size_t fileSize( ftell( pFile_ ) );
    fseek( pFile_, 0, SEEK_SET );

    if ( fileSize )
    {
#ifdef DONT_USE_MMAP
        #pragma omp critical (IO)
        cerr << "Info: Using malloc'ed BWT" << endl;
        // Allocate enough RAM to contain the whole file
        fullFileBuf_ = ( char * ) malloc( fileSize );
        assert( fullFileBuf_ != 0 && "Not enough RAM to load BWT" );
        size_t ret = fread( fullFileBuf_, fileSize, 1, pFile_ );
        assert( ret == 1 );
#else
        #pragma omp critical (IO)
        cerr << "Info: Using mmap'ed BWT" << endl;
        int fd = open( filename.c_str(), O_RDONLY );
        assert( fd >= 0 );
        fullFileBuf_ = ( char * )mmap( NULL, fileSize, PROT_READ, MAP_SHARED /*| MAP_LOCKED | MAP_POPULATE*/, fd, 0 );
        if ( fullFileBuf_ == ( void * ) - 1 )
        {
            perror( "Error: Map failed" );
            assert( false );
            exit( -1 );
        }
        //#define ACTIVATE_LOCKING
#ifdef ACTIVATE_LOCKING
        if ( mlock( fullFileBuf_, fileSize ) != 0 )
        {
            ostringstream oss;
            oss << "Error: Mlock failed for file " << filename << " with code " << errno;
            #pragma omp critical (IO)
            perror( oss.str().c_str() );
        }
#endif // ACTIVATE_LOCKING
#endif // DONT_USE_MMAP
    }
    else
    {
        fullFileBuf_ = 0;
    }
    sizeOfFullFileBuf_ = fileSize;
    posInFullFileBuf_ = 0;

    // This file shouldn't be needed anymore
    fclose( pFile_ );
    pFile_ = NULL;

    for ( unsigned int i( 0 ); i < 256; i++ )
    {
        lengths_[i] = 1 + ( i >> 4 );
        uint j = ( i & 0xF );
        codes_[i] = ( j < alphabetSize ) ? alphabet[j] : notInAlphabet;
    } // ~for i
} // ~ctor

BwtReaderRunLengthRam::BwtReaderRunLengthRam( const BwtReaderRunLengthRam &obj ):
    BwtReaderBase( obj.filename_ ),
    runLength_( 0 ),
    lastChar_( notInAlphabet ),
    currentPos_( 0 ),
#ifndef DONT_USE_MMAP
    mmapLength_( 0 ),
#endif
    isClonedObject_( true )
{
    fullFileBuf_ = obj.fullFileBuf_;
    sizeOfFullFileBuf_ = obj.sizeOfFullFileBuf_;
    posInFullFileBuf_ = 0;

    // This file shouldn't be needed anymore
    fclose( pFile_ );
    pFile_ = NULL;

    for ( unsigned int i( 0 ); i < 256; i++ )
    {
        lengths_[i] = 1 + ( i >> 4 );
        uint j = ( i & 0xF );
        codes_[i] = ( j < alphabetSize ) ? alphabet[j] : notInAlphabet;
    }
} // ~copy ctor


BwtReaderRunLengthRam::~BwtReaderRunLengthRam()
{
    if ( !isClonedObject_ )
#ifdef DONT_USE_MMAP
        free( fullFileBuf_ );
#else
        if ( fullFileBuf_ )
            munmap( fullFileBuf_, mmapLength_ );
#endif
}

void BwtReaderRunLengthRam::rewindFile( void )
{
    // rewind file and set all vars as per constructor
    //    rewind( pFile_ );
    runLength_ = 0;
    lastChar_ = notInAlphabet;
    currentPos_ = 0;
    posInFullFileBuf_ = 0;
} // ~rewindFile

LetterNumber BwtReaderRunLengthRam::tellg( void ) const
{
    return currentPos_;
} // ~tellg


LetterNumber BwtReaderRunLengthRam::readAndCount( LetterCount &c, const LetterNumber numChars )
{
#ifdef DEBUG_RAC
    std::cout << "BR RL readAndCount " << numChars << " chars " << endl;
    std::cout << "Before: " << currentPos_ << " " << ftell( pFile_ ) << " " << c << endl;
#endif

    LetterNumber charsLeft( numChars );
    while ( charsLeft > runLength_ )
    {
        // Below is not great design, at first call of this function it accesses an
        // out-of-range array element. Fortunately it always adds zero to it! :)
        c.count_[whichPile[lastChar_]] += runLength_;
        charsLeft -= runLength_;
#ifdef DEBUG_RAC
        std::cout << "R&C: " << currentPos_ << " " << posInFullFileBuf_ << " " << charsLeft << " " << runLength_ << " " << lastChar_ << " " << c << endl;
#endif
        if ( getRun() == false )
        {
            currentPos_ += ( numChars - charsLeft );
#ifdef DEBUG_RAC
            std::cout << "After (end): " << currentPos_ << " " << posInFullFileBuf_ << " " << c << endl;
#endif
            return ( numChars - charsLeft );
            // assert(1==0);


        } // ~if
    } // ~while

    c.count_[whichPile[lastChar_]] += charsLeft;
    runLength_ -= charsLeft;
    currentPos_ += numChars;
#ifdef DEBUG_RAC
    std::cout << "After (not at end): " << currentPos_ << " " << posInFullFileBuf_ << " " << c << endl;
#endif

    return numChars;
} // ~BwtReaderRunLengthRam::readAndCount( LetterCount& c, const LetterNumber numChars )


LetterNumber BwtReaderRunLengthRam::readAndSend( BwtWriterBase &writer, const LetterNumber numChars )
{
#ifdef DEBUG
    std::cout << "BR RL readAndSend " << numChars << " chars " << endl;
#endif
    LetterNumber charsLeft( numChars );
    while ( charsLeft > runLength_ )
    {
        //      int fred(whichPile[lastChar_]);
        writer.sendRun( lastChar_, runLength_ );
        //      c.count_[whichPile[lastChar_]]+=runLength_;
        charsLeft -= runLength_;
        if ( getRun() == false )
        {
            currentPos_ += ( numChars - charsLeft );
            return ( numChars - charsLeft );
            // assert(1==0);
        } // ~if
    } // ~while

    writer.sendRun( lastChar_, charsLeft );
    //    c.count_[whichPile[lastChar_]]+=charsLeft;
    runLength_ -= charsLeft;
    currentPos_ += numChars;
    return numChars;
} //~BwtReaderRunLengthRam::readAndSend(BwtWriterBase& writer, const LetterNumber numChars)

LetterNumber BwtReaderRunLengthRam::operator()( char *p, LetterNumber numChars )
{
#ifdef DEBUG
    std::cout << "BR RL () :  asked for " << numChars << " " << lastChar_ << " "
              << runLength_ << " " << pFile_ << std::endl;
#endif
    LetterNumber charsLeft( numChars );
    //    return fread( p, sizeof(char), numChars, pFile_ );
    while ( charsLeft > runLength_ )
    {
#ifdef DEBUG
        std::cout << "BR RL () :  setting " << lastChar_ << " "
                  << runLength_ << " " << pFile_ << std::endl;
#endif

        memset( p, lastChar_, runLength_ );
        p += runLength_;

        charsLeft -= runLength_;
        if ( getRun() == false )
        {
            // runLength_=0;
#ifdef DEBUG
            std::cout << "B read " << numChars - charsLeft << " out of "
                      << numChars << std::endl;
#endif
            currentPos_ += ( numChars - charsLeft );
            return ( numChars - charsLeft );
        } // ~if
    } // ~while
#ifdef DEBUG
    std::cout << "BR RL () :  last try - setting " << lastChar_ << " "
              << runLength_ << " " << pFile_ << std::endl;
#endif
    //     runLength_=lengths_[lastChar_];
    //  lastChar_=codes_[lastChar_];
#ifdef DEBUG
    std::cout << "BR RL () :  last try - setting " << lastChar_ << " "
              << runLength_ << " " << pFile_ << std::endl;
#endif

    memset( p, lastChar_, charsLeft );

    runLength_ -= charsLeft;
#ifdef DEBUG
    std::cout << "B delivered " << numChars << " " << charsLeft << " "
              << pFile_ << std::endl;
#endif
    currentPos_ += numChars;
    return numChars;
} // ~operator()

bool BwtReaderRunLengthRam::getRun( void )
{
    if ( posInFullFileBuf_ == sizeOfFullFileBuf_ )
    {
        runLength_ = 0;
        //        cerr << "end reached at " << posInFullFileBuf_ << endl;
        return false;
    }
    else
    {
        const unsigned char c = fullFileBuf_[posInFullFileBuf_];
        runLength_ = lengths_[c];
        lastChar_ = codes_[c];
        ++posInFullFileBuf_;
#ifdef DEBUG
        cout << "Got run: " << runLength_ << " of " << lastChar_ << endl;
#endif
    }

    return true;

} // ~getRun


BwtReaderBase* instantiateBwtPileReader( const string &pileFilename, const string &useShm, const bool keepBwtInRam, const bool forceNotUseIndexClass )
{
    if ( keepBwtInRam )
    {
        Logger::out() << "Info: the option to keep BWT in RAM has been temporarily de-activated" << endl;
    }

    // Detect BWT file type (ASCII/RLE, and which RLE version)
    ifstream bwtFile( pileFilename.c_str() );
    vector<char> buf( 8, 'A' ); // initialised with 'A's to default to ASCII if not enough data
    bwtFile.read( buf.data(), 8 );

    if ( equal( buf.begin(), buf.end(), rleV3Header.begin() ) )
    {
        // RLE_v3 detected
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "BWT file " << pileFilename << " detected as RLE version 3" << endl;
        if ( readWriteCheck( (pileFilename+".idx").c_str(), false, false ) && !forceNotUseIndexClass )
        {
            Logger::out() << "Using indexed BWT file for " << pileFilename << endl;
            return new BwtReaderIndex<BwtReaderRunLengthV3>( pileFilename, useShm );
        }
        else
        {
            return new BwtReaderRunLengthV3( pileFilename );
        }
    }
    else
    {
        bool allAscii = true;
        for ( const char c: buf )
        { 
            switch ( toupper( c ) )
            {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'N':
            case '$':
                break;
            default:
                allAscii = false;
            }
        }

        if (allAscii)
        {
            // ASCII detected
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "BWT file " << pileFilename << " detected as ASCII" << endl;
            return new BwtReaderASCII( pileFilename );
        }
        else
        {
            // old RLE detected
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "BWT file " << pileFilename << " detected as RLE" << endl;
            if ( readWriteCheck( (pileFilename+".idx").c_str(), false, false ) && !forceNotUseIndexClass )
            {
                Logger::out() << "Using indexed BWT file for " << pileFilename << endl;
                return new BwtReaderIndex<BwtReaderRunLength>( pileFilename, useShm );
            }
            else
            {
                return new BwtReaderRunLength( pileFilename );
            }
        }
    }

    assert( false && "unreachable code" );
    return NULL;
}

vector <BwtReaderBase *> instantiateBwtPileReaders( const string &bwtPrefix, const string &useShm )
{
    vector <BwtReaderBase *> inBwt( alphabetSize );

    for ( int i( 0 ); i < alphabetSize; i++ )
    {
        stringstream filenameSS;
        filenameSS << bwtPrefix << "-B0" << i;
        string pileFilename = filenameSS.str().c_str();
        inBwt[i] = instantiateBwtPileReader( pileFilename, useShm, false );
    }

    return inBwt;
}
