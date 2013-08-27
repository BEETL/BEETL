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

#include "BwtReader.hh"

#include "BwtWriter.hh"
#include "CountWords.hh"
#include "LetterCount.hh"

#include <cmath>
#include <cstdlib>
#include <cstring>

using namespace std;


//#define DEBUG_RAC 1

//
// BwtReaderBase member function definitions
//

BwtReaderBase::BwtReaderBase( const string &fileName ) :
    pFile_( fopen( fileName.c_str(), "r" ) )
{
    if ( pFile_ == NULL )
    {
        cerr << "!! BwtReaderBase: failed to open file " << fileName << endl;
        exit( EXIT_FAILURE );
    }
#ifdef DEBUG
    cout << "BwtReaderBase: opened file " << fileName << " " << pFile_ << endl;
#endif
}

BwtReaderBase::~BwtReaderBase()
{
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



LetterNumber BwtReaderASCII::readAndCount( LetterCount &c, const LetterNumber numChars )
{
#ifdef DEBUG
    std::cout << "BR ASCII readAndCount " << numChars << " chars " << endl;
#endif
    LetterNumber charsLeft( numChars ), charsToRead, charsRead;
    while ( charsLeft > 0 )
    {
        charsToRead = ( ( charsLeft > ReadBufferSize ) ? ReadBufferSize : charsLeft );
        charsRead = fread( buf_, sizeof( char ), charsToRead, pFile_ );
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
        numRead = fread( buf_, sizeof ( uchar ), readNextPass, pFile_ );
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
// BwtReaderRunLength member function definitions
//

BwtReaderRunLength::BwtReaderRunLength( const string &fileName ):
    BwtReaderBase( fileName ), runLength_( 0 ),
    pBuf_( buf_ + ReadBufferSize ), pBufMax_( buf_ + ReadBufferSize ),
    lastChar_( notInAlphabet ),
    finished_( false ),
    currentPos_( 0 )
{
    unsigned int j;
    for ( unsigned int i( 0 ); i < 256; i++ )
    {
        lengths_[i] = 1 + ( i >> 4 );
        j = ( i & 0xF );
        codes_[i] = ( j < alphabetSize ) ? alphabet[j] : notInAlphabet;
    } // ~for i
} // ~ctor

void BwtReaderRunLength::rewindFile( void )
{
    // rewind file and set all vars as per constructor
    rewind( pFile_ );
    runLength_ = 0;
    pBuf_ = buf_ + ReadBufferSize;
    pBufMax_ = buf_ + ReadBufferSize;
    lastChar_ = notInAlphabet;
    currentPos_ = 0;
    finished_ = false;
} // ~rewindFile

LetterNumber BwtReaderRunLength::tellg( void ) const
{
    return currentPos_;
} // ~tellg


LetterNumber BwtReaderRunLength::readAndCount( LetterCount &c, const LetterNumber numChars )
{
#ifdef DEBUG_RAC
    std::cout << "BR RL readAndCount " << numChars << " chars " << endl;
    std::cout << "Before: " << currentPos_ << " " << ftell( pFile_ ) << " ";
    c.print();
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
#ifdef DEBUG_RAC
            std::cout << "After (end): " << currentPos_ << " " << ftell( pFile_ ) << " ";
            c.print();
#endif
            return ( numChars - charsLeft );
            // assert(1==0);


        } // ~if
    } // ~while

    c.count_[whichPile[lastChar_]] += charsLeft;
    runLength_ -= charsLeft;
    currentPos_ += numChars;
#ifdef DEBUG_RAC
    std::cout << "After (not at end): " << currentPos_ << " " << ftell( pFile_ ) << " ";
    c.print();
#endif

    return numChars;
} // ~BwtReaderRunLength::readAndCount( LetterCount& c, const LetterNumber numChars )


LetterNumber BwtReaderRunLength::readAndSend( BwtWriterBase &writer, const LetterNumber numChars )
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
} //~BwtReaderRunLength::readAndSend(BwtWriterBase& writer, const LetterNumber numChars)

LetterNumber BwtReaderRunLength::operator()( char *p, LetterNumber numChars )
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
            LetterNumber numRead( fread( buf_, sizeof( uchar ),
                                         ReadBufferSize, pFile_ ) );
            if ( numRead == 0 )
            {
                runLength_ = 0;
                return false;
            }
            else if ( numRead < ReadBufferSize )
            {
                finished_ = true;
                pBufMax_ = buf_ + numRead;
            }
            pBuf_ = buf_;
        } // ~else
    } // ~if
    runLength_ = lengths_[*pBuf_];
    lastChar_ = codes_[*pBuf_];
#ifdef DEBUG
    cout << "Got run: " << runLength_ << " of " << lastChar_ << endl;
#endif
    pBuf_++;

    return true;

} // ~getRun

//
// BwtReaderRunLengthIndex member function definitions
//

BwtReaderRunLengthIndex::BwtReaderRunLengthIndex( const string &fileName ):
    BwtReaderRunLength( fileName ),
    indexFileName_( fileName + ".idx" ),
    isNextIndex_( false ),
    pIndexFile_( NULL )
{
    current_.clear();
    initIndex( current_ );
}


void BwtReaderRunLengthIndex::rewindFile( void )
{
    // rewind file and set all vars as per constructor
    current_.clear();
    initIndex( current_ );
    BwtReaderRunLength::rewindFile();
} // ~rewindFile


LetterNumber BwtReaderRunLengthIndex::readAndCount( LetterCount &c, const LetterNumber numChars )
{
#ifdef DEBUG_RAC
    std::cout << "BR RLI readAndCount " << numChars << " chars " << endl;
    std::cout << "Before: " << currentPos_ << " " << ftell( pFile_ ) << " ";
    c.print();
    std::cout << "Internal start: ";
    current_.print();
#endif
    // bool needToShift(false);
    LetterNumber charsLeft( numChars );

    temp_ = current_;

    if ( ( isNextIndex_ ) && ( charsLeft > nextPos_ - currentPos_ ) )
    {
        //    temp_=currentIndexPos_; // squirrel away current counts vector for final subtraction

        //    currentPos_=currentIndexPos_;
        //    current_=currentIndex;

        //  while ((isNextIndex_)&&(currentPos_+charsLeft>=nextPos_))
        charsLeft += ( currentPos_ - currentIndexPos_ ); // gone back so need to read more
        while ( ( isNextIndex_ ) && ( charsLeft >= nextPos_ - currentIndexPos_ ) )
        {
#ifdef DEBUG_RAC
            std::cout << "BR RLI shift" << endl;
            std::cout << charsLeft << " " << currentPos_ << " " << nextPos_ << " " << currentIndexPos_ << std::endl;
            std::cout << "Shifting (RLI): " << currentPos_ << " " << ftell( pFile_ ) << " ";
            currentIndex_.print();
            next_.print();
#endif

            charsLeft -= ( nextPos_ - currentIndexPos_ );
            currentIndexPos_ = nextPos_;
            currentFilePos_ = nextFilePos_;
            currentIndex_ += next_;

            isNextIndex_ =
                ( fread( &nextPos_, sizeof( LetterNumber ), 1, pIndexFile_ ) == 1 );
            if ( isNextIndex_ )
            {
                assert
                ( fread( &nextFilePos_, sizeof( LetterNumber ), 1, pIndexFile_ ) == 1 );
                assert
                ( fread( &next_, sizeof( LetterCount ), 1, pIndexFile_ ) == 1 );
#ifdef DEBUG_RAC
                std::cout << "BR RLI readAndCount " << numChars << " chars " << endl;
                std::cout << nextPos_ << " " << nextFilePos_ << " ";
                next_.print();
#endif


            } // ~if

        }

        currentPos_ = currentIndexPos_;
        current_ = currentIndex_;
        fseek( pFile_, currentFilePos_, SEEK_SET );
        runLength_ = 0;
        pBuf_ = buf_ + ReadBufferSize;
        pBufMax_ = buf_ + ReadBufferSize;
        //    lastChar_=notInAlphabet;
    } // ~if

#ifdef DEBUG_RAC
    std::cout << "Shifted (RLI): " << currentPos_ << " " << ftell( pFile_ ) << " ";
    current_.print();
#endif


    charsLeft = BwtReaderRunLength::readAndCount( current_, charsLeft );
    c += current_;
    c -= temp_;

#ifdef DEBUG_RAC
    std::cout << "After (RLI) internal: " << currentPos_ << " " << ftell( pFile_ ) << " ";
    current_.print();
    std::cout << "After (RLI) returned: ";
    c.print();
#endif


    return charsLeft;
} // ~BwtReaderRunLengthIndex::readAndCount( LetterCount& c, const LetterNumber numChars )




void BwtReaderRunLengthIndex::buildIndex
( FILE *pIndexFile, const int indexBinSize )
{
    cout << "buildIndex" << endl;
    const int runsPerChunk( indexBinSize );
    int runsThisChunk( 0 );
    LetterCount countsThisChunk;
    LetterNumber runsSoFar( 0 );
    currentPos_ = 0;


    while ( 1 )
    {
        if ( getRun() == false ) break;
        runsSoFar++;
        runsThisChunk++;

        countsThisChunk.count_[whichPile[lastChar_]] += runLength_;
        currentPos_ += runLength_;
        if ( runsThisChunk == runsPerChunk )
        {
            cout << currentPos_ << " " << runsSoFar;
            countsThisChunk.print();

            assert
            ( fwrite( &currentPos_, sizeof( LetterNumber ), 1, pIndexFile ) == 1 );

            assert
            ( fwrite( &runsSoFar, sizeof( LetterNumber ), 1, pIndexFile ) == 1 );

            assert
            ( fwrite( &countsThisChunk, sizeof( LetterCount ), 1, pIndexFile ) == 1 );

            runsThisChunk = 0;
            countsThisChunk.clear();
        }


    }


} // ~buildIndex

void BwtReaderRunLengthIndex::initIndex( const LetterCount &current )
{
    current_.clear();
    currentIndex_.clear();
    currentIndexPos_ = 0;
    if ( pIndexFile_ != NULL ) fclose( pIndexFile_ );
    pIndexFile_ = fopen( indexFileName_.c_str(), "r" );

    if ( pIndexFile_ != NULL )
    {
        isNextIndex_ =
            ( fread( &nextPos_, sizeof( LetterNumber ), 1, pIndexFile_ ) == 1 );
        if ( isNextIndex_ )
        {
            assert
            ( fread( &nextFilePos_, sizeof( LetterNumber ), 1, pIndexFile_ ) == 1 );
            assert
            ( fread( &next_, sizeof( LetterCount ), 1, pIndexFile_ ) == 1 );
        } // ~if
    } // ~if
    //  rewindFile();
} // ~initIndex

// BwtReaderIncrementalRunLength member function definitions
//

extern vector< vector<unsigned char> > ramFiles; // declared in BwtWriter; todo: move those to another a new header file

BwtReaderIncrementalRunLength::BwtReaderIncrementalRunLength( const string &fileName ):
    BwtReaderBase( fileName ), runLength_( 0 ),
    pBuf_( buf_ + ReadBufferSize ), pBufMax_( buf_ + ReadBufferSize ),
    lastChar_( notInAlphabet ),
    lastMetadata_( 0 ),
    finished_( false ),
    currentPos_ ( 0 ),
    posInRamFile_ ( 0 )
{
    //    cout << "BwtReaderIncrementalRunLength: Opening " << fileName << endl;
    if ( fread( &fileNum_, sizeof( fileNum_ ), 1, pFile_ ) != 1 )
    {
        fileNum_ = -1;
        finished_ = true;
    }
    //    cout << "  = file #" << fileNum_ << endl;
    fileNum_ %= 5;
    //    cout << "   => file #" << fileNum_ << endl;
    assert( ( int )ramFiles.size() > fileNum_ );

    uint j;
    for ( uint i( 0 ); i < 256; i++ )
    {
        lengths_[i] = 1 + ( i >> 4 );
        j = ( i & 0xF );
        codes_[i] = ( j < alphabetSize ) ? alphabet[j] : notInAlphabet;
    } // ~for i
} // ~ctor

void BwtReaderIncrementalRunLength::rewindFile( void )
{
    // rewind file and set all vars as per constructor
    rewind( pFile_ );
    runLength_ = 0;
    pBuf_ = buf_ + ReadBufferSize;
    pBufMax_ = buf_ + ReadBufferSize;
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
            LetterNumber numRead( fread( buf_, sizeof( uchar ),
                                         ReadBufferSize, pFile_ ) );
            if ( numRead == 0 )
            {
                runLength_ = 0;
                return false;
            }
            else if ( numRead < ReadBufferSize )
            {
                //                finished_=true;
                pBufMax_ = buf_ + numRead;
            }
            pBuf_ = buf_;
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
    cout << "defragment " << fileNum_ << " : size before= " << sizeBefore << " size after= " << sizeAfter << endl;
}

//
// BwtReaderHuffman member function definitions
//

BwtReaderHuffman::BwtReaderHuffman( const string &fileName ):
    BwtReaderBase( fileName ),
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
    // just make sure everthing is fine
    fseek( pFile_, 0, SEEK_END );
    long fileSize( ftell( pFile_ ) );
    // kill everything if file size is not a multiple of 4byte (~32 bit)
    assert( ( fileSize % sizeof( unsigned int ) ) == 0 ); // read with == 4 byte
    numInts_ = fileSize / sizeof( unsigned int ); // how many ints to read from file
    //cerr << fileName << ": " << fileSize << " bytes/"<< numInts_ << " blocks" << endl;
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
                int i( 0 );
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
