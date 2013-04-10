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

#include "BwtWriter.hh"

#include "LetterCount.hh"
#include "Logger.hh"
#include "Tools.hh"

#include <cstdlib>
#include <cstring>

using namespace std;


//#define DEBUG 1
#define LOCAL_DEBUG 0

//
// BwtWriterBase member function definitions
//

BwtWriterFile::BwtWriterFile( const string &fileName ) : pFile_( fopen( fileName.c_str(), "w" ) )
{
#ifdef DEBUG
    cout << "BwtWriterFile opened file " << fileName << " " << pFile_ << endl;
#endif
    readWriteCheck( fileName.c_str(), 1 ); //    setvbuf( pFile_, NULL, _IOFBF, 262144);
}
BwtWriterFile::~BwtWriterFile()
{

    fclose( pFile_ );
#ifdef DEBUG
    cout << "BwtWriterFile: closed file " << pFile_ << endl;
#endif

}


//
// BwtWriterASCII member function definitions
//

BwtWriterASCII::BwtWriterASCII( const string &fileName ) : BwtWriterFile( fileName )
{
#ifdef DEBUG
    cout << "BW ASCII ctor" << endl;
#endif
}
BwtWriterASCII::~BwtWriterASCII()
{
#ifdef DEBUG
    cout << "BW ASCII dtor" << endl;
#endif

}


void BwtWriterASCII::operator()( const char *p, int numChars )
{
#ifdef DEBUG
    cout << "BW ASCII () - " << *p << " " << numChars << endl;
#endif


    size_t bytesWritten = fwrite( p, sizeof( char ), numChars, pFile_ );
    if ( bytesWritten != ( size_t )numChars )
    {
        cerr << "Unable to write " << numChars
             << " chars. Aborting." << endl;
        exit( EXIT_FAILURE );
    }
} // ~operator()

void BwtWriterASCII::sendRun( char c, int runLength )
{
    for ( int i( 0 ); i < runLength; i++ ) fprintf( pFile_, "%c", c );
}

//
// BwtWriterRunLength member function definitions
//

BwtWriterRunLength::~BwtWriterRunLength()
{
    if ( runLength_ != 0 ) encodeRun( lastChar_, runLength_ );


    if ( pBuf_ != buf_ )
    {

        size_t bytesWritten = fwrite( buf_, sizeof( char ), ( pBuf_ - buf_ ), pFile_ );
        if ( bytesWritten != ( size_t )( pBuf_ - buf_ ) )
        {
            cerr << "Unable to write " << ( pBuf_ - buf_ )
                 << " chars. Aborting." << endl;
            exit( EXIT_FAILURE );
        }

#ifdef REPORT_COMPRESSION_RATIO
        bytesWritten_ += ( LetterCountType )( pBuf_ - buf_ );
#endif
    }

#ifdef REPORT_COMPRESSION_RATIO
    Logger::out( LOG_FOR_DEBUGGING )
            << "BwtWriterRunLength: received "
            << charsReceived_ << " chars, sent "
            << bytesWritten_ << " bytes, compression "
            << ( ( double )8 * bytesWritten_ ) / ( charsReceived_ )
            << " bits per char " << std::endl;
#endif


}

void BwtWriterRunLength::operator()( const char *p, int numChars )
{
#ifdef DEBUG
    std::cout << "BW RL () - " << *p << " " << numChars << " state: " << lastChar_ << " " << runLength_ << endl;
#endif

    for ( int i( 0 ); i < numChars; i++ )
    {
        if ( ( *p ) == lastChar_ )
        {
            runLength_++;
        }
        else
        {
            if ( runLength_ > 0 )
            {
                encodeRun( lastChar_, runLength_ );
            } // ~if
            runLength_ = 1;
            lastChar_ = *p;
        } // ~else
        p++;
    } // ~for
    //    assert(fwrite( p, sizeof(char), numChars, pFile_ )==numChars);
} // ~operator()

void BwtWriterRunLength::sendChar( char c )
{
    *pBuf_ = c;
    if ( ++pBuf_ == pBufMax_ )
    {
#ifdef REPORT_COMPRESSION_RATIO
        bytesWritten_ += ReadBufferSize;
#endif

        size_t bytesWritten = fwrite( buf_, sizeof( char ), ReadBufferSize, pFile_ );
        if ( bytesWritten != ( size_t )ReadBufferSize )
        {
            cerr << "Unable to write " << ReadBufferSize
                 << " chars. Aborting." << endl;
            exit( EXIT_FAILURE );
        }
        pBuf_ = buf_;
    }
}

void BwtWriterRunLength::encodeRun( char c, unsigned int runLength )
{
#ifdef DEBUG
    std::cout << "BW RL encodeRun - sending run " << c << " " << runLength << " " << pFile_
              << std::endl;
#endif
#ifdef REPORT_COMPRESSION_RATIO
    charsReceived_ += runLength;
#endif
    int charIndex( whichPile[( int )c] );

    if ( charIndex == nv )
    {
        cerr << "Char |" << c << "| is not part of the alphabet. Aborting." << endl;
        exit( EXIT_FAILURE );
    }

    uchar outCode( 0xF0 | ( ( uchar )charIndex ) );
    runLength--;
    const unsigned int numMaxChars( runLength >> 4 );
    for ( unsigned int i( 0 ); i < numMaxChars; i++ )
    {
        sendChar( outCode );
    }
    runLength &= ( unsigned int )0xF;

    outCode = ( ( ( uchar )runLength ) << 4 );
    outCode |= charIndex;
    //    assert(((unsigned int)outCode)<256);
    //    assert(fwrite( &outCode, sizeof(char), 1, pFile_ )==1);
    sendChar( outCode );
#ifdef DEBUG
    std::cout << "B sending " << ( unsigned int )outCode << " " << pFile_ << std::endl;
#endif
} // ~encodeRun

void BwtWriterRunLength::sendRun( char c, int runLength )
{
#ifdef DEBUG
    std::cout << "BW RL sendRun - sending run " << c << " " << runLength << " " << endl;
#endif
    if ( runLength != 0 )
    {
        if ( c == lastChar_ )
        {
            runLength_ += runLength;
        }
        else
        {
            if ( runLength_ != 0 ) encodeRun( lastChar_, runLength_ );
            lastChar_ = c;
            runLength_ = runLength;
        }
    }
} // ~sendRun



#ifdef ORIG
void BwtWriterRunLength::sendRun( char c, int runLength )
{
#ifdef DEBUG
    std::cout << "BW RL sendRun - sending run " << c << " " << runLength << " " << endl;
#endif
    if ( c == lastChar_ )
    {
        runLength_ += runLength;
    }
    else
    {
        if ( runLength_ != 0 )
        {
            encodeRun( lastChar_, runLength_ );
            lastChar_ = c;
            runLength_ = runLength;
        }
    }
} // ~sendRun
#endif


//
// BwtWriterIncrementalRunLength member function definitions
//

uint BwtWriterIncrementalRunLength::nextFileNum_ = 0;
vector< vector<unsigned char> > ramFiles( 1000 ); //TODO: make this '1000' dynamic (the array resize just needs to be put in an openmp critical section)

BwtWriterIncrementalRunLength::BwtWriterIncrementalRunLength( const string &fileName )
    : BwtWriterFile( fileName )
    , runLength_( 0 ), pBuf_( buf_ ), pBufMax_( buf_ + ReadBufferSize ), lastChar_( notInAlphabet )
#ifdef REPORT_COMPRESSION_RATIO
    , charsReceived_( 0 ), bytesWritten_( 0 )
#endif
    , fileNumInReader_( 0 )
    , filePosInReader_( 0 )
    , remainingRunLengthInReader_( 0 )
    , lastFileReturnNeeded_( false )
    //    , onHoldUntilNextReturn_data_( 0 )
    , onHoldUntilNextReturn_letter_( 0 )
    , onHoldUntilNextReturn_runLength_( 0 )
    , onHoldUntilNextReturn_metadata_( 0 )
{
    //    assert( ramFiles.size() == nextFileNum_ );
    //    cout << "BwtWriterIncrementalRunLength: Opening " << fileName << endl;

    // todo: remove this hack, which is here to make sure the ramFiles%5 keep pointing to the correct alphabet BWT
    const char *ptr = fileName.c_str();
    while ( *ptr && ( *ptr < '0' || *ptr > '9' ) ) // skip any prefix like in "new_1"
        ++ptr;
    unsigned int letterNum = atoi( ptr );
    /*
        assert( letterNum >= 1 && letterNum <= 5 );
        do {
            fileNum_ = nextFileNum_++;
        }
        while ( (fileNum_%5) != (letterNum-1) );
    */
    extern unsigned int debugCycle;
    extern unsigned int lastDefragCycle;
    if ( lastDefragCycle == 0 )
        fileNum_ = letterNum - 1 + 5 * ( debugCycle - 1 );
    else
        fileNum_ = letterNum - 1 + 5 * ( debugCycle - lastDefragCycle + 1 );


    //    cout << "  = file #" << fileNum_ << endl;
    assert( fwrite( &fileNum_, sizeof( fileNum_ ), 1, pFile_ ) == 1 );
    /*
    #p.ragma omp critical
        if( fileNum_ >= ramFiles.size() )
        {
            ramFiles.resize( fileNum_ + 1 );
        }
    */
    //    ramFileLengths.resize( ramFiles.size() );
    fileNumInReader_ = fileNum_ % ( finalCharCode - 1 ); // todo: improve this modulo
    filePosInReader_ = 0;
}


BwtWriterIncrementalRunLength::~BwtWriterIncrementalRunLength()
{
    if ( runLength_ != 0 ) encodeRun( lastChar_, runLength_ );
    terminateLastInsertion();

    if ( pBuf_ != buf_ )
    {

        size_t bytesWritten = fwrite( buf_, sizeof( char ), ( pBuf_ - buf_ ), pFile_ );
        if ( bytesWritten != ( size_t )( pBuf_ - buf_ ) )
        {
            cerr << "Unable to write " << ( pBuf_ - buf_ )
                 << " chars. Aborting." << endl;
            exit( -1 );
        }

#ifdef REPORT_COMPRESSION_RATIO
        bytesWritten_ += ( LetterCountType )( pBuf_ - buf_ );
#endif
    }

#ifdef REPORT_COMPRESSION_RATIO

#ifndef SEND_DATA_TO_FILES_FOR_DEBUGGING
    bytesWritten_ = ramFiles[fileNum_].size();
#endif //ifndef SEND_DATA_TO_FILES_FOR_DEBUGGING

    Logger::out( LOG_FOR_DEBUGGING )
            << "BwtWriterIncrementalRunLength: received "
            << charsReceived_ << " chars, sent "
            << bytesWritten_ << " bytes, compression "
            << ( ( double )8 * bytesWritten_ ) / ( charsReceived_ )
            << " bits per char " << std::endl;
#endif


}

void BwtWriterIncrementalRunLength::operator()( const char *p, int numChars )
{
#ifdef DEBUG
    std::cout << "BW RL () - " << *p << " " << numChars << " state: " << lastChar_ << " " << runLength_ << endl;
#endif

    for ( int i( 0 ); i < numChars; i++ )
    {
        if ( ( *p ) == lastChar_ )
        {
            runLength_++;
        }
        else
        {
            if ( runLength_ > 0 )
            {
                encodeRun( lastChar_, runLength_ );
            } // ~if
            runLength_ = 1;
            lastChar_ = *p;
        } // ~else
        p++;
    } // ~for
} // ~operator()

void BwtWriterIncrementalRunLength::terminateLastInsertion()
{
    if ( !ramFiles[fileNum_].empty() )
    {
        size_t lastPos = ramFiles[fileNum_].size() - 1;
        unsigned char lastMetadata = ramFiles[fileNum_][lastPos];
        if ( lastMetadata )
        {
            assert( lastMetadata & 0x80 );
            assert( onHoldUntilNextReturn_letter_ == 0 );
            assert( onHoldUntilNextReturn_runLength_ == 0 );
            assert( onHoldUntilNextReturn_metadata_ == 0 );
            return;
        }

        if ( onHoldUntilNextReturn_metadata_ )
        {
            lastMetadata = onHoldUntilNextReturn_metadata_;
            onHoldUntilNextReturn_metadata_ = 0;
        }
        lastMetadata |= 0x80; // return bit

        if ( onHoldUntilNextReturn_letter_ || onHoldUntilNextReturn_runLength_ )
        {
            assert( onHoldUntilNextReturn_runLength_ > 0 );
            ramFiles[fileNum_].push_back( ( onHoldUntilNextReturn_runLength_ - 1 ) << 4 | onHoldUntilNextReturn_letter_ );
            ramFiles[fileNum_].push_back( lastMetadata );
            onHoldUntilNextReturn_letter_ = 0;
            onHoldUntilNextReturn_runLength_ = 0;
        }
        else
        {
            ramFiles[fileNum_][lastPos] = lastMetadata;
        }
    }
}

void BwtWriterIncrementalRunLength::sendChar( unsigned char c, unsigned char metadata )
{
    assert( metadata == 0 );
    uint fnum = fileNumInReader_;
    size_t &fpos = filePosInReader_;

    if ( LOCAL_DEBUG ) clog << "Inserting " << ( unsigned int )c << " in file " << fnum << " pos " << fpos << ", remainingRunLengthInReader=" << remainingRunLengthInReader_ << " fileNum=" << fileNum_ << " ramFiles[fileNum_].size()=" << ramFiles[fileNum_].size() << " ramFiles[fnum].size()=" << ramFiles[fnum].size() << endl;

    if ( lastFileReturnNeeded_ )
    {
        // New insertion starting => we terminate the last insertion properly
        terminateLastInsertion();
        lastFileReturnNeeded_ = false;
        assert( onHoldUntilNextReturn_letter_ == 0 );
        assert( onHoldUntilNextReturn_runLength_ == 0 );
        assert( onHoldUntilNextReturn_metadata_ == 0 );
    }

    // Continuation of an already started insertion
    if ( fileNum_ < 5 || ( !ramFiles[fileNum_].empty() && ( ramFiles[fileNum_].back() & 0x80 ) == 0 ) )
    {
        if ( LOCAL_DEBUG ) clog << " Continuation of an already started insertion" << endl;
        assert( ( fileNum_ < 5 || fpos != 0 ) && "todo: insertion at the start of a file" );
        ramFiles[fileNum_].push_back( c );
        ramFiles[fileNum_].push_back( metadata );
        return;
    }

    // Insertion at the start of a file
    if ( fpos == 0 )
    {
        assert( fnum < 5 && "Only first cycle files may be empty" );
        if ( LOCAL_DEBUG ) clog << " Insertion at the start of a file" << endl;
        unsigned char replacementMetadata = fileNum_ / 5;
        assert( ( replacementMetadata & 0x80 ) == 0 && "Error: there may be too many cycles for this algorithm, which overwrote the Return bit" );
        if ( ramFiles[fnum].empty() )
        {
            ramFiles[fnum].push_back( 0xFF ); // Special byte meaning that we'll ignore this base (but not its associated metadata)
            ramFiles[fnum].push_back( replacementMetadata );
            fpos = 2;

            ramFiles[fileNum_].push_back( c );
            ramFiles[fileNum_].push_back( metadata );
            return;
        }
        else if ( ramFiles[fnum][0] == 0xFF )
        {
            fpos = 2;
        }
        else
        {
            assert( false && "todo" );
        }

    }

    // Insertion of a letter that can be added to the run-length-encoded char
    unsigned char count1 = 1 + ( ramFiles[fnum][fpos - 2] >> 4 );
    unsigned char count2 = 1 + ( c >> 4 );
    unsigned char letter1 = ramFiles[fnum][fpos - 2] & 0x0F;
    unsigned char letter2 = c & 0x0F;
    if ( ( letter1 == letter2 ) && ( count1 + count2 <= 16 ) && ( ( ramFiles[fnum][fpos - 1] & ~0x80 ) != ( fileNum_ / 5 ) ) )
    {
        // we can combine the 2 items
        if ( LOCAL_DEBUG ) clog << " Insertion of a letter that can be added to the run-length-encoded char" << endl;
        ramFiles[fnum][fpos - 2] = ( ( count1 + count2 - 1 ) << 4 ) | letter1;
        return;
    }

    // Insertion of a letter between 2 letters
    if ( remainingRunLengthInReader_ == 0 )
    {
        if ( LOCAL_DEBUG ) clog << " Insertion of a letter between 2 letters" << endl;
        //  - where no insertion was already present
        if ( ramFiles[fnum][fpos - 1] == 0 )
        {
            if ( LOCAL_DEBUG ) clog << "  - where no insertion was already present" << endl;
            unsigned char replacementMetadata = fileNum_ / 5;
            assert( ( replacementMetadata & 0x80 ) == 0 && "Error: there may be too many cycles for this algorithm, which overwrote the Return bit" );
            ramFiles[fnum][fpos - 1] = replacementMetadata;

            ramFiles[fileNum_].push_back( c );
            ramFiles[fileNum_].push_back( metadata );
            return;
        }

        //  - where an insertion was already present, but no return
        if ( ramFiles[fnum][fpos - 1] != 0 && ( ramFiles[fnum][fpos - 1] & 0x80 ) == 0 )
        {
            if ( LOCAL_DEBUG ) clog << "  - where an insertion was already present, but no return" << endl;
            unsigned char replacementMetadata = fileNum_ / 5;
            assert( ( replacementMetadata & 0x80 ) == 0 && "Error: there may be too many cycles for this algorithm, which overwrote the Return bit" );
            if ( ramFiles[fnum][fpos - 1] != replacementMetadata )
            {
                onHoldUntilNextReturn_metadata_ = ramFiles[fnum][fpos - 1];
                ramFiles[fnum][fpos - 1] = replacementMetadata;
            }
            else
            {
                //     - where the insertion point refers to the current newest file - e.g. in the case of an insertion in the middle of a run-length-encoded char followed by readAndSend the rest of the char
                if ( LOCAL_DEBUG ) clog << "    - where the insertion point refers to the current newest file" << endl;
                // in this case we continue the previous insertion
                assert( ramFiles[fileNum_].back() & 0x80 );
                onHoldUntilNextReturn_metadata_ = ramFiles[fileNum_].back() & ~0x80;
                ramFiles[fileNum_].back() = 0;
            }

            ramFiles[fileNum_].push_back( c );
            ramFiles[fileNum_].push_back( metadata );
            return;
        }

        //  - where a return was already present
        if ( ramFiles[fnum][fpos - 1] == 0x80 )
        {
            if ( LOCAL_DEBUG ) clog << "  - where a return was already present" << endl;
            unsigned char replacementMetadata = fileNum_ / 5;
            assert( ( replacementMetadata & 0x80 ) == 0 && "Error: there may be too many cycles for this algorithm, which overwrote the Return bit" );
            ramFiles[fnum][fpos - 1] = replacementMetadata | 0x80; // the return stays

            ramFiles[fileNum_].push_back( c );
            ramFiles[fileNum_].push_back( metadata );
            return;
        }

        //  - where an insertion+return was already present
        if ( ( ramFiles[fnum][fpos - 1] & ~0x80 ) != 0 && ( ramFiles[fnum][fpos - 1] & 0x80 ) != 0 )
        {
            if ( LOCAL_DEBUG ) clog << "  - where an insertion+return was already present" << endl;
            unsigned char replacementMetadata = fileNum_ / 5;
            assert( ( replacementMetadata & 0x80 ) == 0 && "Error: there may be too many cycles for this algorithm, which overwrote the Return bit" );
            if ( ( ramFiles[fnum][fpos - 1] & ~0x80 ) != replacementMetadata )
            {
                onHoldUntilNextReturn_metadata_ = ramFiles[fnum][fpos - 1] & ~0x80;
                ramFiles[fnum][fpos - 1] = replacementMetadata | 0x80; // the return stays
            }
            else
            {
                //     - where the insertion point refers to the current newest file - e.g. in the case of an insertion in the middle of a run-length-encoded char followed by readAndSend the rest of the char
                if ( LOCAL_DEBUG ) clog << "    - where the insertion point refers to the current newest file" << endl;
                // in this case we continue the previous insertion
                assert( ramFiles[fileNum_].back() & 0x80 );
                onHoldUntilNextReturn_metadata_ = ramFiles[fileNum_].back() & ~0x80;
                ramFiles[fileNum_].back() = 0;
            }

            ramFiles[fileNum_].push_back( c );
            ramFiles[fileNum_].push_back( metadata );
            return;
        }

        assert( false && "Should never reach here" );
    }

    // Insertion of a letter in the middle of a run-length-encoded char
    if ( remainingRunLengthInReader_ != 0 )
    {
        if ( LOCAL_DEBUG ) clog << " Insertion of a letter in the middle of a run-length-encoded char" << endl;
        //  - where no insertion was already present
        if ( ramFiles[fnum][fpos - 1] == 0 )
        {
            if ( LOCAL_DEBUG ) clog << "  - where no insertion was already present" << endl;
            unsigned char replacementMetadata = fileNum_ / 5;
            assert( ( replacementMetadata & 0x80 ) == 0 && "Error: there may be too many cycles for this algorithm, which overwrote the Return bit" );

            onHoldUntilNextReturn_letter_ = letter1;
            onHoldUntilNextReturn_runLength_ = remainingRunLengthInReader_;
            assert( remainingRunLengthInReader_ < count1 );
            ramFiles[fnum][fpos - 2] = ( ( count1 - remainingRunLengthInReader_ - 1 ) << 4 ) | letter1;
            ramFiles[fnum][fpos - 1] = replacementMetadata;

            ramFiles[fileNum_].push_back( c );
            ramFiles[fileNum_].push_back( metadata );
            return;
        }

        //  - where an insertion was already present, but no return
        if ( ramFiles[fnum][fpos - 1] != 0 && ( ramFiles[fnum][fpos - 1] & 0x80 ) == 0 )
        {
            if ( LOCAL_DEBUG ) clog << "  - where an insertion was already present, but no return" << endl;
            unsigned char replacementMetadata = fileNum_ / 5;
            assert( ( replacementMetadata & 0x80 ) == 0 && "Error: there may be too many cycles for this algorithm, which overwrote the Return bit" );

            //     - if the insertion point is the current newest file
            if ( ramFiles[fnum][fpos - 1] == replacementMetadata )
            {
                if ( LOCAL_DEBUG ) clog << "     - if the insertion point is the current newest file" << endl;
                assert( ramFiles[fileNum_].back() != 0 ); // we should already have properly inserted the return, with an eventual jump, and we'll recycle those
                onHoldUntilNextReturn_metadata_ = ramFiles[fileNum_].back(); // todo: & ~0x80?
                ramFiles[fileNum_].pop_back();
                assert( ( ramFiles[fileNum_].back() & 0x0F ) == letter1 );
                assert( ( ramFiles[fileNum_].back() >> 4 ) + 1 > remainingRunLengthInReader_ );
                unsigned char newCount = ( ramFiles[fileNum_].back() >> 4 ) + 1 - remainingRunLengthInReader_;
                onHoldUntilNextReturn_letter_ = letter1;
                onHoldUntilNextReturn_runLength_ = remainingRunLengthInReader_;
                ramFiles[fileNum_].pop_back();

                ramFiles[fileNum_].push_back( ( newCount - 1 ) << 4 | letter1 );
                ramFiles[fileNum_].push_back( 0 );

                ramFiles[fileNum_].push_back( c );
                ramFiles[fileNum_].push_back( metadata );
                return;
            }

            //     - if the insertion point is an older file
            else
            {
                if ( LOCAL_DEBUG ) clog << "     - if the insertion point is an older file" << endl;
                onHoldUntilNextReturn_letter_ = letter1;
                onHoldUntilNextReturn_runLength_ = remainingRunLengthInReader_;
                onHoldUntilNextReturn_metadata_ = ramFiles[fnum][fpos - 1];
                assert( remainingRunLengthInReader_ < count1 );
                ramFiles[fnum][fpos - 2] = ( ( count1 - remainingRunLengthInReader_ - 1 ) << 4 ) | letter1;
                ramFiles[fnum][fpos - 1] = replacementMetadata;

                ramFiles[fileNum_].push_back( c );
                ramFiles[fileNum_].push_back( metadata );
                return;
            }
        }

        //  - where a return was already present
        if ( ramFiles[fnum][fpos - 1] == 0x80 )
        {
            if ( LOCAL_DEBUG ) clog << "  - where a return was already present" << endl;
            unsigned char replacementMetadata = fileNum_ / 5;
            assert( ( replacementMetadata & 0x80 ) == 0 && "Error: there may be too many cycles for this algorithm, which overwrote the Return bit" );

            onHoldUntilNextReturn_letter_ = letter1;
            onHoldUntilNextReturn_runLength_ = remainingRunLengthInReader_;
            assert( remainingRunLengthInReader_ < count1 );
            ramFiles[fnum][fpos - 2] = ( ( count1 - remainingRunLengthInReader_ - 1 ) << 4 ) | letter1;
            ramFiles[fnum][fpos - 1] = replacementMetadata | 0x80; // the return stays

            ramFiles[fileNum_].push_back( c );
            ramFiles[fileNum_].push_back( metadata );
            return;
        }

        //  - where an insertion+return was already present
        if ( ( ramFiles[fnum][fpos - 1] & ~0x80 ) != 0 && ( ramFiles[fnum][fpos - 1] & 0x80 ) != 0 )
        {
            if ( LOCAL_DEBUG ) clog << "  - where an insertion+return was already present" << endl;
            unsigned char replacementMetadata = fileNum_ / 5;
            assert( ( replacementMetadata & 0x80 ) == 0 && "Error: there may be too many cycles for this algorithm, which overwrote the Return bit" );

            //     - if the insertion point is the current newest file
            if ( ( ramFiles[fnum][fpos - 1] & ~0x80 ) == replacementMetadata )
            {
                if ( LOCAL_DEBUG ) clog << "     - if the insertion point is the current newest file" << endl;
                assert( ramFiles[fileNum_].back() != 0 ); // we should already have properly inserted the return, with an eventual jump, and we'll recycle those
                onHoldUntilNextReturn_metadata_ = ramFiles[fileNum_].back(); // todo: & ~0x80?
                ramFiles[fileNum_].pop_back();
                assert( ( ramFiles[fileNum_].back() & 0x0F ) == letter1 );
                assert( ( ramFiles[fileNum_].back() >> 4 ) + 1 > remainingRunLengthInReader_ );
                unsigned char newCount = ( ramFiles[fileNum_].back() >> 4 ) + 1 - remainingRunLengthInReader_;
                onHoldUntilNextReturn_letter_ = letter1;
                onHoldUntilNextReturn_runLength_ = remainingRunLengthInReader_;
                ramFiles[fileNum_].pop_back();

                ramFiles[fileNum_].push_back( ( newCount - 1 ) << 4 | letter1 );
                ramFiles[fileNum_].push_back( 0 );

                ramFiles[fileNum_].push_back( c );
                ramFiles[fileNum_].push_back( metadata );
                return;
            }

            //     - if the insertion point is an older file
            else
            {
                if ( LOCAL_DEBUG ) clog << "     - if the insertion point is an older file" << endl;
                onHoldUntilNextReturn_letter_ = letter1;
                onHoldUntilNextReturn_runLength_ = remainingRunLengthInReader_;
                onHoldUntilNextReturn_metadata_ = ramFiles[fnum][fpos - 1] & ~0x80;
                assert( remainingRunLengthInReader_ < count1 );
                ramFiles[fnum][fpos - 2] = ( ( count1 - remainingRunLengthInReader_ - 1 ) << 4 ) | letter1;
                ramFiles[fnum][fpos - 1] = replacementMetadata | 0x80; // the return stays

                ramFiles[fileNum_].push_back( c );
                ramFiles[fileNum_].push_back( metadata );
                return;
            }
        }

        assert( false && "Should never reach here" );
    }
}

void BwtWriterIncrementalRunLength::encodeRun( char c, uint runLength )
{
#ifdef DEBUG
    std::cout << "BW RL encodeRun - sending run " << c << " " << runLength << " " << pFile_
              << std::endl;
#endif
#ifdef REPORT_COMPRESSION_RATIO
    charsReceived_ += runLength;
#endif
    int charIndex( whichPile[( int )c] );

    if ( charIndex == nv )
    {
        cerr << "Char is not part of the alphabet. Aborting." << endl;
        exit( -1 );
    }

    uchar outCode( 0xF0 | ( ( uchar )charIndex ) );
    runLength--;
    const uint numMaxChars( runLength >> 4 );
    for ( uint i( 0 ); i < numMaxChars; i++ )
    {
        sendChar( outCode, 0 );
    }
    runLength &= ( uint )0xF;

    outCode = ( ( ( uchar )runLength ) << 4 );
    outCode |= charIndex;
    //    assert(((uint)outCode)<256);
    //    assert(fwrite( &outCode, sizeof(char), 1, pFile_ )==1);
    sendChar( outCode, 0 );
#ifdef DEBUG
    std::cout << "B sending " << ( uint )outCode << " " << pFile_ << std::endl;
#endif
} // ~encodeRun

void BwtWriterIncrementalRunLength::sendRun( char c, int runLength )
{
    assert( false );
#ifdef DEBUG
    std::cout << "BW RL sendRun - sending run " << c << " " << runLength << " " << endl;
#endif
    if ( runLength != 0 )
    {
        if ( c == lastChar_ )
        {
            runLength_ += runLength;
        }
        else
        {
            if ( runLength_ != 0 ) encodeRun( lastChar_, runLength_ );
            lastChar_ = c;
            runLength_ = runLength;
        }
    }
} // ~sendRun


// sendRunOfPreExistingData is expected to be called by reader's readAndSend
void BwtWriterIncrementalRunLength::sendRunOfPreExistingData( char c, int runLength, int fileNum, size_t posInRamFile, int remainingRunLength )
{
    if ( runLength == 0 )
    {
        return;
    }

    if ( runLength_ != 0 && ( fileNumInReader_ != ( uint )fileNum || filePosInReader_ != posInRamFile || remainingRunLengthInReader_ != remainingRunLength ) )
    {
        encodeRun( lastChar_, runLength_ );
        lastChar_  = notInAlphabet;
        runLength_ = 0;
    }

    // We shouldn't need to write anything
    // For debugging, we just want to check that the data is already present and advance our cursor
    if ( fileNum >= 0 )
    {
        fileNumInReader_ = fileNum;
        filePosInReader_ = posInRamFile;
        remainingRunLengthInReader_ = remainingRunLength;
    }

    lastFileReturnNeeded_ = true;
} // ~sendRunOfPreExistingData




// Huffman implementation

BwtWriterHuffman::~BwtWriterHuffman() // destructor
{
    sendRun( lastChar_, runLength_ ); //gets last normal chars from buffer
    sendRun( notInAlphabet, 1 ); // send termination char
    BwtWriterHuffman::emptyBuffer();
    while ( bitsUsed_ > 0 )
    {
        assert( fwrite( &soFar_.ui, sizeof( unsigned int ), 1, pFile_ ) == 1 );
#ifdef REPORT_COMPRESSION_RATIO
        bytesWritten_ += ( LetterCountType ) sizeof( unsigned int );
#endif
#ifdef DEBUG
        cout << endl;
        for ( unsigned int i( 1 ); i != 0; i <<= 1 )
            cout << ( ( soFar_.ui & i ) ? '1' : '0' );
        cout << endl;
#endif
        bitsUsed_ -= 32;
    }
#ifdef REPORT_COMPRESSION_RATIO
    Logger::out( LOG_FOR_DEBUGGING )
            << "BwtWriterHuffman: received "
            << charsReceived_ << " chars, sent "
            << bytesWritten_ << " bytes, compression "
            << ( ( double )8 * bytesWritten_ ) / ( charsReceived_ )
            << " bits per char " << std::endl;
#endif
} // ~BwtWriterHuffman()

void BwtWriterHuffman::operator()( const char *p, int numChars )
{
    for ( int i( 0 ); i < numChars; i++ )
    {
        if ( ( *p ) == lastChar_ )
        {
            runLength_++;
        }
        else
        {
            if ( runLength_ > 0 )
            {
                sendRun( lastChar_, runLength_ );
            } // ~if

            runLength_ = 1;
            lastChar_ = *p;
        } // ~else
        p++;
    } // ~for
    sendRun( lastChar_, runLength_ );
    runLength_ = 0;
} // ~operator()

void BwtWriterHuffman::sendToken( unsigned long long code, unsigned int length )
{
    toAdd_.ull = code;

    toAdd_.ull <<= bitsUsed_; // left shift to the next free position
    soFar_.ull |= toAdd_.ull; // update so far
    bitsUsed_ += length;

    if ( bitsUsed_ > 32 ) // if we have more than 32bit / 4 byte
    {
        assert( fwrite( &soFar_.ui, sizeof ( unsigned int ), 1, pFile_ ) == 1 );
#ifdef REPORT_COMPRESSION_RATIO
        bytesWritten_ += ( LetterCountType ) sizeof( unsigned int );
#endif
#ifdef DEBUG
        for ( unsigned int i( 1 ); i != 0; i <<= 1 )
            cout << ( ( soFar_.ui & i ) ? '1' : '0' );
        cout << endl;
#endif
        soFar_.ull >>= 32; // shift rest to the right border
        bitsUsed_ -= 32; // update bits used
    } //

} // ~sendToken()

void BwtWriterHuffman::sendRun( char c, int runLength )
{
#ifdef REPORT_COMPRESSION_RATIO
    charsReceived_ += runLength;
#endif
    if ( runLength > 0 )
    {

        for ( int i( 0 ); i < runLength; i++ )
        {

            if ( huffmanBufferPos == huffmanWriterBufferSize - 1 )
            {
                symBuf[huffmanBufferPos] = c;
                processBuffer( huffmanWriterBufferSize );
            }
            else
            {
                symBuf[huffmanBufferPos] = c;
                huffmanBufferPos++;
            } // ~else

        } // ~for
    }
} // ~sendRun

void BwtWriterHuffman::emptyBuffer( void )
{
    processBuffer( huffmanBufferPos );
}

void BwtWriterHuffman::processBuffer( int itemsToPrint )
{
    //  cerr << pFile_ << "PROCESSING BUFFER" << endl;

    if ( itemsToPrint > 0 )
    {
        char localLastChar = 0;
        unsigned int localRunLength = 0;

        for ( int i( 0 ); i < itemsToPrint; i++ )
        {
            //  cerr << pFile_ << " accessing char at " << i << endl;
            if ( symBuf[i] == localLastChar )
            {
                localRunLength++;
            }
            else
            {
                if ( localRunLength > 0 )
                {
                    // get number of this char, 0-5
                    int charIndex( whichPile[( int ) localLastChar] );
                    assert( charIndex != nv ); // crash if not from alphabet

                    if ( localRunLength == 1 ) // single run only
                    {
                        sendToken( singleCharCode[charIndex],
                                   singleCharLength[charIndex] );
                    }// ~if
                    else
                    {
                        sendToken( doubleCharCode[charIndex],
                                   doubleCharLength[charIndex] );
                        sendNum( localRunLength );
                    } //~else
                } // ~if

                localRunLength = 1;
                localLastChar = symBuf[i];
            } // ~else
        } // ~for

        // process last entry of the buffer
        int charIndex( whichPile[( int ) localLastChar] ); // get number of this char, 0-5

        if ( ( int )localLastChar > 0 && whichPile[( int ) localLastChar] < alphabetSize )
        {

            assert( charIndex != nv ); // crash if not from alphabet

            if ( localRunLength == 1 ) // single run only
            {
                sendToken( singleCharCode[charIndex], singleCharLength[charIndex] );
            }// ~if
            else
            {
                sendToken( doubleCharCode[charIndex], doubleCharLength[charIndex] );
                sendNum( localRunLength );
            } //~else
        }
        huffmanBufferPos = 0; // reset counter
    }
} // ~sendRun

void BwtWriterHuffman::sendNum( unsigned int runLength )
{
    if ( runLength < 17 ) // max 16
    {
        runLength--;
        runLength--;
        numBuf_.ui = runLength; // set  new run length
        sendToken( numBuf_.ull, 4 ); // write one token, encoding for the runlength
    }// ~if
    else // larger than 16 -> 2 byte
    {
        runLength -= 17; // substract 16 + 1
        numBuf_.ui = 0xF; // set unsigned int to 16 -> 1111
        sendToken( numBuf_.ull, 4 ); // send binary encoded 16 using 4 bits
        do
        {
            numBuf_.ui = runLength; // set unsigned int to remaining runlength
            numBuf_.ui &= 0x7F; // AND with 111|1111
            // 1st: OR with 1000|0000 2nd multiply with 1 128 if RL > 127 or with 0
            numBuf_.ui |= 0x80 * ( runLength > 0x7F );
            sendToken( numBuf_.ull, 8 );
            runLength >>= 7;
        }// ~while
        while ( runLength != 0 );
    } // ~else
} // sendNum

//
// BwtWriterImplicit member function definitions
//

BwtWriterImplicit::~BwtWriterImplicit()
{
    if ( inSAP_ == true )
    {
        flushSAP();
    }
    else if ( lastChar_ != notInAlphabet )
    {
        pWriter_->sendRun( lastChar_, lastRun_ );
    }
    delete pWriter_;
}

void BwtWriterImplicit::flushSAP( void )
{
    assert( alphabet[firstSAP_] == lastChar_ );
    if ( countSAP_.count_[firstSAP_] > 0 ) pWriter_->sendRun( alphabet[firstSAP_], countSAP_.count_[firstSAP_] );

    for ( int i( 0 ); i < alphabetSize; i++ )
    {
        if ( ( i != firstSAP_ ) && ( countSAP_.count_[i] > 0 ) ) pWriter_->sendRun( alphabet[i], countSAP_.count_[i] );
    }
}


void BwtWriterImplicit::operator()( const char *p, int numChars )
{

    for ( int i( 0 ); i < numChars; i++, p++ )
    {
        if ( islower( *p ) )
        {
            if ( inSAP_ == false )
            {
                countSAP_.clear();
                assert ( lastChar_ != notInAlphabet );
                firstSAP_ = whichPile[( int )lastChar_];
                assert( firstSAP_ != nv );
                countSAP_.count_[firstSAP_] += lastRun_;
                inSAP_ = true;
            } // ~if
            countSAP_ += *p;
        } // ~if
        else
        {
            if ( inSAP_ == true )
            {
                flushSAP();
                inSAP_ = false;
            }
            else if ( lastChar_ != notInAlphabet )
            {
                pWriter_->sendRun( lastChar_, lastRun_ );
            }
            lastChar_ = *p;
            lastRun_ = 1;
        }
    }
}

void BwtWriterImplicit::sendRun( char c, int runLength )
{
    if ( islower( c ) )
    {
        if ( inSAP_ == false )
        {
            countSAP_.clear();
            assert ( lastChar_ != notInAlphabet );
            firstSAP_ = whichPile[( int )lastChar_];
            assert( firstSAP_ != nv );
            countSAP_.count_[firstSAP_] += lastRun_;
            inSAP_ = true;
        } // ~if
        countSAP_.count_[whichPile[( int )c]] += runLength;
    }
    else
    {
        if ( inSAP_ == true )
        {
            flushSAP();
            inSAP_ = false;
        }
        else if ( lastChar_ != notInAlphabet )
        {
            pWriter_->sendRun( lastChar_, lastRun_ );
        }
        lastChar_ = c;
        lastRun_ = runLength;
    }

    //  (*pWriter_).sendRun(toupper(c), runLength);
}

#ifdef XXX
void BwtWriterImplicit::operator()( const char *p, int numChars )
{
    // could be smarter about this
    char c;
    for ( int i( 0 ); i < numChars; i++, p++ )
    {
        c = toupper( *p );
        ( *pWriter_ )( &c, 1 );
    }
}

void BwtWriterImplicit::sendRun( char c, int runLength )
{
    ( *pWriter_ ).sendRun( toupper( c ), runLength );
}
#endif
