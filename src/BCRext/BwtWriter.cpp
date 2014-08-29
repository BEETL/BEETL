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

#include "BwtWriter.hh"

#include "LetterCount.hh"
#include "Tools.hh"
#include "libzoo/util/Logger.hh"

#include <cstdlib>
#include <cstring>

using namespace std;


//#define DEBUG 1
#define LOCAL_DEBUG 0

//
// BwtWriterBase member function definitions
//

BwtWriterFile::BwtWriterFile( const string &fileName ) : pFile_( fopen( fileName.c_str(), "wb" ) )
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

void BwtWriterFile::flush()
{

    fflush( pFile_ );
}

//
// BwtWriterASCII member function definitions
//

BwtWriterASCII::BwtWriterASCII( const string &fileName ) : BwtWriterFile( fileName ), lastChar_( notInAlphabet )
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


void BwtWriterASCII::operator()( const char *p, LetterNumber numChars )
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

    if ( numChars > 0 )
        lastChar_ = p[numChars - 1];
} // ~operator()

void BwtWriterASCII::sendRun( char c, LetterNumber runLength )
{
    if ( runLength )
    {
        for ( LetterNumber i( 0 ); i < runLength; i++ ) fputc( c, pFile_ );
        lastChar_ = c;
    }
}

char BwtWriterASCII::getLastChar()
{
    return lastChar_;
}


//
// BwtWriterRunLengthBase member function definitions
//

BwtWriterRunLengthBase::~BwtWriterRunLengthBase()
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
        bytesWritten_ += ( LetterNumber )( pBuf_ - buf_ );
#endif
    }

#ifdef REPORT_COMPRESSION_RATIO
    Logger_if( LOG_FOR_DEBUGGING ) Logger::out()
            << "BwtWriterRunLengthBase: received "
            << charsReceived_ << " chars, sent "
            << bytesWritten_ << " bytes, compression "
            << ( ( double )8 * bytesWritten_ ) / ( charsReceived_ )
            << " bits per char " << std::endl;
#endif

#ifdef GENERATE_RLE_HISTOGRAM
    for ( auto & kv : histogram_ )
    {
        cout << "histogram:\t" << kv.first.first << "\t" << kv.first.second << "\t" << kv.second << endl;
    }
#endif
}

void BwtWriterRunLengthBase::operator()( const char *p, LetterNumber numChars )
{
#ifdef DEBUG
    std::cout << "BW RL () - " << *p << " " << numChars << " state: " << lastChar_ << " " << runLength_ << endl;
#endif

    for ( LetterNumber i( 0 ); i < numChars; i++ )
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

void BwtWriterRunLengthBase::flushBuffer()
{
    size_t bufLength = pBuf_ - buf_;
#ifdef REPORT_COMPRESSION_RATIO
    bytesWritten_ += bufLength;
#endif

    size_t bytesWritten = fwrite( buf_, sizeof( char ), bufLength, pFile_ );
    if ( bytesWritten != ( size_t )bufLength )
    {
        cerr << "Unable to write " << bufLength
             << " chars. Aborting." << endl;
        exit( EXIT_FAILURE );
    }
    pBuf_ = buf_;
}

void BwtWriterRunLengthBase::sendChar( char c )
{
    *pBuf_ = c;
    if ( ++pBuf_ == pBufMax_ )
        flushBuffer();
}

void BwtWriterRunLengthBase::encodeRun( char c, LetterNumber runLength )
{
#ifdef DEBUG
    std::cout << "BW RL encodeRun - sending run " << c << " " << runLength << " " << pFile_
              << std::endl;
#endif
#ifdef GENERATE_RLE_HISTOGRAM
    ++histogram_[ make_pair( c, runLength ) ];
#endif
#ifdef REPORT_COMPRESSION_RATIO
    charsReceived_ += runLength;
#endif
    const unsigned char charIndex( whichPile[( int )c] );

    if ( charIndex == nv )
    {
        cerr << "Char |" << c << "| is not part of the alphabet. Aborting." << endl;
        exit( EXIT_FAILURE );
    }

    assert( charIndex >> baseFieldWidthInBits_ == 0 );
    uchar outCode = charIndex | lengthFieldMask_;
    runLength--;
    const LetterNumber numMaxChars( runLength >> lengthFieldWidthInBits_ );
    for ( LetterNumber i( 0 ); i < numMaxChars; i++ )
    {
        sendChar( outCode );
    }
    runLength -= numMaxChars << lengthFieldWidthInBits_;
    assert( runLength >> lengthFieldWidthInBits_ == 0 );
    outCode = charIndex | ( static_cast<uchar>( runLength ) << baseFieldWidthInBits_ );

    sendChar( outCode );
#ifdef DEBUG
    std::cout << "B sending " << ( unsigned int )outCode << " " << pFile_ << std::endl;
#endif
} // ~encodeRun

void BwtWriterRunLengthBase::sendRun( char c, LetterNumber runLength )
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


char BwtWriterRunLengthBase::getLastChar()
{
    return lastChar_;
}

void BwtWriterRunLengthBase::flush()
{
    if ( runLength_ != 0 )
    {
        encodeRun( lastChar_, runLength_ );
        lastChar_ = notInAlphabet;
        runLength_ = 0;
    }
    flushBuffer();

    fflush( pFile_ );
}


//
// BwtWriterRunLengthV2 member function definitions
//
BwtWriterRunLengthV2::BwtWriterRunLengthV2( const string &fileName )
    : BwtWriterRunLengthBase( fileName, 3 )
{
    symbolForRunLength1ForPile_.resize( alphabetSize );
    maxEncodedRunLengthForPile_.resize( alphabetSize );

    LetterNumber bytecode = 0;
    symbolForRunLength1ForPile_[ whichPile['A'] ] = bytecode;
    bytecode += ( maxEncodedRunLengthForPile_[ whichPile['A'] ] = 63 );

    symbolForRunLength1ForPile_[ whichPile['C'] ] = bytecode;
    bytecode += ( maxEncodedRunLengthForPile_[ whichPile['C'] ] = 63 );

    symbolForRunLength1ForPile_[ whichPile['G'] ] = bytecode;
    bytecode += ( maxEncodedRunLengthForPile_[ whichPile['G'] ] = 63 );

    symbolForRunLength1ForPile_[ whichPile['T'] ] = bytecode;
    bytecode += ( maxEncodedRunLengthForPile_[ whichPile['T'] ] = 63 );

    symbolForRunLength1ForPile_[ whichPile['N'] ] = bytecode;
    bytecode += ( maxEncodedRunLengthForPile_[ whichPile['N'] ] = 1 );

    symbolForRunLength1ForPile_[ whichPile['$'] ] = bytecode;
    bytecode += ( maxEncodedRunLengthForPile_[ whichPile['$'] ] = 3 );

    assert( bytecode == 256 );

    // Send file header (inspired by the PNG format)
    fputc( 137, pFile_ ); // non-ASCII to avoid confusion with text files
    fputc( 'B', pFile_ );
    fputc( 'W', pFile_ );
    fputc( 'T', pFile_ );
    fputc( 13, pFile_ ); // \r\n sequence to check for invalid dos/unix format conversions
    fputc( 10, pFile_ );
    fputc( 26, pFile_ ); // Ctrl-Z, making some text viewers stop here

    fputc( 2, pFile_ ); // Format version number

    /*
        // Send conversion table: 256 entries with { base: 1 char, run length: 1 byte }
        for (int i=1; i<=63; ++i)
        {
            fputc( 'A', pFile_ );
            fputc( i, pFile_ );
        }
        for (int i=1; i<=63; ++i)
        {
            fputc( 'C', pFile_ );
            fputc( i, pFile_ );
        }
        for (int i=1; i<=63; ++i)
        {
            fputc( 'G', pFile_ );
            fputc( i, pFile_ );
        }
        for (int i=1; i<=63; ++i)
        {
            fputc( 'T', pFile_ );
            fputc( i, pFile_ );
        }
        for (int i=1; i<=1; ++i)
        {
            fputc( 'N', pFile_ );
            fputc( i, pFile_ );
        }
        for (int i=1; i<=3; ++i)
        {
            fputc( '$', pFile_ );
            fputc( i, pFile_ );
        }
    */
}

BwtWriterRunLengthV2::~BwtWriterRunLengthV2()
{
    // We need to override this call from the base class' destructor, in order to use our own encodeRun
    if ( runLength_ != 0 )
    {
        encodeRun( lastChar_, runLength_ );
        runLength_ = 0;
    }
}

void BwtWriterRunLengthV2::encodeRun( char c, LetterNumber runLength )
{
    assert( runLength > 0 );
#ifdef DEBUG
    std::cout << "BW RL encodeRun - sending run " << c << " " << runLength << " " << pFile_
              << std::endl;
#endif
#ifdef GENERATE_RLE_HISTOGRAM
    ++histogram_[ make_pair( c, runLength ) ];
#endif
#ifdef REPORT_COMPRESSION_RATIO
    charsReceived_ += runLength;
#endif
    const unsigned char charIndex( whichPile[( int )c] );

    if ( charIndex == nv )
    {
        cerr << "Char |" << c << "| is not part of the alphabet. Aborting." << endl;
        exit( EXIT_FAILURE );
    }

    const uchar symbolForRunLength1 = symbolForRunLength1ForPile_[charIndex];
    const LetterNumber maxEncodedRunLength = maxEncodedRunLengthForPile_[charIndex];

    uchar outCode = symbolForRunLength1 + maxEncodedRunLength - 1;
    LetterNumber countAtFullLength = runLength / maxEncodedRunLength;
    for ( LetterNumber i = 0; i < countAtFullLength; ++i )
    {
        sendChar( outCode );
    }

    runLength %= maxEncodedRunLength;
    if ( runLength > 0 )
    {
        outCode = symbolForRunLength1 + runLength - 1;
        sendChar( outCode );
    }
}


//
// BwtWriterRunLengthV3 member function definitions
//
BwtWriterRunLengthV3::BwtWriterRunLengthV3( const string &fileName )
    : BwtWriterRunLengthBase( fileName, 3 )
{
    symbolForRunLength1ForPile_.resize( alphabetSize );
    maxEncodedRunLengthForPile_.resize( alphabetSize );

    // Send file header (inspired by the PNG format)
    fputc( 'B', pFile_ );
    fputc( 'W', pFile_ );
    fputc( 'T', pFile_ );
    fputc( 13, pFile_ ); // \r\n sequence to check for invalid dos/unix format conversions
    fputc( 10, pFile_ );
    fputc( 26, pFile_ ); // Ctrl-Z, making some text viewers stop here and being non-ASCII to avoid confusion with text files

    // Format version number, on 2 bytes to help identify endianness problems
    uint16_t formatVersion = 3;
    assert( fwrite( &formatVersion, sizeof( formatVersion ), 1, pFile_ ) == 1 );

    // Initialise conversion table: enough ranges to cover 256 entries, following the format { base: 1 char, range length: 1 byte, first run length: 2 bytes }
    uint16_t bytecode = 0;
    bytecode = initialiseCodeRange( 'A', 58, 1, bytecode );
    bytecode = initialiseCodeRange( 'C', 58, 1, bytecode );
    bytecode = initialiseCodeRange( 'G', 58, 1, bytecode );
    bytecode = initialiseCodeRange( 'T', 58, 1, bytecode );
    bytecode = initialiseCodeRange( 'N', 4, 1, bytecode );
    bytecode = initialiseCodeRange( '$', 4, 1, bytecode );
    bytecode = initialiseCodeRange( '+', 16, 0, bytecode );
    assert( bytecode == 256 );
}

uint16_t BwtWriterRunLengthV3::initialiseCodeRange( const uint8_t base, const uint8_t rangeLength, const uint16_t firstRunLength, const uint8_t firstBytecode )
// returns next available bytecode, as uint16 to reach value 256
{
    assert( fwrite( &base, sizeof( base ), 1, pFile_ ) == 1 );
    assert( fwrite( &rangeLength, sizeof( rangeLength ), 1, pFile_ ) == 1 );
    assert( fwrite( &firstRunLength, sizeof( firstRunLength ), 1, pFile_ ) == 1 );

    if (base != '+')
    { // Common case for "normal" bases
        symbolForRunLength1ForPile_[ whichPile[base] ] = firstBytecode;
        maxEncodedRunLengthForPile_[ whichPile[base] ] = rangeLength;
    }
    else
    { // Special case for continuation symbol
        firstContinuationSymbol_ = firstBytecode;
        maxEncodedRunLengthMultiplierForContinuationSymbol_ = rangeLength;
    }

    return (uint16_t)firstBytecode + rangeLength;
}

BwtWriterRunLengthV3::~BwtWriterRunLengthV3()
{
    // We need to override this call from the base class' destructor, in order to use our own encodeRun
    if ( runLength_ != 0 )
    {
        encodeRun( lastChar_, runLength_ );
        runLength_ = 0;
    }
}

void BwtWriterRunLengthV3::encodeRun( char c, LetterNumber runLength )
{
    assert( runLength > 0 );
    LetterNumber runLengthMinus1 = runLength - 1;
#ifdef DEBUG
    std::cout << "BW RL encodeRun - sending run " << c << " " << runLength << " " << pFile_
              << std::endl;
#endif
#ifdef GENERATE_RLE_HISTOGRAM
    ++histogram_[ make_pair( c, runLength ) ];
#endif
#ifdef REPORT_COMPRESSION_RATIO
    charsReceived_ += runLength;
#endif
    const unsigned char charIndex( whichPile[( int )c] );

    if ( charIndex == nv )
    {
        cerr << "Char |" << c << "| is not part of the alphabet. Aborting." << endl;
        exit( EXIT_FAILURE );
    }

    const uchar symbolForRunLength1 = symbolForRunLength1ForPile_[charIndex];
    const LetterNumber maxEncodedRunLength = maxEncodedRunLengthForPile_[charIndex];

    // First char, encoded according to table
    LetterNumber symbolOffset = runLengthMinus1 % maxEncodedRunLength;
    uchar outCode = symbolForRunLength1 + symbolOffset;
    sendChar( outCode );
    runLengthMinus1 /= maxEncodedRunLength;

    // Subsequent chars, encoded in base `maxEncodedRunLengthMultiplierForContinuationSymbol_`, little endian
    while (runLengthMinus1 > 0)
    {
        symbolOffset = runLengthMinus1 % maxEncodedRunLengthMultiplierForContinuationSymbol_;
        outCode = firstContinuationSymbol_ + symbolOffset;
        sendChar( outCode );
        runLengthMinus1 /= maxEncodedRunLengthMultiplierForContinuationSymbol_;
    }
}


//
// BwtWriterIncrementalRunLength member function definitions
//

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
    int firstDigitPos = fileName.size() - 1;
    while ( firstDigitPos >= 0 && ( fileName[firstDigitPos] >= '0' && fileName[firstDigitPos] <= '9' ) ) // skip any prefix like in "new_1"
        --firstDigitPos;
    ++firstDigitPos;
    unsigned int letterNum = atoi( fileName.c_str() + firstDigitPos );

    extern unsigned int debugCycle;
    extern unsigned int lastDefragCycle;
    assert( letterNum > 0 && "This class doesn't store pile 0" );
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
    fileNumInReader_ = fileNum_ % ( alphabetSize - 1 ); // todo: improve this modulo
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
        bytesWritten_ += ( LetterNumber )( pBuf_ - buf_ );
#endif
    }

#ifdef REPORT_COMPRESSION_RATIO

#ifndef SEND_DATA_TO_FILES_FOR_DEBUGGING
    bytesWritten_ = ramFiles[fileNum_].size();
#endif //ifndef SEND_DATA_TO_FILES_FOR_DEBUGGING

    Logger_if( LOG_FOR_DEBUGGING ) Logger::out()
            << "BwtWriterIncrementalRunLength: received "
            << charsReceived_ << " chars, sent "
            << bytesWritten_ << " bytes, compression "
            << ( ( double )8 * bytesWritten_ ) / ( charsReceived_ )
            << " bits per char " << std::endl;
#endif


}

void BwtWriterIncrementalRunLength::operator()( const char *p, LetterNumber numChars )
{
#ifdef DEBUG
    std::cout << "BW RL () - " << *p << " " << numChars << " state: " << lastChar_ << " " << runLength_ << endl;
#endif

    for ( LetterNumber i( 0 ); i < numChars; i++ )
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

    if ( LOCAL_DEBUG )
    {
        clog << "Inserting " << ( unsigned int )c << " in file " << fnum << " pos " << fpos << ", remainingRunLengthInReader=" << remainingRunLengthInReader_ << " fileNum=" << fileNum_ << " ramFiles[fileNum_].size()=" << ramFiles[fileNum_].size() << " ramFiles[fnum].size()=" << ramFiles[fnum].size() << endl;
    }

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
                assert( ( ramFiles[fileNum_].back() >> 4 ) + 1 > ( int )remainingRunLengthInReader_ );
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
                assert( ( ramFiles[fileNum_].back() >> 4 ) + 1 > ( int )remainingRunLengthInReader_ );
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

void BwtWriterIncrementalRunLength::encodeRun( char c, LetterNumber runLength )
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
    const LetterNumber numMaxChars( runLength >> 4 );
    for ( LetterNumber i( 0 ); i < numMaxChars; i++ )
    {
        sendChar( outCode, 0 );
    }
    runLength &= ( LetterNumber )0xF;

    outCode = ( ( ( uchar )runLength ) << 4 );
    outCode |= charIndex;
    //    assert(((uint)outCode)<256);
    //    assert(fwrite( &outCode, sizeof(char), 1, pFile_ )==1);
    sendChar( outCode, 0 );
#ifdef DEBUG
    std::cout << "B sending " << ( uint )outCode << " " << pFile_ << std::endl;
#endif
} // ~encodeRun

void BwtWriterIncrementalRunLength::sendRun( char c, LetterNumber runLength )
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
void BwtWriterIncrementalRunLength::sendRunOfPreExistingData( char c, LetterNumber runLength, int fileNum, size_t posInRamFile, LetterNumber remainingRunLength )
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

#ifdef ACTIVATE_HUFFMAN

BwtWriterHuffman::~BwtWriterHuffman() // destructor
{
    sendRun( lastChar_, runLength_ ); //gets last normal chars from buffer
    sendRun( notInAlphabet, 1 ); // send termination char
    BwtWriterHuffman::emptyBuffer();
    while ( bitsUsed_ > 0 )
    {
        assert( fwrite( &soFar_.ui, sizeof( unsigned int ), 1, pFile_ ) == 1 );
#ifdef REPORT_COMPRESSION_RATIO
        bytesWritten_ += ( LetterNumber ) sizeof( unsigned int );
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
    Logger_if( LOG_FOR_DEBUGGING ) Logger::out()
            << "BwtWriterHuffman: received "
            << charsReceived_ << " chars, sent "
            << bytesWritten_ << " bytes, compression "
            << ( ( double )8 * bytesWritten_ ) / ( charsReceived_ )
            << " bits per char " << std::endl;
#endif
} // ~BwtWriterHuffman()

void BwtWriterHuffman::operator()( const char *p, LetterNumber numChars )
{
    for ( LetterNumber i( 0 ); i < numChars; i++ )
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

void BwtWriterHuffman::sendToken( unsigned long long code, LetterNumber length )
{
    toAdd_.ull = code;

    toAdd_.ull <<= bitsUsed_; // left shift to the next free position
    soFar_.ull |= toAdd_.ull; // update so far
    bitsUsed_ += length;

    if ( bitsUsed_ > 32 ) // if we have more than 32bit / 4 byte
    {
        assert( fwrite( &soFar_.ui, sizeof ( unsigned int ), 1, pFile_ ) == 1 );
#ifdef REPORT_COMPRESSION_RATIO
        bytesWritten_ += ( LetterNumber ) sizeof( unsigned int );
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

void BwtWriterHuffman::sendRun( char c, LetterNumber runLength )
{
#ifdef REPORT_COMPRESSION_RATIO
    charsReceived_ += runLength;
#endif
    if ( runLength > 0 )
    {

        for ( LetterNumber i( 0 ); i < runLength; i++ )
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
        LetterNumber localRunLength = 0;

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

void BwtWriterHuffman::sendNum( LetterNumber runLength )
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

#endif //ifdef ACTIVATE_HUFFMAN


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


void BwtWriterImplicit::operator()( const char *p, LetterNumber numChars )
{

    for ( LetterNumber i( 0 ); i < numChars; i++, p++ )
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

void BwtWriterImplicit::sendRun( char c, LetterNumber runLength )
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
void BwtWriterImplicit::operator()( const char *p, LetterNumber numChars )
{
    // could be smarter about this
    char c;
    for ( LetterNumber i( 0 ); i < numChars; i++, p++ )
    {
        c = toupper( *p );
        ( *pWriter_ )( &c, 1 );
    }
}

void BwtWriterImplicit::sendRun( char c, LetterNumber runLength )
{
    ( *pWriter_ ).sendRun( toupper( c ), runLength );
}
#endif
