
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

#include "Alphabet.hh"
#include "BwtReader.hh"
#include "BwtWriter.hh"
#include "CountWords.hh"
#include "LetterCount.hh"

#include <cmath>
#include <cstdlib>
#include <cstring>

using namespace std;


//#define DEBUG 1

int main( int numArgs, const char *args[] )
{
    if ( numArgs<3 )
    {
        cerr << "Usage: " << args[0] << " indexBinSize <BWT files to index>" << endl;
        exit( EXIT_FAILURE );
    }
    string indexFileName;
    FILE *pFile;
    int indexBinSize( atoi( args[1] ) );
    assert( indexBinSize>0 );

    for ( int i( 2 ); i<numArgs; i++ )
    {
        BwtReaderRunLengthIndex reader( args[i] );
        indexFileName=( string )args[i];
        indexFileName+=".idx";
        cerr << "Will build index " << indexFileName << " from " << args[i]
             << endl;
        pFile=fopen( indexFileName.c_str(),"w" );
        if ( pFile==NULL )
        {
            cerr<< "Problem opening file " << indexFileName <<" for writing"
                << endl;
            exit( EXIT_FAILURE );
        }
        reader.buildIndex( pFile, indexBinSize );
        fclose ( pFile );
    }

}



#ifdef OLD

//
// BwtReaderBase member function definitions
//

BwtReaderBase::BwtReaderBase( const string &fileName ) :
    pFile_( fopen( fileName.c_str(),"r" ) )
{
    if ( pFile_==NULL )
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



unsigned int BwtReaderBase::readAndCount( LetterCount &c )
{
    // will call the relevant virtual function
    return readAndCount( c, maxLetterCountType );
} // ~readAndCount

unsigned int BwtReaderBase::readAndSend( BwtWriterBase &writer )
{
    // will call the relevant virtual function
    return readAndSend( writer,  2000000000 ); // TBD numChars->LetterCountType
} // ~readAndCount


//
// BwtReaderASCII member function definitions
//

void BwtReaderASCII::rewindFile( void )
{
    rewind( pFile_ );
    currentPos_=0;
} // ~rewindFile

LetterCountType BwtReaderASCII::tellg( void ) const
{
    return currentPos_;
} // ~tellg



unsigned int BwtReaderASCII::readAndCount( LetterCount &c, const LetterCountType numChars )
{
#ifdef DEBUG
    std::cout << "BR ASCII readAndCount " << numChars << " chars " << endl;
#endif
    LetterCountType charsLeft( numChars ), charsToRead, charsRead;
    while ( charsLeft>0 )
    {
        charsToRead=( ( charsLeft>ReadBufferSize )?ReadBufferSize:charsLeft );
        charsRead=fread( buf_, sizeof( char ), charsToRead, pFile_ );
#ifdef DEBUG
        std::cout << "Reading " << charsRead << " chars ";
#endif
        for ( LetterCountType i( 0 ); i<charsRead; i++ )
        {
#ifdef DEBUG
            std::cout << buf_[i];
#endif
            c.count_[whichPile[( int )buf_[i]]]++;
        }
#ifdef DEBUG
        std::cout << std::endl;
#endif
        charsLeft-=charsRead;
        if ( charsRead<charsToRead )
        {
            // did not get everything asked for! return num of chars actually found
            currentPos_+=( numChars-charsLeft );
            return ( numChars-charsLeft );
        } // ~if
    } // ~while
    currentPos_+=numChars;
    return numChars;
} // ~int BwtReaderASCII::readAndCount( LetterCount& c, const int numChars )

unsigned int BwtReaderASCII::readAndSend( BwtWriterBase &writer, const int numChars )
{
#ifdef DEBUG
    std::cout << "BR ASCII readAndSend " << numChars << " chars " << endl;
#endif

    unsigned int totalRead = 0;
    // read readbufferzise bytes, if there are less bytes
    // ordered only fetch the last missing bytes
    unsigned int readNextPass =
        ( ( numChars-totalRead )<ReadBufferSize )?( numChars-totalRead ):ReadBufferSize;

    // std::cout << "Reading " << numChars << " chars " << endl;


    while ( totalRead<( unsigned int )numChars )
    {
        unsigned int numRead = 0;

        // try to read buffersize byte from file
        numRead = fread( buf_, sizeof ( uchar ), readNextPass, pFile_ );
        totalRead+= numRead;
        if ( numRead == 0 ) break;

        readNextPass =
            ( ( numChars-totalRead )<ReadBufferSize )?( numChars-totalRead ):ReadBufferSize;
        //std::cout << "next pass " << numRead << " chars " << endl;

        //    writer( buf_, numRead );
#define XXX 1
#ifdef XXX
        unsigned int charsLeft = numRead;

        for ( unsigned int counter = 0; counter < charsLeft; counter++ )
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
    runLength_=0; // next call to his function will then again start a new run
#endif
    return totalRead;

} // ~int BwtReaderASCII::readAndSend( BwtWriterBase& writer, const int numChars )

int BwtReaderASCII::operator()( char *p, int numChars )
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
    pBuf_( buf_+ReadBufferSize ),pBufMax_( buf_+ReadBufferSize ),
    lastChar_( notInAlphabet ),
    finished_( false )
{
    unsigned int j;
    for ( unsigned int i( 0 ); i<256; i++ )
    {
        lengths_[i]=1+( i>>4 );
        j=( i&0xF );
        codes_[i]=( j<alphabetSize )?alphabet[j]:notInAlphabet;
    } // ~for i
} // ~ctor

void BwtReaderRunLength::rewindFile( void )
{
    // rewind file and set all vars as per constructor
    rewind( pFile_ );
    runLength_=0;
    pBuf_=buf_+ReadBufferSize;
    pBufMax_=buf_+ReadBufferSize;
    lastChar_=notInAlphabet;
    currentPos_=0;
    finished_=false;
} // ~rewindFile

LetterCountType BwtReaderRunLength::tellg( void ) const
{
    return currentPos_;
} // ~tellg

// TODO: return type should be unsigned in order to avoid compiler warnings
unsigned int BwtReaderRunLength::readAndCount( LetterCount &c, const LetterCountType numChars )
{
#ifdef DEBUG
    std::cout << "BR RL readAndCount " << numChars << " chars " << endl;
#endif
    LetterCountType charsLeft( numChars );
    while ( charsLeft>runLength_ )
    {
        // Below is not great design, at first call of this function it accesses an
        // out-of-range array element. Fortunately it always adds zero to it! :)
        c.count_[whichPile[lastChar_]]+=runLength_;
        charsLeft-=runLength_;
        if ( getRun()==false )
        {
            currentPos_+=( numChars-charsLeft );
            return ( numChars-charsLeft );
            // assert(1==0);
        } // ~if
    } // ~while

    c.count_[whichPile[lastChar_]]+=charsLeft;
    runLength_-=charsLeft;
    currentPos_+=numChars;
    return numChars;
} // ~BwtReaderRunLength::readAndCount( LetterCount& c, const int numChars )

unsigned int BwtReaderRunLength::readAndSend( BwtWriterBase &writer, const int numChars )
{
#ifdef DEBUG
    std::cout << "BR RL readAndSend " << numChars << " chars " << endl;
#endif
    unsigned int charsLeft( numChars );
    while ( charsLeft>runLength_ )
    {
        //      int fred(whichPile[lastChar_]);
        writer.sendRun( lastChar_, runLength_ );
        //      c.count_[whichPile[lastChar_]]+=runLength_;
        charsLeft-=runLength_;
        if ( getRun()==false )
        {
            currentPos_+=( numChars-charsLeft );
            return ( numChars-charsLeft );
            // assert(1==0);
        } // ~if
    } // ~while

    writer.sendRun( lastChar_, charsLeft );
    //    c.count_[whichPile[lastChar_]]+=charsLeft;
    runLength_-=charsLeft;
    currentPos_+=numChars;
    return numChars;
} //~BwtReaderRunLength::readAndSend(BwtWriterBase& writer, const int numChars)

int BwtReaderRunLength::operator()( char *p, int numChars )
{
#ifdef DEBUG
    std::cout << "BR RL () :  asked for " << numChars << " " << lastChar_ << " "
              << runLength_ << " " << pFile_ << std::endl;
#endif
    unsigned int charsLeft( numChars );
    //    return fread( p, sizeof(char), numChars, pFile_ );
    while ( charsLeft>runLength_ )
    {
#ifdef DEBUG
        std::cout << "BR RL () :  setting " << lastChar_ << " "
                  << runLength_ << " " << pFile_ << std::endl;
#endif

        memset( p, lastChar_, runLength_ );
        p+=runLength_;

        charsLeft-=runLength_;
        if ( getRun()==false )
        {
            // runLength_=0;
#ifdef DEBUG
            std::cout << "B read " << numChars-charsLeft << " out of "
                      << numChars << std::endl;
#endif
            currentPos_+=( numChars-charsLeft );
            return ( numChars-charsLeft );
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

    runLength_-=charsLeft;
#ifdef DEBUG
    std::cout << "B delivered " << numChars << " " << charsLeft << " "
              << pFile_ << std::endl;
#endif
    currentPos_+=numChars;
    return numChars;
} // ~operator()

bool BwtReaderRunLength::getRun( void )
{
    if ( pBuf_==pBufMax_ )
    {
        if ( finished_ ) return false;
        else
        {
            unsigned int numRead( fread( buf_, sizeof( uchar ),
                                         ReadBufferSize, pFile_ ) );
            if ( numRead==0 ) return false;
            else if ( numRead<ReadBufferSize )
            {
                finished_=true;
                pBufMax_=buf_+numRead;
            }
            pBuf_=buf_;
        } // ~else
    } // ~if
    runLength_=lengths_[*pBuf_];
    lastChar_=codes_[*pBuf_];
#ifdef DEBUG
    cout << "Got run: " << runLength_ << " of " << lastChar_ << endl;
#endif
    pBuf_++;

    return true;

} // ~getRun

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
    assert( ( fileSize%sizeof( unsigned int ) )==0 ); // read with == 4 byte
    numInts_=fileSize/sizeof( unsigned int ); // how many ints to read from file
    //cerr << fileName << ": " << fileSize << " bytes/"<< numInts_ << " blocks" << endl;
    fseek( pFile_, 0, SEEK_SET );
    //init arrays

    for ( int i=0; i<huffmanBufferSize; i++ )
    {
        symBuf[i]=0;
        runBuf[i]=0;
    }
    soFar_.ull=0;
    toAdd_.ull=0;

    // init the token lookup table, does not need to be cleared when the file
    // is rewind since its more or less static
    // TODO: hardcode the token table? stays the same for each program call

    unsigned int codeMask;
    for ( unsigned int i( 0 ); i<numTokens; i++ )
    {
        tokenTable_[i]=0xFF;
        for ( unsigned int j( 0 ); j<numSingleCodes; j++ )
        {
            codeMask=( 1<<singleCharLength[j] )-1; // (*2^3==) -1
            if ( ( i&codeMask )==singleCharCode[j] )
            {
                assert ( tokenTable_[i]==0xFF );
                tokenTable_[i]=( j<<1 );
                //   cerr << "TT @ " << i << " is "<< itoa(tokenTable_[i],2) << endl;
            }
        } // ~for j
        for ( unsigned int j( 0 ); j<numDoubleCodes; j++ )
        {
            codeMask=( 1<<doubleCharLength[j] )-1;
            if ( ( i&codeMask )==doubleCharCode[j] )
            {
                assert ( tokenTable_[i]==0xFF );
                tokenTable_[i]=( ( j<<1 )|0x1 );
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
    runLength_=0;
    lastChar_=notInAlphabet;
    currentPos_=0;
    finished_=false;
    bitsUsed_=0;
    soFar_.ull=0;
    numSymbols_=0;
    queueCounter_=0;
    maxSymbols_=0;
    intCounter_=0;
    firstRun_ = true; // fixes 3 Bit huffman error (BTL-17)
    nearlyFinished_ = false;

    for ( int i=0; i<huffmanBufferSize; i++ )
    {
        symBuf[i]=0;
        runBuf[i]=0;
    }
} // ~rewindFile

LetterCountType BwtReaderHuffman::tellg( void ) const
{
    return currentPos_;
} // ~tellg

unsigned int BwtReaderHuffman::readAndCount( LetterCount &c, const LetterCountType numChars )
{

    LetterCountType charsLeft( numChars );
    while ( charsLeft>runLength_ )
    {
        // Below is not great design, at first call of this function it accesses an
        // out-of-range array element. Fortunately it always adds zero to it! :)
        c.count_[whichPile[lastChar_]]+=runLength_;
        charsLeft-=runLength_;
        if ( getRun()==false )
        {
            currentPos_+=( numChars-charsLeft );
            return ( numChars-charsLeft );
        } // ~if
    } // ~while

    c.count_[whichPile[lastChar_]]+=charsLeft;
    runLength_-=charsLeft;
    currentPos_+=numChars;
    return numChars;
} // ~BwtReaderHuffman::readAndCount( LetterCount& c, const int numChars )

unsigned int BwtReaderHuffman::readAndSend( BwtWriterBase &writer, const int numChars )
{
    if ( numChars == 0 )
    {
        return numChars;   // exit directy
    }

    unsigned int charsLeft( numChars );
    while ( charsLeft>runLength_ )
    {

        writer.sendRun( lastChar_, runLength_ );
        charsLeft-=runLength_;

        if ( getRun()==false )
        {
            currentPos_+=( numChars-charsLeft );
            return ( numChars-charsLeft );
        } // ~if
    } // ~while
    writer.sendRun( lastChar_, charsLeft );
    runLength_-=charsLeft;
    currentPos_+=numChars;
    return numChars;
} //~BwtReaderHuffman::readAndSend(BwtWriterBase& writer, const int numChars)


bool BwtReaderHuffman::getRun( void )
{
    numSymbols_=-1; // needed for loops

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
            toAdd_.ull=0; // init
            elementsRead = fread( &toAdd_.ui, sizeof ( unsigned int ), 1, pFile_ );

            if ( elementsRead == 1 )
            {
                toAdd_.ull <<= bitsUsed_; // left shift of how many bits used, start 0
                soFar_.ull |= toAdd_.ull; // glue both 32bit words together
                bitsUsed_ += 32; // we have succesfully read 32 bit
                intCounter_++; // and we have used one int for that
            }
            firstRun_=false;
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
int BwtReaderHuffman::operator()( char *p, int numChars )
{
    assert( 1 == 0 );
    return -1;
} // ~operator()

#endif
