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

#include "BCRext.hh"

#include "BwtReader.hh"
#include "BwtWriter.hh"
#include "LetterCount.hh"
#include "ReadBuffer.hh"
#include "SeqReader.hh"

#include "Config.hh"
#include "Filename.hh"
#include "Timer.hh"
#include "Tools.hh"
#include "Types.hh"
#include "config.h"
#include "libzoo/util/Logger.hh"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#ifdef HAVE_POSIX_FADVISE
#  include <fcntl.h>
#endif


#define BCREXT_ID "@(#) $Id: BCRext.cpp,v 1.8 2011/12/20 14:33:24 tjakobi Exp $"


using namespace std;


//typedef ReadBufferASCII ReadBuffer;
#ifdef USE_4_BITS_PER_BASE
typedef ReadBuffer4Bits ReadBuffer;
#else
#ifdef USE_PREFIX_ONLY
typedef ReadBuffer4Bits ReadBuffer;
#else
typedef ReadBufferASCII ReadBuffer;
#endif
#endif


// added by Tobias, interface to new Beetl executable
BCRext::BCRext( bool huffman, bool runlength,
                bool ascii, bool implicitSort,
                bool seqFile, const string &inFile, const string &prefix ) :

    // set tool flags
    useHuffmanEncoder_( huffman ),
    useRunlengthEncoder_( runlength ),
    useAsciiEncoder_( ascii ),
    useImplicitSort_( implicitSort ),
    useSeqFile_ ( seqFile ),
    inFile_( inFile ),
    prefix_( prefix )
{
    // Notes
    if ( implicitSort && ( huffman || runlength ) )
    {
        cout << "-> Note: -sap mode needs ASCII intermediate files," << endl
             << "-> will revert to requested compression type for final output" << endl;
    }
}

// called from beetl class, replaces main method
// args have been set through constructor before
void BCRext::run( void )
{

    Timer  timer;

    //string prefix = (string)"[" + (string)args[0] + (string)"]: ";

    string prefix = ( string )"[" + ( string )"BCRext" + ( string )"]: ";
    cerr << prefix << "time now is " << timer.timeNow();
    // Tony 13.6.12 - BCREXT_ID is not informative now we are in git world
    //  cerr << prefix << "software version is " << BCREXT_ID << endl;


    const string fileStem( "tmp1" );
    const string fileStemTemp( "tmp2" );


    string tmpIn = fileStem;
    string tmpOut = fileStemTemp;
    string tmpSwap;
    string fileName;

    // output streams for sequences - 1 per pile
    vector <FILE *> outSeq( alphabetSize );
    // output streams for positions of suffixes in pile - 1 per pile
    vector <FILE *> outPtr( alphabetSize );
    // output stream for updated BWTs - 1 per pile
    vector <BwtWriterBase *> outBwt( alphabetSize );
    // output streams for original read numberings - 1 per pile
    vector <FILE *> outNum( alphabetSize );

    // input streams for BWT from previous iter - 1 per pile
    // this one used to compute the counts to do backward search
    vector <BwtReaderBase *> inBwt( alphabetSize );

    // input streams for BWT from previous iter - 1 per pile
    // this one used to read the BWT chunk by chunk to allow new
    // chars to be interspersed
    vector <BwtReaderBase *> inBwt2( alphabetSize );

    //  FILE* inPtr;
    //  FILE* inNum;
    BwtWriterBase *outDollarBwt;

    vector<char> bwtBuf;
    // extra byte accounts for a fact that inserted BWT characters
    // are appended to the buffer after the interspersed characters so as to
    // make a single write to file
    bwtBuf.resize( bwtBufferSize + 1 );

    if (    whichPile[( int ) '$'] != 0 ||
            whichPile[( int ) 'A'] != 1 ||
            whichPile[( int ) 'C'] != 2 ||
            whichPile[( int ) 'G'] != 3 ||
            whichPile[( int ) alphabet[4]] != 4 ||
            whichPile[( int ) alphabet[5]] != 5
#ifdef USE_EXTRA_CHARACTER_Z
            || whichPile[( int ) notInAlphabet] != 6
#endif
       )
    {
        cerr << "Something seems to be wrong with the alphabet table!" << endl;
        exit( EXIT_FAILURE );
    }


    Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << prefix << "Using alphabet = " << alphabet << ", size = " << alphabetSize << endl;
    //  cerr << BUFSIZ << endl

    if ( alphabetSize >= 10 )
    {
        cerr << "Alphabet sizes larger than 9 are not supported yet." << endl;
        exit( EXIT_FAILURE );
    }

    SequenceNumber seqNum( 0 );


    const LetterNumber sameAsPrevFlag( ( ( LetterNumber )1 ) << ( ( 8 * sizeof( LetterNumber ) ) - 1 ) );
    const LetterNumber sameAsPrevMask( ~sameAsPrevFlag );

    //  cout << sameAsPrevFlag << " "<< sameAsPrevMask << " "<< (sameAsPrevMask&sameAsPrevFlag) << " " << (((LetterNumber)1)<<63) << endl;

    LetterNumber seqPtr;
    int thisPile, lastPile;
    LetterNumber posInPile;
    //  char inChar;
    LetterNumber charsToGrab;// charsLeft, charsThisBatch;


    Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << prefix << "Will read sequences from file " << inFile_ << endl;

    // read first sequence to determine read size
    SeqReaderFile *seqReader( SeqReaderFile::getReader( fopen( inFile_.c_str(), "rb" ) ) );
    const char *seqBuf = seqReader->thisSeq();

    const int seqSize( strlen( seqBuf ) - 1 ); // -1 compensates for \n at end
    Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << prefix << "Assuming all sequences are of length " << seqSize << endl;
    //  inFile.seekg(0,ios::beg);
    //  rewind(inSeq);

    if ( ( seqSize % 2 ) == 1 ) // if odd
    {
        //    cout << "ODD" << endl;
        tmpIn = fileStem;
        tmpOut = fileStemTemp;
    } // ~if
    else
    {
        //    cout << "EVEN" << endl;
        tmpIn = fileStemTemp;
        tmpOut = fileStem;
    } // ~else

    fileName = TmpFilename( fileStem, "-B0", 0 ).str();
    readWriteCheck( fileName.c_str(), true );
    //outDollarBwt = fopen( fileName.c_str(), "w" );
    if ( useImplicitSort_ || useAsciiEncoder_ )
        outDollarBwt = new BwtWriterASCII( fileName.c_str() );
#ifdef ACTIVATE_HUFFMAN
    else if ( useHuffmanEncoder_ )
        outDollarBwt = new BwtWriterHuffman( fileName.c_str() );
#endif
    else if ( useRunlengthEncoder_ )
        outDollarBwt = new BwtWriterRunLengthV3( fileName.c_str() );
    else
        assert( false );

    for ( int j( 1 ); j < alphabetSize; j++ )
    {

        fileName = TmpFilename( tmpIn, "-S0", j ).str();
        readWriteCheck( fileName.c_str(), true );
        outSeq[j] = fopen( fileName.c_str(), "w" );


        fileName = TmpFilename( tmpIn, "-P0", j ).str();
        readWriteCheck( fileName.c_str(), true );
        outPtr[j] = fopen( fileName.c_str(), "w" );


#ifdef TRACK_SEQUENCE_NUMBER
        fileName = TmpFilename( tmpIn, "-N0", j ).str();
        readWriteCheck( fileName.c_str(), true );
        outNum[j] = fopen( fileName.c_str(), "w" );
#endif

        fileName = TmpFilename( tmpIn, "-B0", j ).str();

        if ( useImplicitSort_ || useAsciiEncoder_ )
            outBwt[j] = new BwtWriterASCII( fileName.c_str() );
#ifdef ACTIVATE_HUFFMAN
        else if ( useHuffmanEncoder_ )
            outBwt[j] = new BwtWriterHuffman( fileName.c_str() );
#endif
        else if ( useRunlengthEncoder_ )
            outBwt[j] = new BwtWriterRunLengthV3( fileName.c_str() );
        else
            assert( false );

    } // ~for

    LetterCount dollars;
    //  vector<LetterCount> alreadyInPile(alphabetSize);
    //  vector<LetterCount> countedThisIter(alphabetSize);
    LetterCountEachPile alreadyInPile;
    LetterCountEachPile countedThisIter;
    LetterCountEachPile newCharsThisIter;

    // TBD Rationalize count names wrt first and subsequent iterations
    LetterCount addedSoFar;
    LetterCount outputSoFar;

    LetterCount prevCharsOutputThisIter;
    LetterCount newCharsAddedThisIter;


    //  LetterCount readSoFar[alphabetSize];

    // First iteration
    // - move original sequence into piles based on last character
    // - work out BWT corresponding to 0-suffixes and 1-suffixes
    // TBD check for invalid chars, do qual masking

    ReadBuffer readBuffer( seqSize, -1, -1, -1 );

    if ( readBuffer.blockSize_ <= seqSize + 1 )
    {
        cerr << "ReadBuffer blocksize is too small (" << readBuffer.blockSize_ << "). Aborting." << endl;
        exit( EXIT_FAILURE );
    }

    // copy first sequence over
    strcpy( readBuffer.seqBufBase_, seqBuf );
    do
    {
        thisPile = whichPile[( int )readBuffer.seqBufBase_[seqSize - 1]];

        // zero is terminator so should not be present
        if ( thisPile < 0 )
        {
            cerr << "Pile must not be < 0. Aborting. At char |" << readBuffer.seqBufBase_[seqSize - 1] << "|" << endl;
            exit( EXIT_FAILURE );
        }

        if ( thisPile > alphabetSize )
        {
            cerr << "Pile must not be > alphabet size. Aborting." << endl;
            exit( EXIT_FAILURE );
        }

#ifdef DEBUG
        cout << readBuffer.seqBufBase_ << endl;
#endif

        // count characters and output first N chars of BWT
        // (those preceding terminator characters)
        dollars.count_[thisPile]++;

        /*
                if ( fwrite( readBuffer.seqBufBase_ + seqSize - 1, sizeof ( char ), 1, outDollarBwt ) != 1 )
                {
                    cerr << "Could not write to Dollar Pile. Aborting." << endl;
                    exit( EXIT_FAILURE );
                }
        */
        ( *outDollarBwt )( readBuffer.seqBufBase_ + seqSize - 1, 1 );

        //    fprintf( outSeq[thisPile], "%s", readBuffer.seqBufBase_);
        readBuffer.convertFromASCII();
        readBuffer.sendTo( outSeq[thisPile] );

#ifdef TRACK_SEQUENCE_NUMBER
        assert( fwrite( &seqNum, sizeof( SequenceNumber ),
                        1, outNum[thisPile] ) == 1 );
        //    seqNum++;
#endif

        // create BWT corresponding to 1-suffixes

        if ( whichPile[( int )readBuffer.seqBufBase_[seqSize - 2]] < 0 ||
             whichPile[( int )readBuffer.seqBufBase_[seqSize - 2]] > alphabetSize  )
        {
            cerr << "Trying to write non alphabet character to pile. Aborting." << endl;
            exit( EXIT_FAILURE );
        }

        countedThisIter[thisPile].count_[whichPile[( int )readBuffer.seqBufBase_[seqSize - 2]]]++;
        //    assert(fwrite( readBuffer.seqBufBase_+seqSize-2, sizeof(char), 1, outBwt[thisPile] )==1);

        seqPtr = *( addedSoFar.count_ + thisPile );

        if ( useImplicitSort_ && ( addedSoFar.count_[thisPile] != 0 ) )
        {
            //      cout << thisPile << " " << addedSoFar.count_[thisPile] << " 1\n";
            seqPtr |= sameAsPrevFlag; // TBD replace if clause with sum
            //      *(readBuffer.seqBufBase_+seqSize-2)+=32;//tolower(*(readBuffer.seqBufBase_+seqSize-2));
            *( readBuffer.seqBufBase_ + seqSize - 2 ) = tolower( *( readBuffer.seqBufBase_ + seqSize - 2 ) );
        }

        ( *outBwt[thisPile] )( readBuffer.seqBufBase_ + seqSize - 2, 1 );


        if ( fwrite( &seqPtr, sizeof( LetterNumber ),
                     1, outPtr[thisPile] ) != 1 )
        {
            cerr << "Could not write to pointer pile. Aborting." << endl;
            exit( EXIT_FAILURE );
        }
        addedSoFar.count_[thisPile]++;
        seqNum++;

        seqReader->readNext( readBuffer.seqBufBase_ );
    } // ~while
    while ( !seqReader->allRead() );

    delete outDollarBwt;

    Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << prefix << "Read " << seqNum << " sequences" << endl;
    for ( int i( 1 ); i < alphabetSize; i++ )
    {
        fclose( outSeq[i] );
        fclose( outPtr[i] );
#ifdef TRACK_SEQUENCE_NUMBER
        fclose( outNum[i] );
#endif
        delete outBwt[i];
        //    fclose(outBwt[i]);
    }
    //    return (0);

    LetterCount lastSAPInterval;
    LetterNumber thisSAPInterval;

    //  ReadBuffer buffer(seqSize);

    // Main loop
    for ( int i( 2 ); i <= seqSize; i++ )
    {
        thisSAPInterval = 0;
        lastSAPInterval.clear();

        cout << "Starting iteration " << i << ", time now: " << timer.timeNow();
        cout << "Starting iteration " << i << ", usage: " << timer << endl;

        // don't do j=0 - this is the $ sign which is done already
        for ( int j( 1 ); j < alphabetSize; j++ )
        {
            // prep the output files

            fileName = TmpFilename( tmpOut, "-S0", j ).str();
            readWriteCheck( fileName.c_str(), true );
            outSeq[j] = fopen( fileName.c_str(), "w" );

            fileName = TmpFilename( tmpOut, "-P0", j ).str();
            readWriteCheck( fileName.c_str(), true );
            outPtr[j] = fopen( fileName.c_str(), "w" );

#ifdef TRACK_SEQUENCE_NUMBER
            fileName = TmpFilename( tmpOut, "-N0", j ).str();
            readWriteCheck( fileName.c_str(), true );
            outNum[j] = fopen( fileName.c_str(), "w" );
#endif

            fileName = TmpFilename( tmpOut, "-B0", j ).str();




            if ( ( useImplicitSort_ && ( i != seqSize ) ) || useAsciiEncoder_ )
                outBwt[j] = new BwtWriterASCII( fileName.c_str() );
#ifdef ACTIVATE_HUFFMAN
            else if ( useHuffmanEncoder_ )
                outBwt[j] = new BwtWriterHuffman( fileName.c_str() );
#endif
            else if ( useRunlengthEncoder_ )
                outBwt[j] = new BwtWriterRunLengthV3( fileName.c_str() );
            else
                assert( false );

            if ( useImplicitSort_ && ( i == seqSize ) )
            {
                BwtWriterBase *p( new BwtWriterImplicit( outBwt[j] ) );
                outBwt[j] = p; // ... and the deception is complete!!!
            } // ~if


#ifdef DEBUG
            cout << "Prepping output file " << tmpOut << endl;
#endif

            setvbuf( outSeq[j], NULL, _IOFBF, 262144 );
            //  setvbuf( outPtr[j], NULL, _IOFBF, 65536);
            //  setvbuf( outNum[j], NULL, _IOFBF, 65536);
            //  setvbuf( outBwt[j], NULL, _IOFBF, 65536);

            // prep the input files
            fileName = TmpFilename( tmpIn, "-B0", j ).str();
            // select the proper input module
            if ( useImplicitSort_ || useAsciiEncoder_ )
            {
                inBwt[j] = new BwtReaderASCII( fileName.c_str() );
                inBwt2[j] = new BwtReaderASCII( fileName.c_str() );
            }
#ifdef ACTIVATE_HUFFMAN
            else if ( useHuffmanEncoder_ )
            {
                inBwt[j] = new BwtReaderHuffman( fileName.c_str() );
                inBwt2[j] = new BwtReaderHuffman( fileName.c_str() );
            }
#endif
            else if ( useRunlengthEncoder_ )
            {
                inBwt[j] = new BwtReaderRunLengthV3( fileName.c_str() );
                inBwt2[j] = new BwtReaderRunLengthV3( fileName.c_str() );
            }
            else
                assert( false );

#ifdef DEBUG
            cout << "Prepping input file " << tmpIn << endl;
#endif

        } // ~for j

        addedSoFar.clear();
        outputSoFar.clear();

        prevCharsOutputThisIter.clear();
        newCharsAddedThisIter.clear();

#ifdef DEBUG
        cout << "already in pile" << endl;
        alreadyInPile.print();
        cout << "counted this iter" << endl;
        countedThisIter.print();
#endif
        countedThisIter.clear();
        newCharsThisIter.clear();

#ifdef DEBUG
        cout << "Count in dollars pile: ";
        dollars.print();
#endif

        // don't do j=0; $ sign done already
        for ( int j( 1 ); j < alphabetSize; j++ )
        {
            int fdSeq, fdNum, fdPtr;

#ifndef TRACK_SEQUENCE_NUMBER
            fdNum = 0;
#endif

            // read each input file in turn

            fileName = TmpFilename( tmpIn, "-S0", j ).str();
            readWriteCheck( fileName.c_str(), false );
            fdSeq = open( fileName.c_str(), O_RDONLY, 0 );

            fileName = TmpFilename( tmpIn, "-P0", j ).str();
            readWriteCheck( fileName.c_str(), false );
            fdPtr = open( fileName.c_str(), O_RDONLY, 0 );


#ifdef TRACK_SEQUENCE_NUMBER
            fileName = TmpFilename( tmpIn, "-N0", j ).str();
            readWriteCheck( fileName.c_str(), false );
            fdNum = open( fileName.c_str(), O_RDONLY, 0 );

#ifdef HAVE_POSIX_FADVISE
            assert( posix_fadvise( fdNum, 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_NOREUSE | POSIX_FADV_WILLNEED ) != -1 );
#endif
#endif

#ifdef HAVE_POSIX_FADVISE
            assert( posix_fadvise( fdSeq, 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_NOREUSE | POSIX_FADV_WILLNEED ) != -1 );
            assert( posix_fadvise( fdPtr, 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_NOREUSE | POSIX_FADV_WILLNEED ) != -1 );
#endif

#ifdef USE_PREFIX_ONLY
            ReadBufferPrefix buffer( seqSize, i, fdSeq, fdNum, fdPtr );
#else
            ReadBuffer buffer( seqSize, fdSeq, fdNum, fdPtr );
#endif

            while ( buffer.getNext( seqNum, seqPtr ) )
            {
                bool thisSAPValue = ( ( seqPtr & sameAsPrevFlag ) != 0 );
                if ( thisSAPValue )
                {
                    seqPtr &= sameAsPrevMask;
                }
                else
                {
                    thisSAPInterval++;
                }

                thisPile = buffer[seqSize - i];

                //thisPile=whichPile[seqBuff[seqSize-i]];

                if ( thisPile < 0 )
                {
                    cerr << "Pile must not be < 0. Aborting." << endl;
                    exit( EXIT_FAILURE );
                }
                lastPile = buffer[seqSize - i + 1];

                //lastPile=whichPile[seqBuff[seqSize-i+1]];
                if ( lastPile < 0 )
                {
                    cerr << "Pile must not be < 0. "  << endl;
                    exit( EXIT_FAILURE );
                }



#ifdef DEBUG
                cout << ( ( thisSAPValue ) ? '1' : '0' ) << " " << thisSAPInterval << " " << seqPtr << " " << seqNum << " " << thisPile << " " << lastPile << endl;
                cout << "Read in " << seqPtr << " " << seqNum << " " << thisPile << " " << lastPile << endl;
                for ( int ZZ( 0 ); ZZ < seqSize; ZZ++ )
                {
                    cout << "ZZ " << ZZ << endl;
                    cout << alphabet[buffer[ZZ]] << endl;
                }
                cout << endl;
#endif

                // *** work out position in new pile ***

                // sum contents of lexicographically smaller piles
                // ... probably possible to speed this up by storing cumulative sums
#ifdef DEBUG
                cout << "already in pile" << endl;
                alreadyInPile.print();
#endif

                // posInPile=0;
                posInPile = dollars.count_[thisPile];
                //   cout << posInPile << " " << thisPile << " " << lastPile << endl;
                for ( int k( 1 ); k < lastPile; k++ )
                {
                    posInPile += alreadyInPile[k].count_[thisPile];
                    //   cout << posInPile << endl;
                } // ~for k

#ifdef DEBUG
                cout << "posInPile starts at " << posInPile << endl;
                cout << "counting in pile " << alphabet[lastPile] << endl;
#endif

                // count all chars prior to seqPtr in lastPile
                // read seqPtr too, but don't add to countedThisIter
                charsToGrab = seqPtr - addedSoFar.count_[lastPile]; //+1;
#ifdef DEBUG
                cout << "charsToGrab " << charsToGrab << endl;
#endif

                // Should now always read at least 1 byte
                if ( charsToGrab < 0 )
                {
                    cerr << "Tried to grap < 0 chars. Aborting." << endl;
                    exit( EXIT_FAILURE );
                }

                LetterNumber readCountChars = inBwt[lastPile]->readAndCount
                                              ( countedThisIter[lastPile], charsToGrab );

                if ( readCountChars != charsToGrab )
                {
                    cerr << "BWT readAndCount returned only " << readCountChars
                         << " chars. Expected " << charsToGrab
                         << " chars. Aborting." << endl;
                    exit( EXIT_FAILURE );
                }

                inBwt[lastPile]->readAndCount( newCharsThisIter[lastPile], 1 );


                addedSoFar.count_[lastPile] = seqPtr + 1;


                posInPile += countedThisIter[lastPile].count_[thisPile];



#ifdef DEBUG
                cout << "counted this iter" << endl;
                countedThisIter.print();
#endif

                // *** add char into new pile ***

                // read and output bytes up to insertion point

                charsToGrab = posInPile - prevCharsOutputThisIter.count_[thisPile];

                LetterNumber readSendChars = inBwt2[thisPile]->readAndSend
                                             ( *outBwt[thisPile], charsToGrab );

                if ( readSendChars != charsToGrab )
                {
                    cerr << "BWT readAndSend returned only " << readSendChars
                         << " chars. Expected " << charsToGrab
                         << " chars. Aborting." << endl;
                    exit( EXIT_FAILURE );
                }

                // bwtBuf[0]=(seqSize-i-1>=0)?baseNames[buffer[seqSize-i-1]]:'$';
                bwtBuf[0] = ( seqSize - i - 1 >= 0 ) ? alphabet[buffer[seqSize - i - 1]] : alphabet[0];
                // if (thisSAPValue==true) bwtBuf[0]+=32;//=tolower(bwtBuf[0]);



                prevCharsOutputThisIter.count_[thisPile] = posInPile;

                // pointer into new pile must be offset by number of new entries added
                // seqPtr+=newCharsAddedThisIter.count_[thisPile];
                seqPtr = posInPile + newCharsAddedThisIter.count_[thisPile];

                if ( useImplicitSort_ )
                {
                    if ( lastSAPInterval.count_[thisPile] == thisSAPInterval )
                    {
                        bwtBuf[0] = tolower( bwtBuf[0] );
                        //     bwtBuf[0]+=32;
                        seqPtr |= sameAsPrevFlag;
                        //   cout << thisSAPInterval << endl;
                    }
                    else
                    {
                        //   cout << thisSAPInterval << " " << lastSAPInterval.count_[thisPile] << endl;
                        lastSAPInterval.count_[thisPile] = thisSAPInterval;
                    }
                }

                ( *outBwt[thisPile] )( bwtBuf.data(), 1 );

                if ( fwrite( &seqPtr, sizeof ( LetterNumber ),
                             1, outPtr[thisPile] ) != 1 )
                {
                    cerr << "BWT readAndSend returned only " << readSendChars
                         << " chars. Expected " << charsToGrab
                         << " chars. Aborting." << endl;
                    exit( EXIT_FAILURE );
                }


#ifdef DEBUG
                cout << "adding pointer " << seqPtr << " to pile "
                     << alphabet[thisPile] << endl;
#endif
                // then the offset itself is updated
                newCharsAddedThisIter.count_[thisPile]++;

                // do radix sort
                // fprintf( outSeq[thisPile], "%s\n", seqBuff);
                // assert(fwrite( seqBuff, sizeof(char),
                //       1+seqSize, outSeq[thisPile] )==1+seqSize);
                buffer.sendTo( outSeq[thisPile] );


#ifdef TRACK_SEQUENCE_NUMBER
                assert( fwrite( &seqNum, sizeof( SequenceNumber ),
                                1, outNum[thisPile] ) == 1 );
#endif

            } // ~while


            close( fdSeq );
            close( fdPtr );
#ifdef TRACK_SEQUENCE_NUMBER
            close( fdNum );
#endif

        } // ~for j

        Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "All new characters inserted, usage: " << timer << endl;
        for ( int j( 1 ); j < alphabetSize; j++ )
        {

            while ( inBwt[j]->readAndCount
                    ( countedThisIter[j], ReadBufferSize ) == ReadBufferSize );

        } // ~for j

        Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "finishing off BWT strings" << endl;
        for ( int j( 1 ); j < alphabetSize; j++ )
        {

            while ( inBwt2[j]->readAndSend( *outBwt[j], ReadBufferSize )
                    == ReadBufferSize );

        } // ~for

#ifdef DEBUG
        cout << "final value of counted this iter" << endl;
        countedThisIter.print();
        cout << "final value of new chars this iter" << endl;
        newCharsThisIter.print();
#endif

        alreadyInPile += newCharsThisIter;

        for ( int j( 1 ); j < alphabetSize; j++ )
        {
            fclose( outSeq[j] );
            fclose( outPtr[j] );
#ifdef TRACK_SEQUENCE_NUMBER
            fclose( outNum[j] );
#endif
            delete ( outBwt[j] );
            delete ( inBwt[j] );
            delete ( inBwt2[j] );
        } // ~for j

        tmpSwap = tmpIn;
        tmpIn = tmpOut;
        tmpOut = tmpSwap;

        Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "finished iteration " << i  << ", usage: " << timer << endl;

        //    assert(i<2);
    } // ~for i (main iteration)

    string fileTypes( "BPS" );
#ifdef REMOVE_TEMPORARY_FILES
    for ( int j( 1 ); j < alphabetSize; j++ )
    {
        for ( unsigned int i( 0 ); i < fileTypes.size(); i++ )
        {
            char s[4] = "-B0";
            s[1] = fileTypes[i];
            fileName = TmpFilename( tmpOut, s, j ).str();
            if ( remove( fileName.c_str() ) != 0 )
            {
                cerr << "Warning: failed to clean up temporary file " << fileName
                     << endl;
            } // ~if
            else
            {
                Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Removed temporary file " << fileName << endl;
            } // ~if
        } // ~for i
    } // ~for j
#endif

    // Move files to final output directory
    for ( int j( 0 ); j < alphabetSize; j++ )
    {
        for ( unsigned int i( 0 ); i < fileTypes.size(); i++ )
        {
            if ( j == 0 && i >= 1 ) continue; // "-S00" and "-P00" don't exist
            char s[4] = "-B0";
            s[1] = fileTypes[i];
            fileName = TmpFilename( fileStem, s, j ).str();
            Filename filenameOut( prefix_, s, j );
            if ( safeRename( fileName, filenameOut.str() ) != 0 )
            {
                cerr << "Warning: failed to rename temporary file \"" << fileName << "\" to final output " << endl;
            }
        }
    }

    Logger::out() << "Final output files are named " << prefix_ << "-Bxx and similar" << endl;
}


