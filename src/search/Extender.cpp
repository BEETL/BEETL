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

#include "Extender.hh"

#include "BCRext.hh"
#include "BCRexternalBWT.hh"
#include "Timer.hh"
#include "config.h"
#include "libzoo/util/Logger.hh"
#include "libzoo/util/TemporaryFilesManager.hh"

#include <cstring>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;


vector<IntervalRecord> inputIntervals;


Extender::Extender(
    const ExtendParameters &extendParams
)
    : extendParams_( extendParams )
{
    // input checks: BWT-end-pos file needed if sequence numbers output requested
    if ( extendParams_["sequence numbers output filename"].isSet() )
    {
        string bwtPrefix = extendParams_["bwt filename prefix"];
        string endPosFilename = bwtPrefix + "-end-pos";
        if ( !readWriteCheck( endPosFilename.c_str(), false, false ) )
        {
            cerr << "Error: " << endPosFilename << " is needed to extract sequence numbers (it can be generated during BWT construction with beetl-bwt --generate-end-pos-file)." << endl;
            exit( 1 );
        }
    }
}

void Extender::fillRangeStore( RangeStoreExternal &rangeStore, const LetterCountEachPile &countsPerPile, const LetterCountEachPile &countsCumulative )
{
    string intervalsFileName = extendParams_.getStringValue( "intervals filename" );
    ifstream intervalsFile( intervalsFileName.c_str() );
    IntervalReader intervalReader( intervalsFile );

    inputIntervals = intervalReader.readFullFileAsVector();
    std::sort( inputIntervals.begin(), inputIntervals.end(), IntervalRecord::bwtPositionCompare );

    LetterNumber minPos = 0;
    for ( IntervalRecord & rec : inputIntervals )
    {
        Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Processing interval: " << rec << endl;

        LetterNumber pos = rec.position;
        LetterNumber count = rec.count;
        if ( pos < minPos )
        {
            LetterNumber diff = minPos - pos;
            if ( count > diff )
            {
                pos = minPos;
                count -= diff;
            }
            else
                count = 0;
        }

        if ( count )
        {
            LetterNumber start = 0, end = 0;
            for ( int i( 0 ); i < alphabetSize; i++ )
            {
                start = end;
                end += countsPerPile[whichPile[( int )rec.kmer[0]]].count_[i];

                LetterNumber intersectionStart = max( pos, start );
                LetterNumber intersectionEnd = min( pos + count, end );

                if ( intersectionStart < intersectionEnd )
                {
                    rangeStore.addRange(
                        Range(
                            rec.kmer,
                            intersectionStart,
                            intersectionEnd - intersectionStart,
                            false, // unused
                            true, // hasUserData
                            &rec
                        ),
                        whichPile[( int )rec.kmer[0]],
                        i,
                        "", //subset
                        1 //initial cycle
                    );
                }
            }
            minPos = pos + count;
        }
    }
}

void Extender::run()
{
    Timer  timer;
    bool compressIntermediateBwts = true;

    vector <BwtReaderBase *> inBwt( alphabetSize );

    LetterCountEachPile countsPerPile, countsCumulative;

    string subset_ = "";

    for ( int i( 0 ); i < alphabetSize; i++ )
    {
        stringstream fileNameSS;
        fileNameSS << extendParams_.getStringValue( "bwt filename prefix" ) << "-B0" << i;
        string fileName = fileNameSS.str().c_str();
        if ( compressIntermediateBwts == true )
        {
            if ( extendParams_["use indexing"].isSet() )
            {
                cout << "Using indexed BWT file" << endl;
                inBwt[i] = new BwtReaderRunLengthIndex( fileName );
            }
            else
            {
                cout << "Using non-indexed BWT file" << endl;
                inBwt[i] = new BwtReaderRunLength( fileName );
            }
        }
        else
            inBwt[i] = new BwtReaderASCII( fileName );

        inBwt[i]->readAndCount( countsPerPile[i] );
    }

    countsCumulative = countsPerPile;

    for ( int i( 1 ); i < alphabetSize; i++ )
        countsCumulative[i] += countsCumulative[i - 1];

    cout << "countsPerPile:" << endl;
    countsPerPile.print();
    cout << "countsCumulative:" << endl;
    countsCumulative.print();

    RangeStoreExternal r;
    fillRangeStore( r, countsPerPile, countsCumulative );


    string currentWord;

    LetterCount countsSoFar;
    LetterNumber currentPos;
    LetterNumber numRanges, numSingletonRanges;

    r.clear();
    for ( int cycle( 1 ); ; ++cycle )
    {
        cout << "cycle: " << cycle << endl;

        Logger_if( LOG_SHOW_IF_VERBOSE )
        {
            Logger::out() << "   time now: " << timer.timeNow();
            Logger::out() << "   usage: " << timer << endl;
        }

#ifdef PROPAGATE_SEQUENCE
        currentWord.resize( cycle + 2 );
#endif

        numRanges = 0;
        numSingletonRanges = 0;
        r.setCycleNum( cycle );

        for ( int i( 1 ); i < alphabetSize; ++i )
        {
            inBwt[i]->rewindFile();
            currentPos = 0;
            countsSoFar.clear();
            countsSoFar += countsCumulative[i - 1];

            for ( int j( 1 ); j < alphabetSize; ++j )
            {
                r.setPortion( i, j );

                OneBwtBackTracker backTracker(
                    inBwt[i],
                    currentPos,
                    r,
                    countsSoFar,
                    subset_,
                    cycle + 1,
                    true, // propagate breakpoints until $ sign
                    true // skip-already-processed-intervals deactivated
                );

                ExtenderIntervalHandler intervalHandler;
                Range rangeObject( true );
                backTracker.process( i, currentWord, intervalHandler, rangeObject );

                numRanges += backTracker.numRanges_;
                numSingletonRanges += backTracker.numSingletonRanges_;

            } // ~for j
            //cerr << "Done i " << i <<endl;
        }     // ~for i
        r.clear();

        Logger_if( LOG_SHOW_IF_VERBOSE )
        {
            Logger::out() << "Finished cycle " << cycle << ": ranges=" << numRanges << " singletons=" << numSingletonRanges << endl;
        }

        if ( numRanges == 0 ) break;

    } // ~for c

    for ( int i = 0; i < alphabetSize; i++ )
        delete inBwt[i];

    // Dollar BWT positions output
    if ( extendParams_["dollar positions output filename"].isSet() )
    {
        string outputFilename = extendParams_["dollar positions output filename"];
        ofstream ofs( outputFilename.c_str() );
        IntervalWriter writer( ofs );
        for ( auto rec : inputIntervals )
        {
            writer.writeV2( rec );
        }
    }

    // Sequence numbers output
    if ( extendParams_["sequence numbers output filename"].isSet() )
    {
        string endPosFilename = extendParams_.getStringValue( "bwt filename prefix" ) + "-end-pos";
        ifstream endPosFile( endPosFilename.c_str() );
        SequenceNumber numDollarEntries = 0;
        endPosFile.read( reinterpret_cast< char * >( &numDollarEntries ), sizeof( SequenceNumber ) );
        uint8_t subSequenceCount = 0;
        endPosFile.read( reinterpret_cast< char * >( &subSequenceCount ), sizeof( uint8_t ) );
        uint8_t hasRevComp = 0;
        endPosFile.read( reinterpret_cast< char * >( &hasRevComp ), sizeof( uint8_t ) );
        assert( endPosFile.good() );
        numDollarEntries *= subSequenceCount * ( hasRevComp ? 2 : 1 );

        string outputFilename = extendParams_["sequence numbers output filename"];
        ofstream ofs( outputFilename.c_str() );
        IntervalWriter writer( ofs );
        for ( auto rec : inputIntervals )
        {
            for ( auto dollarPos : rec.dollarSignPositions )
            {
                //assert( dollarPos < numDollarEntries );
                if ( dollarPos >= numDollarEntries )
                {
                    cout << "Warning: dollarPos " << dollarPos << " >= numDollarEntries " << numDollarEntries << endl;
                    //                    continue;
                    dollarPos %= numDollarEntries;
                }
                endPosFile.seekg( sizeof( SequenceNumber ) + 2 * sizeof( uint8_t ) + ( dollarPos ) * ( sizeof( SequenceNumber ) + sizeof( uint8_t ) ) );

                SequenceNumber seqN;
                endPosFile.read( reinterpret_cast< char * >( &seqN ), sizeof( SequenceNumber ) );
                //                LetterNumber posN;
                //                endPosFile.read( reinterpret_cast< char * >( &posN ), sizeof( LetterNumber ) );
                //                AlphabetSymbol pileN;
                //                endPosFile.read( reinterpret_cast< char * >( &pileN ), sizeof( AlphabetSymbol ) );
                uint8_t subSequenceNum;
                endPosFile.read( reinterpret_cast< char * >( &subSequenceNum ), sizeof( uint8_t ) );

                ofs << seqN << " # " << rec.kmer << " (subSequence " << ( int )subSequenceNum << ")\n";
            }
        }
    }
}
