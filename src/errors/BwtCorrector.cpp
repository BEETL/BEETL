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

#include "BwtCorrector.hh"

#include "BCRext.hh"
#include "BCRexternalBWT.hh"
#include "Timer.hh"
#include "config.h"
#include "parameters/BwtParameters.hh"
#include "libzoo/util/Logger.hh"
#include "libzoo/util/TemporaryFilesManager.hh"

#include <cstring>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;
using namespace BeetlBwtParameters;

int BwtCorrector::getMinSupport( int cycle )
{
    //Decide on minimum support to use based on cycle.
    //Always use "min support" parameter if it is set,
    //otherwise re-calculate based on cycle as done in HiTEC algorithm.
    if ( minSupport_ == 0 )
    {
        HiTECStats stats(
            errorRate_,
            genomeLength_,
            numberOfReads_,
            readLength_
        );

        return stats.CalculateSupport(
                   max<int>( stats.Calculate_wM(), cycle + 1 )
               );
    }
    else
        return minSupport_;
}

ErrorStore BwtCorrector::findErrors()
{
    const bool propagateSequence = false;//( compareParams_ ? ( *compareParams_ )["propagate sequence"] : false );
    ErrorStore result;

    Timer  timer;
    vector <BwtReaderBase *> inBwt( alphabetSize );

    LetterCountEachPile countsPerPile, countsCumulative;

    RangeStoreExternal r;

    int numCycles( readLength_ );
    //    int minOcc( numberOfReads_ );

    inBwt = instantiateBwtPileReaders( indexPrefix_, correctorParams_->getStringValue( "use shm" ) );

    for ( int i( 0 ); i < alphabetSize; i++ )
    {
        inBwt[i]->  readAndCount( countsPerPile[i] );
    }
    countsCumulative = countsPerPile;
    for ( int i( 1 ); i < alphabetSize; i++ )
        countsCumulative[i] += countsCumulative[i - 1];

    countsCumulative.print();

    string currentWord = "xx";
    for ( int i( 1 ); i < alphabetSize; ++i )
    {
        if ( propagateSequence )
            currentWord[1] = alphabet[i];
        for ( int j( 1 ); j < alphabetSize; ++j )
        {
            if ( propagateSequence )
                currentWord[0] = alphabet[j];

            if ( countsPerPile[i].count_[j] != 0 )
                r.addRange(
                    ErrorCorrectionRange(
                        currentWord,
                        countsCumulative[i - 1].count_[j],
                        countsPerPile[i].count_[j],
                        false
                    )
                    , j,
                    i,
                    subset_,
                    1
                );
        } // ~for j
    } // ~for i

    LetterCount countsSoFar;
    LetterNumber currentPos;
    LetterNumber numRanges, numSingletonRanges;

    r.clear();
    for ( int c( 0 ); c < numCycles; ++c )
    {
        int minimumSupport = getMinSupport( c );

        cout << "cycle: " << c << endl;

        Logger_if( LOG_SHOW_IF_VERBOSE )
        {
            Logger::out() << "   time now: " << timer.timeNow();
            Logger::out() << "   usage: " << timer << endl;
        }

        numRanges = 0;
        numSingletonRanges = 0;
        r.setCycleNum( c + 1 );

        for ( int i( 1 ); i < alphabetSize; ++i )
        {
            string thisWord( c + 3, 'x' );

            inBwt[i]->rewindFile();
            currentPos = 0;
            countsSoFar.clear();
            countsSoFar += countsCumulative[i - 1];
#ifdef PROB_NOT_NEEDED
            BwtReaderRunLengthIndex *pRun;

            pRun = dynamic_cast<BwtReaderRunLengthIndex *>( inBwt[i] );
            if ( pRun != NULL )
                pRun->initIndex( countsSoFar );
            else
                inBwtA[i]->rewindFile();
#endif

            for ( int j( 1 ); j < alphabetSize; ++j )
            {
                r.setPortion( i, j );

                OneBwtBackTracker backTracker(
                    inBwt[i],
                    currentPos,
                    r,
                    countsSoFar,
                    //                    numCycles,
                    subset_,
                    c + 2,
                    false, // doesn't propagate to end of reads
                    true, // skip-already-processed-intervals deactivated
                    propagateSequence,
                    endPosFile_
                );

                BwtCorrectorIntervalHandler intervalHandler( result, minWitnessLength_, minimumSupport, c + 2 );
                ErrorCorrectionRange rangeObject;
                backTracker.process( i, thisWord, intervalHandler, rangeObject );

                numRanges += backTracker.numRanges_;
                numSingletonRanges += backTracker.numSingletonRanges_;

            } // ~for j
            //cerr << "Done i " << i <<endl;
        }     // ~for i
        r.clear();
        //    return 0; // %%%
        Logger_if( LOG_SHOW_IF_VERBOSE )
        {
            Logger::out() << "Finished cycle " << c << ": ranges=" << numRanges << " singletons=" << numSingletonRanges << endl;
        }

        if ( numRanges == 0 ) break;

    } // ~for c
    for ( int i = 0; i < alphabetSize; i++ )
        delete inBwt[i];

    return result;
}

void BwtCorrector::showExecutionPlan()
{
    cout << "Run algorithm with..." << endl
         << "   Number of reads = " << numberOfReads_ << endl
         << "   Read length = " << readLength_ << endl
         << "   Genome length = " << genomeLength_ << endl
         << "   Error rate = " << errorRate_ << endl
         << "   Min witness length = " << minWitnessLength_ << endl;
}

void BwtCorrector::run()
{
    //run the backtracking algorithm and create the error objects...
    ErrorStore errors = findErrors();
    vector<ErrorInfo> errorLocations;

    int numErrors = errors.size();
    cout << "There were " << numErrors << " errors discovered..." << endl;

    //we have completed the error info, we do not need the bwt-positions... so just load them into a vector...
    for ( map<LetterNumber, ErrorInfo>::iterator it = errors.begin(); it != errors.end(); ++it )
        errorLocations.push_back(
            ErrorInfo(
                it->second.seqNum,
                it->second.readEnd,
                it->second.lastCycle,
                it->second.firstCycle,
                it->second.corrector
            )
        );

    //the ErrorInfo objects have their seqNum fields set to the position of the target read as if it was
    //sorted alphabetically, so map it back to the original ordering (as it appeared in the FAST* file )
    //using the -end-pos file and the SetReadNumbersToOriginal method...
    string endPosFileName = indexPrefix_ + "-end-pos";
    std::sort( errorLocations.begin(), errorLocations.end(), ErrorInfo::SortByRead );
    ErrorInfo::SetReadNumbersToOriginal( ( char * )endPosFileName.c_str(), errorLocations );

    //the ErrorInfo objects still have read numbers that refer to any of the original reads or their reverse complements
    //however we are only interested in the original reads, so we map errors discovered in the reverse complements back
    //to the originals with ConvertRCCorrectionsToOriginal...
    ErrorInfo::ConvertRCCorrectionsToOriginal( errorLocations, numberOfReads_, readLength_ );

    std::sort( errorLocations.begin(), errorLocations.end(), ErrorInfo::SortByRead );

    ErrorInfo::CorrectionsToCsv( outputFile_, errorLocations );

    cout << "Done!" << endl;
}
