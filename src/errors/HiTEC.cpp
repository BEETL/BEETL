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

#include "HiTEC.hh"

#include "BCRext.hh"
#include "BCRexternalBWT.hh"
#include "countWords/Range.hh"

#include "Timer.hh"
#include "config.h"

#include "parameters/BwtParameters.hh"
#include "libzoo/util/Logger.hh"
#include "shared/SeqReader.hh"
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#ifdef HAVE_POSIX_FADVISE
#endif

using namespace std;
using namespace BeetlBwtParameters;

HiTEC::HiTEC(
    const string &readsFile,
    double errorRate,
    double genomeLength
)
{
    setRequiredAttributes( readsFile, errorRate, genomeLength );
    setImpliedAttributes();
    setRecommendedAttributes();
}

void HiTEC::setWitnessLength( int witnessLength )
{
    wM_ = witnessLength;
    wm_ = witnessLength;
    for ( int i = 0; i < HITEC_MAX_ITERATIONS; i++ )
        witnessLengthSequence_[i] = witnessLength;
}

void HiTEC::setCorrectThreshold( int minOccurrences )
{
    minSupport_ = minOccurrences;
}

RangeStoreExternal *HiTEC::getPutativeErrors( int witnessLength )
{
    makeBwtAndLcp();
    RangeStoreExternal *r;
    r = new RangeStoreExternal( tempFilesPrefix_ + ".error.ranges.1", tempFilesPrefix_ + ".error.ranges.2" );
    int errorCount = 0;
    LetterCountEachPile cumulative = cumulativePileCounts();
    //can start at 1 because we are not interested in the $...
    for ( int i = 1; i < alphabetSize; i++ )
    {
        WitnessReader *wr = new WitnessReader(
            lcpFileName( i ),
            bwtFileName( i ),
            witnessLength,
            minSupport_,
            compressIntermediateBwts_
        );

        LetterCount lc;

        while ( wr->nextWitnessBlock( lc ) )
        {
            vector<int> correctLetters;
            vector<int> errorLetters;
            correctLetters.clear();
            errorLetters.clear();
            for ( int charNo = 0; charNo < alphabetSize; charNo++ )
                if ( lc.count_[charNo] >= minSupport_ )
                    correctLetters.push_back( charNo );
                else if ( lc.count_[charNo] >= 1 )
                    errorLetters.push_back( charNo );

            if ( correctLetters.size() == 1 )
            {
                for ( int errorCharIndex = 0; errorCharIndex < errorLetters.size(); errorCharIndex++ )
                {
                    int errorCharNo = errorLetters[errorCharIndex];

                    stringstream correctLetterString;
                    correctLetterString << alphabet[correctLetters[0]];

                    int positionInNextPile = wr->TotalCountSoFar().count_[errorCharNo] + cumulative[i - 1].count_[errorCharNo] - lc.count_[errorCharNo];
                    r->addRange(
                        Range(
                            correctLetterString.str(), //should be set to the accumulated string so far,
                            //but not vital for functionality - we are just going to hack this for the minute by putting the
                            //string of the correct character for this error-interval...
                            positionInNextPile,
                            lc.count_[errorCharNo],
                            false
                        )
                        , errorCharNo,
                        i,
                        "",
                        -1
                    );
                    errorCount += lc.count_[errorCharNo];
                }
            }
            else if ( correctLetters.size() > 1 )
            {
                for ( int errorCharIndex = 0; errorCharIndex < errorLetters.size(); errorCharIndex++ )
                    errorCount += lc.count_[errorLetters[errorCharIndex]];
            }
        }
        delete wr;
    }
    return r;
}

vector<HiTECErrorLocation> HiTEC::findErrorsInSortedReads( int witnessLength, RangeStoreExternal *errorIntervals )
{
    vector<HiTECErrorLocation> result;

    Timer  timer;
    vector <BwtReaderBase *> inBwt( alphabetSize );

    LetterCountEachPile countsPerPile, countsCumulative;

    int numCycles( readLength_ );
    int minOcc( numberOfReads_ );

    for ( int i( 0 ); i < alphabetSize; i++ )
    {
        string fileName = bwtFileName( i );
        if ( compressIntermediateBwts_ == true )
            inBwt[i] = new BwtReaderRunLengthIndex( fileName );
        else
            inBwt[i] = new BwtReaderASCII( fileName );
        inBwt[i]->readAndCount( countsPerPile[i] );
    }

    countsCumulative = countsPerPile;

    for ( int i( 1 ); i < alphabetSize; i++ )
        countsCumulative[i] += countsCumulative[i - 1];

    countsCumulative.print();

    string thisWord;

    LetterCount countsSoFar;
    LetterNumber currentPos;
    LetterNumber numRanges, numSingletonRanges;

    errorIntervals->clear();
    for ( int c( 0 ); c < numCycles; ++c )
    {
#ifdef PROPAGATE_PREFIX
        thisWord.resize( c + 3 );
#endif
        numRanges = 0;
        numSingletonRanges = 0;
        errorIntervals->swap();

        for ( int i( 1 ); i < alphabetSize; ++i )
        {
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
                errorIntervals->setPortion( i, j );

                OneBwtBackTracker backTracker(
                    inBwt[i],
                    currentPos,
                    ( *errorIntervals ),
                    countsSoFar,
                    minOcc,
                    numCycles,
                    subset_,
                    c + 2,
                    doesPropagateBkptToSeqNumInSet_,
                    true
                );

                ErrorIntervalHandler intervalHandler( witnessLength, c, result );
                backTracker( i, thisWord, intervalHandler );

                numRanges += backTracker.numRanges_;
                numSingletonRanges += backTracker.numSingletonRanges_;

            } // ~for j
            //cerr << "Done i " << i <<endl;
        }     // ~for i
        errorIntervals->clear();
        //    return 0; // %%%
        if ( numRanges == 0 ) break;

    } // ~for c

    //remember to close the files!
    for ( int i = 0; i < alphabetSize; i++ )
        delete inBwt[i];
    return result;
}

void CorrectFastaInPlace( const string &fileName, vector<HiTECErrorLocation> &errors )
{
    fstream file ( fileName.c_str(), fstream::in | fstream::out );
    int numErrors = errors.size();
    int currentError = 0;
    int currentRead = -1;
    string line;
    while ( getline( file, line ) && currentError < numErrors )
    {
        if ( line[0] == '>' )
            currentRead++;
        int seekFromPos = 0;
        while ( currentRead == errors[currentError].readNum )
        {
            file.seekp( errors[currentError].positionInRead - seekFromPos, ios_base::cur );
            //file.readsome(read, errors[currentError].positionInRead);
            file.put( errors[currentError].corrector );
            //just read in rest of line...
            seekFromPos = errors[currentError].positionInRead + 1;
            currentError++;
        }
        getline( file, line );
    }
    file.close();
}

void HiTEC::run( void )
{
    int iteration = 0;
    int totalCorrections = 0;
    double proportionCorrected;
    ofstream correctionsFile;
    correctionsFile.open( "corrections.csv" );
    correctionsFile << "iteration read position correction" << endl; //the columns of the output csv file.
    do
    {
        cout << "HiTEC Iteration " << iteration << "..." << endl;
        RangeStoreExternal *errorIntervals = getPutativeErrors( witnessLengthSequence_[iteration] );

        vector<HiTECErrorLocation> errorLocations = findErrorsInSortedReads(
                    witnessLengthSequence_[iteration],
                    errorIntervals
                );

        HiTECErrorLocation::SetReadNumbersToOriginal( ( char * )"hitec.temp-end-pos", errorLocations );

        delete errorIntervals;

        int correctedThisIteration = errorLocations.size();
        totalCorrections += correctedThisIteration;

        //interpret the reverse complement errors as errors in original reads...
        HiTECErrorLocation::ConvertRCCorrectionsToOriginal( errorLocations, numberOfReads_, readLength_ );

        //sort errors by their original read number and position...
        std::sort( errorLocations.begin(), errorLocations.end(), HiTECErrorLocation::ErrorLocationSorter );

        cout << correctedThisIteration << " Errors Found (showing first 20):" << endl;
        for ( int errNo = 0; errNo < min( correctedThisIteration, 20 ); errNo++ )
            errorLocations[errNo].print();

        for ( int errNo = 0; errNo < correctedThisIteration; errNo++ )
        {
            correctionsFile
                    << iteration << " "
                    << errorLocations[errNo].readNum << " "
                    << errorLocations[errNo].positionInRead << " "
                    << errorLocations[errNo].corrector << endl;
        }

        cout << endl;

        CorrectFastaInPlace( readsFile_, errorLocations );
        //delete errorLocations;
        proportionCorrected = ( double )correctedThisIteration / ( double )( readLength_ * numberOfReads_ );
        iteration++;
        //clear out files...
        removeUsedFiles();

        //clear out errors...

    }
    while ( iteration < HITEC_MAX_ITERATIONS && proportionCorrected >= proportionCorrectedThreshold_ );

    correctionsFile.close();
    cout << "Total corrections made: " << totalCorrections << endl;
}

void HiTEC::showExecutionPlan()
{
    cout << "Launching HiTEC with following parameters:" << endl;
    cout << "   Error Rate: " << errorRate_ << endl;
    cout << "   Genome Length: " << genomeLength_ << endl;
    cout << "   Number Of Reads: " << numberOfReads_ << endl;
    cout << "   Read Length: " << readLength_ << endl;
    cout << "   wm = " << wm_ << endl;
    cout << "   wM = " << wM_ << endl;
    cout << "   Minimum support for witness to be deemed correct = " << minSupport_ << endl;
    cout << "   Witness lengths to be used at successive iterations of algorithm:" << endl;
    cout << "   Number Of HiTEC Iterations: " << HITEC_MAX_ITERATIONS << endl;
    cout << "       ";
    for ( int i = 0; i < HITEC_MAX_ITERATIONS; i++ )
        cout << witnessLengthSequence_[i] << " ";
    cout << endl;
}

void HiTEC::makeBwtAndLcp( void )
{
    BwtParameters *params;

    int bcrMode = 0 ;

    params = new BwtParameters;
    params->commitDefaultValues();

    ( *params )["intermediate format"] = INTERMEDIATE_FORMAT_ASCII;
    ( *params )["output format"] = compressIntermediateBwts_ ? OUTPUT_FORMAT_RLE : OUTPUT_FORMAT_ASCII;
    ( *params )["generate LCP"] = 1;
    ( *params )["algorithm"] = "bcr";
    ( *params )["output filename"] = tempFilesPrefix_;
    ( *params )["generate endPosFile"] = 1;
    ( *params )["add reverse complement"] = 1;

    char *readsFileChars = new char[readsFile_.size() + 1];
    strcpy( readsFileChars, readsFile_.c_str() );

    char *tempPrefixChars = new char[tempFilesPrefix_.size() + 1];
    strcpy( tempPrefixChars, tempFilesPrefix_.c_str() );

    BCRexternalBWT bwt(
        readsFileChars,
        tempPrefixChars,
        bcrMode,
        compressionASCII,
        params
    );

    delete params;
}

string HiTEC::lcpFileName( int letterNumber )
{
    stringstream fileNameSS;
    fileNameSS << tempFilesPrefix_ << "-L0" << letterNumber;
    return fileNameSS.str().c_str();
}

string HiTEC::bwtFileName( int letterNumber )
{
    stringstream fileNameSS;
    fileNameSS << tempFilesPrefix_ << "-B0" << letterNumber;
    return fileNameSS.str().c_str();
}

LetterCountEachPile HiTEC::cumulativePileCounts()
{
    vector <BwtReaderBase *> inBwt( alphabetSize );
    LetterCountEachPile countsPerPile, countsCumulative;

    for ( int i( 0 ); i < alphabetSize; i++ )
    {
        string fileName = bwtFileName( i );
        if ( compressIntermediateBwts_ == true )
            inBwt[i] = new BwtReaderRunLengthIndex( fileName );
        else
            inBwt[i] = new BwtReaderASCII( fileName );
        inBwt[i]->readAndCount( countsPerPile[i] );
    }

    countsCumulative = countsPerPile;

    for ( int i( 1 ); i < alphabetSize; i++ )
        countsCumulative[i] += countsCumulative[i - 1];

    return countsCumulative;
}

void HiTEC::setRequiredAttributes(
    const string &readsFile,
    double errorRate,
    double genomeLength
)
{
    cout << "Initialising HiTEC variables..." << endl;
    genomeLength_ = genomeLength;
    errorRate_ = errorRate;
    readsFile_ = readsFile;
    doesPropagateBkptToSeqNumInSet_ = true;
    subset_ = "";
    tempFilesPrefix_ = "hitec.temp";
    workingReadsFile_ = "reads.hitec.temp";
    compressIntermediateBwts_ = true;
    proportionCorrectedThreshold_ = 0.0001;
}

void HiTEC::setImpliedAttributes()
{
    char *readsFileChars = new char[readsFile_.size() + 1];
    strcpy( readsFileChars, readsFile_.c_str() );
    FILE *tempReads = fopen( readsFileChars, "r" );

    SeqReaderFasta fasta( tempReads );

    readLength_ = fasta.length();

    numberOfReads_ = 0;
    while ( !fasta.allRead() )
    {
        fasta.readNext();
        numberOfReads_++;
    }

    fclose( tempReads );

    cout << "finished setting the required attributes..." << endl;
}

void HiTEC::setRecommendedAttributes()
{
    HiTECStats stats(
        errorRate_,
        genomeLength_,
        numberOfReads_,
        readLength_
    );
    cout << "setting hitec stats with error rate = " << errorRate_ << ", genome length: " << genomeLength_ << ", num reads: " << numberOfReads_ << ", " << readLength_ << endl;


    wm_ = stats.Calculate_wm();
    wM_ = stats.Calculate_wM();
    minSupport_ = stats.CalculateSupport( wM_ );

    witnessLengthSequence_[0] = wm_ + 1;
    witnessLengthSequence_[1] = wM_ + 1;
    witnessLengthSequence_[2] = wM_ + 1;
    witnessLengthSequence_[3] = wm_;
    witnessLengthSequence_[4] = wM_;
    witnessLengthSequence_[5] = wM_;
    witnessLengthSequence_[6] = wm_ - 1;
    witnessLengthSequence_[7] = wM_ - 1;
    witnessLengthSequence_[8] = wM_ - 1;
}

void HiTEC::removeUsedFiles()
{
    //todo! less hacky!
    string usedFiles[15] =
    {
        "hitec.temp-B00",
        "hitec.temp-B03",
        "hitec.temp-B06",
        "hitec.temp-L01",
        "hitec.temp-L04",
        "hitec.temp-B01",
        "hitec.temp-B04",
        "hitec.temp-end-pos",
        "hitec.temp-L02",
        "hitec.temp-L05",
        "hitec.temp-B02",
        "hitec.temp-B05",
        "hitec.temp-L00",
        "hitec.temp-L03",
        "hitec.temp-L06"
    };

    for ( int i = 0; i < 15; i++ )
    {
        if ( remove( usedFiles[i].c_str() ) != 0 )
        {
            cerr << "Error removing temporary file " << usedFiles[i] << endl;
            exit( 1 );
        }
    }
}
