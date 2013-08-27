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

#include "CountWords.hh"

#include "Timer.hh"
#include "config.h"
#include "libzoo/util/Logger.hh"

using namespace std;


CountWords::CountWords( bool inputACompressed,
                        bool inputBCompressed, char whichHandler,
                        int paramN, int paramK, const vector<string> &setA,
                        const vector<string> &setB, const vector<string> &setC,
                        const string &ncbiTax, bool testDB, unsigned int minWordLen, string prefix, string subset
                        , const bool doesPropagateBkptToSeqNumInSetA
                        , const bool doesPropagateBkptToSeqNumInSetB
                        , const bool noComparisonSkip
                      ) :
    setA_( setA ), setB_( setB ), setC_( setC ), subset_( subset )
    , doesPropagateBkptToSeqNumInSetA_( doesPropagateBkptToSeqNumInSetA )
    , doesPropagateBkptToSeqNumInSetB_( doesPropagateBkptToSeqNumInSetB )
    , noComparisonSkip_( true ) // Always set to 1 for BEETL release 0.6, as this feature isn't tested // real value: noComparisonSkip
{

#ifdef XXX
    // get memory allocated
    setA_ = new char[setA.length() + 1];
    setB_ = new char[setB.length() + 1];

    // copt to char * in order to get valid c strings
    setA.copy( setA_, setA.length() );
    setB.copy( setB_, setB.length() );

    // append \0 to obtain a valid escaped c string
    setA_[setA.length()] = '\0';
    setB_[setB.length()] = '\0';
#endif
    /*Information for metagenomics search.
      only needed if countWords is used as a metagenome classifier
    */

    //ncbiInfo_ should be the information about which files in the database have which taxonomy
    ncbiInfo_ = ncbiTax;
    // set tool flags
    inputACompressed_ = inputACompressed;
    inputBCompressed_ = inputBCompressed;
    whichHandler_ = whichHandler;
    //flag if the database of the metagenomics information should also be tested
    //the testing of the database uses a lot of disk space so it should be handled carefully
    testDB_ = testDB;

    paramN_ = paramN;
    paramK_ = paramK;
    //only used in the metagenomics part. tests of taxonomy only after a certain suffix length is reached
    minWordLen_ = minWordLen;

    tmpPrefix_ = prefix;
}

void CountWords::run( void )
{
    Timer  timer;
    // input streams for BWT from previous iter - 1 per pile
    // Used to compute the counts to do backward search
    vector <BwtReaderBase *> inBwtA( alphabetSize );
    vector <BwtReaderBase *> inBwtB( alphabetSize );
    //  RangeStoreRAM rA,rB;
    RangeStoreExternal rA( tmpPrefix_ + "count_A_1", tmpPrefix_ + "count_A_2" ), rB( tmpPrefix_ + "count_B_1", tmpPrefix_ + "count_B_2" );
    LetterCountEachPile countsPerPileA, countsCumulativeA;
    LetterCountEachPile countsPerPileB, countsCumulativeB;
    bool referenceMode( false ); // if false, use splice junction mode
    bool metagenomeMode( false );

    //Christina:
    //now three different possibilities for count words, so a little parsing is needed here
    referenceMode = ( whichHandler_ == 'r' ) ? true : false;

    metagenomeMode = ( whichHandler_ == 'm' ) ? true : false;

    if ( referenceMode && testDB_ )
    {
        cerr << "WARNING Database test in reference Mode not possible" << endl
             << "ABORTING Dastabase test" << endl;
        testDB_ = false;
    }
    int numCycles( paramK_ );
    int minOcc( paramN_ );

    if ( metagenomeMode )
    {
        loadFileNumToTaxIds( ncbiInfo_ );
    }
    vector<FILE *> mergeCSet;
    mergeCSet.resize( alphabetSize );
    #pragma omp parallel for
    for ( int i = 0; i < alphabetSize; i++ )
    {
        if ( metagenomeMode )
        {
            //open the files containing the file numbers, in those files there should be a fileNumber for each BWTpostion
            //indicating from which the this suffix came
            FILE *mergC = fopen( setC_[i].c_str(), "r" );
            if ( mergC == NULL )
            {
                #pragma omp critical
                {
                    cout << "Could not open File \"" << setC_[i] << "\"" << endl;
                }
            }
            assert( mergC != NULL );
            mergeCSet[i] = fopen( setC_[i].c_str(), "r" );
        }
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
        {
            #pragma omp critical
            {
                Logger::out() << "i=" << i << ", alphabet[i]=" << alphabet[i] << endl;
                Logger::out() << "setA[i]=" << setA_[i] << ", setB[i]=" << setB_[i] << endl;
            }
        }
        if ( inputACompressed_ == true )
            inBwtA[i] = new BwtReaderRunLengthIndex( setA_[i] );
        else
            inBwtA[i] = new BwtReaderASCII( setA_[i] );
        inBwtA[i]->readAndCount( countsPerPileA[i] );

        if ( inputBCompressed_ == true )
            inBwtB[i] = new BwtReaderRunLengthIndex( setB_[i] );
        else
            inBwtB[i] = new BwtReaderASCII( setB_[i] );
        //    inBwtB[i]= new BwtReaderASCII(args[i+alphabetSize+nextArg]);
        inBwtB[i]->readAndCount( countsPerPileB[i] );

        Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
        {
            #pragma omp critical
            {
                countsPerPileA[i].print();
                countsPerPileB[i].print();
            }
        }
    }

    countsCumulativeA = countsPerPileA;
    countsCumulativeB = countsPerPileB;
    for ( int i( 1 ); i < alphabetSize; i++ )
    {
        countsCumulativeA[i] += countsCumulativeA[i - 1];
        countsCumulativeB[i] += countsCumulativeB[i - 1];
    }

    Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
    {
        countsCumulativeA.print();
        countsCumulativeB.print();
    }
    string thisWord;
    const int dontKnowIndex( whichPile[( int )dontKnowChar] );

    // sort out first iter
    for ( int i( 1 ); i < alphabetSize; ++i )
    {
        for ( int j( 1 ); j < alphabetSize; ++j )
        {
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                Logger::out() << i << " " << j << " " << matchFlag << endl;
                Logger::out() << ( countsCumulativeA[i - 1].count_[j]
                                   | ( matchFlag * ( countsPerPileB[i].count_[j] != 0 ) ) )
                              << " " << ( countsCumulativeB[i - 1].count_[j]
                                          | ( matchFlag * ( countsPerPileA[i].count_[j] != 0 ) ) ) << endl;
            }
#ifdef PROPAGATE_PREFIX
            thisWord.clear();
            thisWord += alphabet[j];
            thisWord += alphabet[i];
#else
            /*
                        if (!subset.empty() && subset[subset.size()-1] != alphabet[i])
                        {
                            Logger::out() << "  skipped" << endl;
                            continue;
                        }
            */
#endif

            if ( ( i != dontKnowIndex ) && ( j != dontKnowIndex ) )
            {
                // get rid of any ranges with N in them
                if ( countsPerPileA[i].count_[j] != 0 )
                    rA.addRange( j, i, thisWord,
                                 ( countsCumulativeA[i - 1].count_[j]
                                   | ( matchFlag * ( LetterNumber )( countsPerPileB[i].count_[j] != 0 ) ) ),
                                 countsPerPileA[i].count_[j], false, subset_, 1 );
                if ( countsPerPileB[i].count_[j] != 0 )
                    rB.addRange( j, i, thisWord,
                                 ( countsCumulativeB[i - 1].count_[j]
                                   | ( matchFlag * ( countsPerPileA[i].count_[j] != 0 ) ) ),
                                 countsPerPileB[i].count_[j], false, subset_, 1 );
            } // ~if
        } // ~for j
    } // ~for i


    LetterCount countsSoFarA, countsSoFarB;
    LetterNumber currentPosA, currentPosB;
    LetterNumber numRanges, numSingletonRanges;
    //  Range thisRangeA,thisRangeB;
    rA.clear();
    rB.clear();
    for ( int c( 0 ); c < numCycles; ++c )
    {
        cerr << "Starting iteration " << c << ", time now: " << timer.timeNow();
        cerr << "Starting iteration " << c << ", usage: " << timer << endl;

#ifdef PROPAGATE_PREFIX
        thisWord.resize( c + 3 );
#endif

        numRanges = 0;
        numSingletonRanges = 0;

        rA.swap();
        rB.swap();

        //    pTemp=pNext;pNext=pThis;pThis=pTemp;
        for ( int i( 1 ); i < alphabetSize; ++i )
        {
            inBwtA[i]->rewindFile();
            inBwtB[i]->rewindFile();
            currentPosA = 0;
            currentPosB = 0;
            countsSoFarA.clear();
            countsSoFarA += countsCumulativeA[i - 1];
            countsSoFarB.clear();
            countsSoFarB += countsCumulativeB[i - 1];

#ifdef PROB_NOT_NEEDED
            BwtReaderRunLengthIndex *pRun;

            pRun = dynamic_cast<BwtReaderRunLengthIndex *>( inBwtA[i] );
            if ( pRun != NULL )
                pRun->initIndex( countsSoFarA );
            else
                inBwtA[i]->rewindFile();
            //cout <<"fred" << endl;
            pRun = dynamic_cast<BwtReaderRunLengthIndex *>( inBwtB[i] );
            if ( pRun != NULL )
                pRun->initIndex( countsSoFarB );
            else
                inBwtB[i]->rewindFile();
#endif

            for ( int j( 1 ); j < alphabetSize; ++j )
            {
                Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "positions - A: " << currentPosA << " B: " << currentPosB << endl;
                rA.setPortion( i, j );
                rB.setPortion( i, j );

                BackTracker backTracker( inBwtA[i], inBwtB[i]
                                         , currentPosA, currentPosB
                                         , rA, rB, countsSoFarA, countsSoFarB
                                         , minOcc, numCycles, subset_, c + 2
                                         , doesPropagateBkptToSeqNumInSetA_
                                         , doesPropagateBkptToSeqNumInSetB_
                                         , noComparisonSkip_
                                       );
                if ( referenceMode == true )
                {
                    IntervalHandlerReference intervalHandler( minOcc );
                    backTracker( i, thisWord, intervalHandler );
                }
                else if ( metagenomeMode == true )
                {
                    IntervalHandlerMetagenome intervalHandler( minOcc, mergeCSet, fileNumToTaxIds_, testDB_, minWordLen_, paramK_ );
                    backTracker( i, thisWord, intervalHandler );
                }
                else
                {
                    IntervalHandlerSplice intervalHandler( minOcc );
                    backTracker( i, thisWord, intervalHandler );
                }
                numRanges += backTracker.numRanges_;
                numSingletonRanges += backTracker.numSingletonRanges_;

            } // ~for j
            //cerr << "Done i " << i <<endl;
        }     // ~for i
        rA.clear();
        rB.clear();
        cerr << "Finished cycle " << c << ": ranges=" << numRanges
             << " singletons=" << numSingletonRanges
             << " usage: " << timer << endl;

        //    return 0; // %%%
        if ( numRanges == 0 ) break;
    } // ~for c
} // countWords::run()

/*
  @taxIdNames , file containing the taxonomic information about the database
  loads the taxonomic information so it is known which file number corresponds to which taxa
  fills up the vector fileNumToTaxIds;
  This should be called before the comparison starts
*/

void CountWords::loadFileNumToTaxIds( const string &taxIdNames )
{
    ifstream taxas( taxIdNames.c_str(), ios::in );
    string line;
    //each line should contain the fileNumber followed by the taxIds split up with one space character
    while ( taxas.good() )
    {
        vector<int> taxIDs;
        getline( taxas, line );
        if ( line.compare( "" ) != 0 )
        {
            unsigned int fileNum = atoi( line.substr( 0, line.find( " " ) ).c_str() );
            line = line.substr( line.find( " " ) + 1, line.length() );
            while ( line.find( " " ) != string::npos )
            {

                taxIDs.push_back( atoi( line.substr( 0, line.find( " " ) ).c_str() ) );
                line = line.substr( line.find( " " ) + 1, line.length() );
            }
            taxIDs.push_back( atoi( line.c_str() ) );
            //test if all TaxIds were found
            if ( taxIDs.size() < taxLevelSize )
            {
                cerr << "Tax Ids have not enough taxonomic Information. Only  " << taxIDs.size() << " could be found " << endl
                     << "Will add unknown taxa until size is right" << endl;
                for ( unsigned int i( taxIDs.size() - 1 ) ; i < taxLevelSize; i++ )
                    taxIDs.push_back( 0 );
            }
            if ( taxIDs.size() > taxLevelSize )
                cerr << "Tax Ids have to much taxonomic information. "
                     << "Please note, that the taxonomic information about one file should be shown as: " << endl
                     << "FileNumber Superkingdom Phylum Order Family Genus Species Strain " << endl;
            fileNumToTaxIds_.push_back( taxIDs );
            unsigned int test = fileNum + 1;
            if ( test != fileNumToTaxIds_.size() )
                cout << "Wrong filenumber " << fileNum << " " << fileNumToTaxIds_.size() << endl;
        }
    }
    //cout << " fineNumToTaxIds " << fileNumToTaxIds.size() <<endl;
}

