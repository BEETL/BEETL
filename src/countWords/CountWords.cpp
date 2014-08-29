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

#include "CountWords.hh"

#include "Timer.hh"
#include "Tools.hh"
#include "config.h"
#include "libzoo/util/Logger.hh"

#include <fcntl.h>
#include <sstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _OPENMP
# include <omp.h>
#endif //ifdef _OPENMP

using namespace std;

//#define DEBUG__FORCE_WRITE_COUNT
//#define DEBUG__SKIP_ITERATION0
//#define FIRST_CYCLE 14


CountWords::CountWords( bool inputACompressed,
                        bool inputBCompressed, char whichHandler,
                        int paramN, int paramK, const vector<string> &setA,
                        const vector<string> &setB, const vector<string> &setC,
                        const string &ncbiTax, bool testDB, unsigned int minWordLen, string subset
                        , const CompareParameters *compareParams
                      )
    : numCycles_( paramK )
    , minOcc_( paramN )
    , setA_( setA )
    , setB_( setB )
    , setC_( setC )
    , subset_( subset )
    , compareParams_( compareParams )
    , doPauseBetweenCycles_( compareParams_ ? ( *compareParams_ )["pause between cycles"] : false )
    , doesPropagateBkptToSeqNumInSetA_( compareParams_ ? ( *compareParams_ )["generate seq num A"] : false )
    , doesPropagateBkptToSeqNumInSetB_( compareParams_ ? ( *compareParams_ )["generate seq num B"] : false )
    , noComparisonSkip_( compareParams_ ? ( *compareParams_ )["no comparison skip"] : false )
    , bwtInRam_( compareParams_ ? ( *compareParams_ )["BWT in RAM"] : false )
    , propagateSequence_( compareParams_ ? ( *compareParams_ )["propagate sequence"] : false )
    , outputDirectory_( compareParams_ ? ( *compareParams_ )["output directory"] : string( "BeetlCompareOutput" ) )
    , inBwtA_( alphabetSize )
    , inBwtB_( alphabetSize )
    , fsizeRatio_( 0 )
{

    if ( compareParams_ == NULL )
    {
        // Legacy mode, used only by the OldBeetl executable: We create a CompareParameters structure with default values
        compareParams_ = new CompareParameters;

        // Three different modes, so a little parsing is needed (tumour-normal is not a legacy mode)
        switch ( whichHandler )
        {
            case 's':
                mode_ = BeetlCompareParameters::MODE_SPLICE;
                break;
            case 'r':
                mode_ = BeetlCompareParameters::MODE_REFERENCE;
                break;
            case 'm':
                mode_ = BeetlCompareParameters::MODE_METAGENOMICS;
                break;
            default:
                assert( false && "unexpected compare mode" );
        }
    }
    else
        mode_ = static_cast<enum BeetlCompareParameters::Mode>( compareParams_->getValue( "mode" ) );

    // set tool flags
    inputACompressed_ = inputACompressed;
    inputBCompressed_ = inputBCompressed;

    // Check that output directory is empty. Create it if necessary.
    mkdir( outputDirectory_.c_str(), 0750 );
    if ( isDirectoryEmpty( outputDirectory_ ) == 0 )
    {
        Logger::error() << "Error: Output directory not empty" << endl;
        assert( false && "output directory not empty" );
    }

    /*
      Information for metagenomics search.
      Only needed if countWords is used as a metagenome classifier
    */
    //ncbiInfo_ should be the information about which files in the database have which taxonomy
    ncbiInfo_ = ncbiTax;
    //flag if the database of the metagenomics information should also be tested
    //the testing of the database uses a lot of disk space so it should be handled carefully
    testDB_ = testDB;
    //only used in the metagenomics part. tests of taxonomy only after a certain suffix length is reached
    minWordLen_ = minWordLen;

    if ( testDB_ && mode_ != BeetlCompareParameters::MODE_METAGENOMICS )
    {
        cerr << "WARNING Database test is only available in metagenome Mode" << endl
             << "ABORTING Database test" << endl;
        testDB_ = false;
    }
}

void CountWords::run( void )
{
    Timer  timer;
    LetterCountEachPile countsPerPileA;
    LetterCountEachPile countsPerPileB;
    const bool useThreadsForSubsets = true;
    int cyclesToSkipComparisonFor = -1;
    int previousComparisonDeactivationLength = 0;

#ifdef _OPENMP
    // Use nested openmp parallelisation
    omp_set_nested( 1 );
#endif

    // Metagenomics-specific stuff
    if ( mode_ == BeetlCompareParameters::MODE_METAGENOMICS )
        initialiseMetagomeMode();

    vector<RangeStoreExternal *> rangeStoresA;
    vector<RangeStoreExternal *> rangeStoresB;
    vector<int> pileToThread;
    if ( !useThreadsForSubsets )
    {
        rangeStoresA.push_back( new RangeStoreExternal( propagateSequence_, "Intervals_setA" ) );
        rangeStoresB.push_back( new RangeStoreExternal( propagateSequence_, "Intervals_setB" ) );
        pileToThread.resize( alphabetSize, 0 ); // All piles point to thread 0
    }
    else
    {
        for ( int i = 0; i < 4; ++i )
        {
            ostringstream oss;
            oss << "Intervals_subset" << i;
            rangeStoresA.push_back( new RangeStoreExternal( propagateSequence_, oss.str() + "_setA" ) );
            rangeStoresB.push_back( new RangeStoreExternal( propagateSequence_, oss.str() + "_setB" ) );
        }
        pileToThread.resize( alphabetSize, 0 );
        pileToThread[whichPile['A']] = 0; // pile 'A' processed by thread 0
        pileToThread[whichPile['C']] = 1; // pile 'C' processed by thread 1
        pileToThread[whichPile['G']] = 2; // pile 'G' processed by thread 2
        pileToThread[whichPile['T']] = 3; // pile 'T' processed by thread 3
        pileToThread[whichPile['N']] = 1; // pile 'N' processed by thread 1
#ifdef USE_EXTRA_CHARACTER_Z
        pileToThread[whichPile['Z']] = 2; // pile 'Z' processed by thread 2
#endif
    }

    //    cout << "OMP threadCount: " << omp_get_num_threads() << endl;
    //omp_set_num_threads( 2 * alphabetSize );

    // We reproduce the same nested parallel distribution here than during the real processing loops, in order to load the BWT data in the most appropriate NUMA nodes
    const int subsetThreadCount = rangeStoresA.size();
#if defined(_OPENMP)
    omp_set_num_threads( subsetThreadCount * alphabetSize );
#endif

    #pragma omp parallel for
    for ( int threadNum = 0; threadNum < subsetThreadCount * alphabetSize; ++threadNum )
    {
        int subsetThreadNum = threadNum / alphabetSize; // for ( int subsetThreadNum = 0; subsetThreadNum < subsetThreadCount; ++subsetThreadNum )
        int i = threadNum % alphabetSize; // for ( int i = 0; i < alphabetSize; ++i )

        int datasetAorB;
        switch ( subsetThreadNum )
        {
                // set A goes to same node as subsets 0&1: node 0
                // We distribute the read&Count as follows: ACGT=>subset 0's cpus, others=>subset 1's cpus
            case 0:
                if ( alphabet[i] == 'A' || alphabet[i] == 'C' || alphabet[i] == 'G' || alphabet[i] == 'T' )
                    datasetAorB = 0;
                else
                    continue;
                break;
            case 1:
                if ( alphabet[i] == 'A' || alphabet[i] == 'C' || alphabet[i] == 'G' || alphabet[i] == 'T' )
                    continue;
                else
                    datasetAorB = 0;
                break;

                // set B goes to same node as subsets 2&3: node 1
                // We distribute the read&Count as follows: ACGT=>subset 2's cpus, others=>subset 3's cpus
            case 2:
                if ( alphabet[i] == 'A' || alphabet[i] == 'C' || alphabet[i] == 'G' || alphabet[i] == 'T' )
                    datasetAorB = 1;
                else
                    continue;
                break;
            case 3:
                if ( alphabet[i] == 'A' || alphabet[i] == 'C' || alphabet[i] == 'G' || alphabet[i] == 'T' )
                    continue;
                else
                    datasetAorB = 1;
                break;

            default:
                continue;
        }
        /*
                omp_set_num_threads( alphabetSize );
                #pragma omp parallel for
                for ( int i = 0; i < alphabetSize; ++i )
        */
        {
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                int pid = -1;
                int numThreads = -1;
                int processor = -1;
                readProcSelfStat( pid, numThreads, processor );
                #pragma omp critical (IO)
                cerr << "Reading BWT: subsetThreadNum=" << subsetThreadNum << " i=" << i << " pid=" << pid << " numThreads=" << numThreads << " processor=" << processor << endl;
            }
            if ( !isDistributedProcessResponsibleForPile( i ) )
                continue;

            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                #pragma omp critical (IO)
                {
                    Logger::out() << "i=" << i << ", alphabet[i]=" << alphabet[i] << ": ";
                    if ( datasetAorB == 0 )
                        Logger::out() << "setA[i]=" << setA_[i] << endl;
                    else
                        Logger::out() << "setB[i]=" << setB_[i] << endl;
                }
            }

            if ( datasetAorB == 0 )
            {
                inBwtA_[i] = instantiateBwtPileReader( setA_[i], compareParams_->getStringValue( "use shm" ), bwtInRam_ );
#ifndef DEBUG__SKIP_ITERATION0
                Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
                {
                    int pid = -1;
                    int numThreads = -1;
                    int processor = -1;
                    readProcSelfStat( pid, numThreads, processor );
                    #pragma omp critical (IO)
                    Logger::out() << "ReadAndCount set=A i=" << i << " pid=" << pid << " numThreads=" << numThreads << " processor=" << processor << endl;
                }
                inBwtA_[i]->readAndCount( countsPerPileA[i] );
#endif //ifndef DEBUG__SKIP_ITERATION0
            }
            else
            {
                inBwtB_[i] = instantiateBwtPileReader( setB_[i], compareParams_->getStringValue( "use shm" ), bwtInRam_ );
#ifndef DEBUG__SKIP_ITERATION0
                Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
                {
                    int pid = -1;
                    int numThreads = -1;
                    int processor = -1;
                    readProcSelfStat( pid, numThreads, processor );
                    #pragma omp critical (IO)
                    Logger::out() << "ReadAndCount set=B i=" << i << " pid=" << pid << " numThreads=" << numThreads << " processor=" << processor << endl;
                }
                inBwtB_[i]->readAndCount( countsPerPileB[i] );
#endif //ifndef DEBUG__SKIP_ITERATION0
            }

            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                #pragma omp critical (IO)
                {
                    if ( datasetAorB == 0 )
                        Logger::out() << "countsPerPileA[i=" << i << "]: " << countsPerPileA[i] << endl;
                    else
                        Logger::out() << "countsPerPileB[i=" << i << "]: " << countsPerPileB[i] << endl;
                }
            }

#ifndef DEBUG__SKIP_ITERATION0
            // Share counts via files if distributed processing
#ifndef DEBUG__FORCE_WRITE_COUNT
            if ( ( *compareParams_ )["4-way distributed"].isSet() )
#endif
            {
                if ( datasetAorB == 0 )
                {
                    stringstream ss;
                    ss << "counts.A." << i;
                    ofstream ofs( ss.str() );
                    ofs << countsPerPileA[i] << endl;
                }
                else
                {
                    stringstream ss;
                    ss << "counts.B." << i;
                    ofstream ofs( ss.str() );
                    ofs << countsPerPileB[i] << endl;
                }
            }
#endif //ifndef DEBUG__SKIP_ITERATION0
        }
    }

    // Free bwt-B00, especially useful when keeping BWT in RAM
    delete inBwtA_[0];
    inBwtA_[0] = NULL;
    delete inBwtB_[0];
    inBwtB_[0] = NULL;

    cerr << "Finished initialisation cycle, usage: " << timer << endl;

    if ( doPauseBetweenCycles_ )
        pauseBetweenCycles();

    cerr << "Starting iteration " << 0 << ", time now: " << timer.timeNow();
    cerr << "Starting iteration " << 0 << ", usage: " << timer << endl;

    // Gather counts via files if distributed processing
#ifndef DEBUG__SKIP_ITERATION0
    if ( ( *compareParams_ )["4-way distributed"].isSet() )
#endif
    {
        for ( int i = 0; i < alphabetSize; i++ )
        {
            {
                stringstream ss;
                ss << "counts.A." << i;
                ifstream ifs( ss.str() );
                ifs >> countsPerPileA[i];
            }
            {
                stringstream ss;
                ss << "counts.B." << i;
                ifstream ifs( ss.str() );
                ifs >> countsPerPileB[i];
            }

            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                #pragma omp critical (IO)
                {
                    Logger::out() << "countsPerPileA[i=" << i << "]: " << countsPerPileA[i] << endl;
                    Logger::out() << "countsPerPileB[i=" << i << "]: " << countsPerPileB[i] << endl;
                }
            }
        }
    }

    // Calculating fsizeRatio
    {
        LetterNumber pile0LengthA = 0;
        LetterNumber pile0LengthB = 0;
        for ( int i( 0 ); i < alphabetSize; i++ )
        {
            pile0LengthA += countsPerPileA[0].count_[i];
            pile0LengthB += countsPerPileB[0].count_[i];
        }
        fsizeRatio_ =  double( pile0LengthA ) / double( pile0LengthB );
    }

    countsCumulativeA_ = countsPerPileA;
    countsCumulativeB_ = countsPerPileB;
    for ( int i( 1 ); i < alphabetSize; i++ )
    {
        countsCumulativeA_[i] += countsCumulativeA_[i - 1];
        countsCumulativeB_[i] += countsCumulativeB_[i - 1];
    }

    Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
    {
        Logger::out() << "cumulative counts A: " << endl;
        countsCumulativeA_.print();
        Logger::out() << "cumulative counts B: " << endl;
        countsCumulativeB_.print();
    }

    /*
        #pragma omp parallel for
        for ( int subsetNum = 0; subsetNum < 4; subsetNum++ )
        {
            string subsetValues[] = { "A", "C", "G", "T" };
            string threadSubset = subsetValues[subsetNum];
    */

#ifndef DEBUG__SKIP_ITERATION0
    const int dontKnowIndex( whichPile[( int )dontKnowChar] );

    // sort out first iter
    string currentWord = "xx";
    for ( int i( 1 ); i < alphabetSize; ++i )
    {
        const int subsetThreadNum = pileToThread[i];
        if ( propagateSequence_ )
            currentWord[1] = alphabet[i];
        for ( int j( 1 ); j < alphabetSize; ++j )
        {
            Logger_if( LOG_FOR_DEBUGGING )
            {
                #pragma omp critical (IO)
                {
                    Logger::out() << i << " " << j << " " << matchFlag << endl;
                    Logger::out() << ( countsCumulativeA_[i - 1].count_[j]
                    | ( matchFlag * ( countsPerPileB[i].count_[j] != 0 ) ) )
                    << " " << ( countsCumulativeB_[i - 1].count_[j]
                    | ( matchFlag * ( countsPerPileA[i].count_[j] != 0 ) ) ) << endl;
                }
            }
            if ( propagateSequence_ )
                currentWord[0] = alphabet[j];
            /*

                                    if (!subset.empty() && subset[subset.size()-1] != alphabet[i])
                                    {
                                        Logger::out() << "  skipped" << endl;
                                        continue;
                                    }
            */

            if ( ( i != dontKnowIndex ) && ( j != dontKnowIndex ) ) // don't process ranges with N in them
            {
                if ( countsPerPileA[i].count_[j] != 0 )
                    rangeStoresA[subsetThreadNum]->addRange( Range(  currentWord,
                            ( countsCumulativeA_[i - 1].count_[j]
                              | ( matchFlag * ( LetterNumber )( countsPerPileB[i].count_[j] != 0 ) ) ),
                            countsPerPileA[i].count_[j], false )
                            , j, i, subset_, 1 );
                if ( countsPerPileB[i].count_[j] != 0 )
                    rangeStoresB[subsetThreadNum]->addRange( Range( currentWord,
                            ( countsCumulativeB_[i - 1].count_[j]
                              | ( matchFlag * ( countsPerPileA[i].count_[j] != 0 ) ) ),
                            countsPerPileB[i].count_[j], false )
                            , j, i, subset_, 1 );
            } // ~if
        } // ~for j
    } // ~for i
#endif //ifndef DEBUG__SKIP_ITERATION0

    // Get ready for next cycle (flush current files and delete next cycle's output files)
    for ( auto rangeStore : rangeStoresA )
        rangeStore->clear();
    for ( auto rangeStore : rangeStoresB )
        rangeStore->clear();

    cerr << "Finished cycle 0, usage: " << timer << endl;

    if ( doPauseBetweenCycles_ )
        pauseBetweenCycles();

#ifdef DEBUG__SKIP_ITERATION0
    //    omp_set_num_threads( 1 ); // for debugging
#endif
#ifdef FIRST_CYCLE
    const int firstCycle = FIRST_CYCLE;
#else
    const int firstCycle = 1;
#endif
    for ( int cycle = firstCycle; cycle <= numCycles_; ++cycle )
    {
        cerr << "Starting iteration " << cycle << ", time now: " << timer.timeNow();
        cerr << "Starting iteration " << cycle << ", usage: " << timer << endl;

        numRanges_ = 0;
        numSingletonRanges_ = 0;
        numNotSkippedA_ = 0;
        numNotSkippedB_ = 0;
        numSkippedA_ = 0;
        numSkippedB_ = 0;

        /*
                #pragma omp parallel for
                for ( int subsetThreadNum = 0; subsetThreadNum < subsetThreadCount; ++subsetThreadNum )
                {

                    rangeStoresA[subsetThreadNum]->setCycleNum( cycle );
                    rangeStoresB[subsetThreadNum]->setCycleNum( cycle );

                    CountWords_parallelSubsetThread(
                        subsetThreadNum
                        , cycle
                        , *( rangeStoresA[subsetThreadNum] )
                        , *( rangeStoresB[subsetThreadNum] )
                    );
                }
        */
        // Sequential initialisation before the big parallel loop
        for ( int subsetThreadNum = 0; subsetThreadNum < subsetThreadCount; ++subsetThreadNum )
        {
            rangeStoresA[subsetThreadNum]->setCycleNum( cycle );
            rangeStoresB[subsetThreadNum]->setCycleNum( cycle );
        }

        #pragma omp parallel for
        for ( int threadNum = 0; threadNum < subsetThreadCount * alphabetSize; ++threadNum )
        {
            const int subsetThreadNum = threadNum / alphabetSize; // for ( int subsetThreadNum = 0; subsetThreadNum < subsetThreadCount; ++subsetThreadNum )
            //const int i = threadNum % alphabetSize; // for ( int i = 0; i < alphabetSize; ++i )

            CountWords_parallelSubsetThread(
                threadNum
                , cycle
                , *( rangeStoresA[subsetThreadNum] )
                , *( rangeStoresB[subsetThreadNum] )
            );
        }


        cerr << "Finished cycle " << cycle << ": ranges=" << numRanges_
             << " singletons=" << numSingletonRanges_
             << " usage: " << timer << endl;
        cerr
                << " skippedA=" << numSkippedA_
                << " skippedB=" << numSkippedB_
                << " notSkippedA=" << numNotSkippedA_
                << " notSkippedB=" << numNotSkippedB_
                << endl;

        if ( ( *compareParams_ )["no comparison skip"] == false )
        {
            if ( cyclesToSkipComparisonFor == -1 )
            {
                // if #intervals stopped increasing [by more than 2%], then it's time to start disabling comparisons
                if ( numNotSkippedA_ < numRanges_ * 1.02 )
                {
                    double ratio = numSkippedA_ / ( double )numNotSkippedA_;
                    if ( ratio > 0.42 )
                        cyclesToSkipComparisonFor = 0;
                    else if ( ratio > 0.32 )
                        cyclesToSkipComparisonFor = 1;
                    else if ( ratio > 0.25 ) //(0.3*3/(0.7*3+0.49*2+0.34))
                        cyclesToSkipComparisonFor = 2;
                    else if ( ratio > 0.2 )
                        cyclesToSkipComparisonFor = 4;
                    else if ( ratio > 0.1 )
                        cyclesToSkipComparisonFor = 6;
                    else if ( ratio > 0.01 )
                        cyclesToSkipComparisonFor = 10;
                    else if ( ratio > 0.001 )
                        cyclesToSkipComparisonFor = 20;
                    else
                        cyclesToSkipComparisonFor = 100;
                    previousComparisonDeactivationLength = cyclesToSkipComparisonFor;
                }
            }
            else
            {
                if ( noComparisonSkip_ )
                {
                    // it was deactivated => time to reactivate?
                    --cyclesToSkipComparisonFor;
                }
                else
                {
                    // it was active => time to deactivate?
                    double ratio = numSkippedA_ / ( double )numNotSkippedA_ / previousComparisonDeactivationLength;
                    if ( ratio > 0.42 )
                        cyclesToSkipComparisonFor = 0;
                    else if ( ratio > 0.32 )
                        cyclesToSkipComparisonFor = 1;
                    else if ( ratio > 0.25 ) //(0.3*3/(0.7*3+0.49*2+0.34))
                        cyclesToSkipComparisonFor = 2;
                    else if ( ratio > 0.2 )
                        cyclesToSkipComparisonFor = 4;
                    else if ( ratio > 0.1 )
                        cyclesToSkipComparisonFor = 6;
                    else if ( ratio > 0.01 )
                        cyclesToSkipComparisonFor = 10;
                    else if ( ratio > 0.001 )
                        cyclesToSkipComparisonFor = 20;
                    else
                        cyclesToSkipComparisonFor = 100;
                    previousComparisonDeactivationLength = cyclesToSkipComparisonFor;
                }
            }
            if ( cyclesToSkipComparisonFor == 0 && numNotSkippedA_ <= numSingletonRanges_ )
            {
                // If the number of interval decreases, don't reactivate comparison skips
                Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Interval count decreasing. ";
                cyclesToSkipComparisonFor = 1;
            }
            if ( cyclesToSkipComparisonFor > 0 )
            {
                Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Disabling comparison skips for " << cyclesToSkipComparisonFor << " iterations" << endl;
                noComparisonSkip_ = true;
            }
            else
            {
                Logger_if( LOG_SHOW_IF_VERBOSE )
                {
                    if ( noComparisonSkip_ )
                        Logger::out() << "Re-activating comparison skips" << endl;
                    else
                        Logger::out() << "Keeping comparison skips active" << endl;
                }
                noComparisonSkip_ = false;
            }
        }

        if ( doPauseBetweenCycles_ )
            pauseBetweenCycles();

        for ( auto rangeStore : rangeStoresA )
            rangeStore->clear( false );
        for ( auto rangeStore : rangeStoresB )
            rangeStore->clear( false );

        if ( numRanges_ == 0 ) break;
    } // ~for c

    // Clean up
    for ( auto rangeStore : rangeStoresA )
        delete rangeStore;
    for ( auto rangeStore : rangeStoresB )
        delete rangeStore;
    rangeStoresA.clear();
    rangeStoresB.clear();
    for ( auto & bwtReader : inBwtA_ )
    {
        delete bwtReader;
        bwtReader = 0;
    }
    for ( auto & bwtReader : inBwtB_ )
    {
        delete bwtReader;
        bwtReader = 0;
    }

    // Metagenomics-specific stuff
    if ( mode_ == BeetlCompareParameters::MODE_METAGENOMICS )
        releaseMetagomeMode();
} // countWords::run()


void CountWords::CountWords_parallelSubsetThread(
    const int threadNum
    , const int cycle
    , RangeStoreExternal &rangeStoreA
    , RangeStoreExternal &rangeStoreB
)
{
    string currentWord( cycle + 2, 'x' );
    int subsetThreadNum = threadNum / alphabetSize; // for ( int subsetThreadNum = 0; subsetThreadNum < subsetThreadCount; ++subsetThreadNum )
    int i = threadNum % alphabetSize; // for ( int i = 0; i < alphabetSize; ++i )
    /*
        {
            int pid = -1;
            int numThreads = -1;
            int processor = -1;
            readProcSelfStat( pid, numThreads, processor );
            #pragma omp critical (IO)
            cerr << "CountWords_parallelSubsetThread header cycle=" << cycle << " subsetThreadNum=" << subsetThreadNum << " pid=" << pid << " numThreads=" << numThreads << " processor=" << processor << endl;
            usleep( 100 );
        }

        omp_set_num_threads( alphabetSize );
        #pragma omp parallel for
        for ( int i = 0; i < alphabetSize; ++i )
    */
    {
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
        {
            int pid = -1;
            int numThreads = -1;
            int processor = -1;
            readProcSelfStat( pid, numThreads, processor );
            #pragma omp critical (IO)
            Logger::out() << "CountWords_parallelSubsetThread cycle=" << cycle << " subsetThreadNum=" << subsetThreadNum << " i=" << i << " pid=" << pid << " numThreads=" << numThreads << " processor=" << processor << endl;
        }

        if ( i == 0 ) return;
        if ( !isDistributedProcessResponsibleForPile( i ) )
            return;

        BwtReaderBase *inBwtA = inBwtA_[i]->clone();
        BwtReaderBase *inBwtB = inBwtB_[i]->clone();
        inBwtA->rewindFile();
        inBwtB->rewindFile();

        LetterCount countsSoFarA, countsSoFarB;
        LetterNumber currentPosA, currentPosB;

        RangeStoreExternal parallel_rA( rangeStoreA );
        RangeStoreExternal parallel_rB( rangeStoreB );

        inBwtA->rewindFile();
        inBwtB->rewindFile();
        currentPosA = 0;
        currentPosB = 0;
        countsSoFarA.clear();
        countsSoFarA += countsCumulativeA_[i - 1];
        countsSoFarB.clear();
        countsSoFarB += countsCumulativeB_[i - 1];

        for ( int j( 1 ); j < alphabetSize; ++j )
        {
            Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "positions - A: " << currentPosA << " B: " << currentPosB << endl;
            parallel_rA.setPortion( i, j );
            parallel_rB.setPortion( i, j );

            TwoBwtBackTracker backTracker( inBwtA, inBwtB
                                           , currentPosA, currentPosB
                                           , parallel_rA, parallel_rB, countsSoFarA, countsSoFarB
                                           , minOcc_, numCycles_, subset_, cycle + 1
                                           , doesPropagateBkptToSeqNumInSetA_
                                           , doesPropagateBkptToSeqNumInSetB_
                                           , noComparisonSkip_
                                           , propagateSequence_
                                         );
            switch ( mode_ )
            {
                case BeetlCompareParameters::MODE_REFERENCE:
                {
                    IntervalHandlerReference intervalHandler( minOcc_ );
                    backTracker.process( i, currentWord, intervalHandler );
                }
                break;
                case BeetlCompareParameters::MODE_METAGENOMICS:
                {
                    IntervalHandlerMetagenome intervalHandler( minOcc_, setC_, mmappedCFiles_, fileNumToTaxIds_, testDB_, minWordLen_, numCycles_ );
                    intervalHandler.createOutputFile( subsetThreadNum, i, j, cycle + 1, outputDirectory_ );
                    backTracker.process( i, currentWord, intervalHandler );
                }
                break;
                case BeetlCompareParameters::MODE_SPLICE:
                {
                    IntervalHandlerSplice intervalHandler( minOcc_ );
                    backTracker.process( i, currentWord, intervalHandler );
                }
                break;
                case BeetlCompareParameters::MODE_TUMOUR_NORMAL:
                {
                    IntervalHandlerTumourNormal intervalHandler( minOcc_, fsizeRatio_ );
                    intervalHandler.createOutputFile( subsetThreadNum, i, j, cycle + 1, outputDirectory_ );
                    backTracker.process( i, currentWord, intervalHandler );
                }
                break;
                default:
                    assert( false && "Unexpected mode" );
            }

            #pragma omp atomic
            numRanges_ += backTracker.numRanges_;
            #pragma omp atomic
            numSingletonRanges_ += backTracker.numSingletonRanges_;
            #pragma omp atomic
            numNotSkippedA_ += backTracker.numNotSkippedA_;
            #pragma omp atomic
            numNotSkippedB_ += backTracker.numNotSkippedB_;
            #pragma omp atomic
            numSkippedA_ += backTracker.numSkippedA_;
            #pragma omp atomic
            numSkippedB_ += backTracker.numSkippedB_;

            parallel_rA.deleteInputPortion( i, j );
            parallel_rB.deleteInputPortion( i, j );
        } // ~for j
        //cerr << "Done i " << i <<endl;
        parallel_rA.clear( false );
        parallel_rB.clear( false );


        delete inBwtA;
        delete inBwtB;
    }     // ~for i
}


/*
  @taxIdNames , file containing the taxonomic information about the database
  loads the taxonomic information so it is known which file number corresponds to which taxa
  fills up the vector fileNumToTaxIds;
  This should be called before the comparison starts
*/

void CountWords::loadFileNumToTaxIds( const string &taxIdNames )
{
    ifstream taxas( taxIdNames.c_str(), ios::in );
    assert( taxas.good() && "Error opening taxonomy file" );

    string line;
    //each line should contain the fileNumber followed by the taxIds split up with one space character
    while ( getline( taxas, line ) )
    {
        if ( line.empty() || line[0] == '#' )
            continue;

        string originalLine( line );
        if ( !line.empty() )
        {
            vector<int> taxIDs;
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
                cerr << "Tax Ids don't have enough taxonomic Information. Only " << taxIDs.size() << " could be found (" << originalLine << ")" << endl
                     << "Will add unknown taxa until size is right" << endl;
            else if ( taxIDs.size() > taxLevelSize )
                cerr << "Tax Ids have too much taxonomic information (" << originalLine << ")" << endl
                     << "Please note, that the taxonomic information about one file should be shown as: " << endl
                     << "FileNumber Superkingdom Phylum Order Family Genus Species Strain " << endl;
            taxIDs.resize( taxLevelSize );
            fileNumToTaxIds_.push_back( taxIDs );
            unsigned int test = fileNum + 1;
            if ( test != fileNumToTaxIds_.size() )
                cout << "Wrong filenumber " << fileNum << " " << fileNumToTaxIds_.size() << endl;
        }
    }
    //cout << " fineNumToTaxIds " << fileNumToTaxIds.size() <<endl;
}

void CountWords::initialiseMetagomeMode()
{
    loadFileNumToTaxIds( ncbiInfo_ );

    if ( ( *compareParams_ )["mmap C files"] == true )
    {
        mmappedCFiles_.resize( alphabetSize );
        for ( int i = 0; i < alphabetSize; i++ )
        {
            int fd = open( setC_[i].c_str(), O_RDONLY );
            if ( fd == -1 )
            {
                #pragma omp critical (IO)
                {
                    cout << "ERROR: Could not open File \"" << setC_[i] << "\"" << endl;
                }
            }
            assert( sizeof( off_t ) >= 8 && "64 bits architecture required to hold large files" );
            off_t fileSize = lseek( fd, 0, SEEK_END );
            mmappedCFiles_[i] = ( char * )mmap( NULL, fileSize, PROT_READ, MAP_SHARED , fd, 0 );
            if ( mmappedCFiles_[i] == ( void * ) - 1 )
            {
                perror ( "Error mmap " );
                mmappedCFiles_[i] = NULL;
            }
        }
    }
}

void CountWords::releaseMetagomeMode()
{
    for ( unsigned int i = 0; i < alphabetSize; i++ )
    {
        if ( mmappedCFiles_.size() > i && mmappedCFiles_[i] != NULL )
        {
            int fd = open( setC_[i].c_str(), O_RDONLY );
            off_t fileSize = lseek( fd, 0, SEEK_END );
            munmap( mmappedCFiles_[i], fileSize );
        }
    }
}

bool CountWords::isDistributedProcessResponsibleForPile( const int pile )
{
    if ( ( *compareParams_ )["4-way distributed"].isSet() )
    {
        // Delete own barrier file
        // Wait until all barrier files are deleted

        int processNum = ( *compareParams_ )["4-way distributed"];
        switch ( processNum )
        {
            case 0: // A
                if ( pile != 1 ) return false;
                break;
            case 1: // C+[^ACGNT]
                if ( pile != 2 && pile != 0 && pile <= 5 ) return false;
                break;
            case 2: // G+N
                if ( pile != 3 && pile != 4 ) return false;
                break;
            case 3: // T
                if ( pile != 5 ) return false;
                break;
        }

        #pragma omp critical (IO)
        cout << "process num = " << processNum << " processing pile " << pile << endl;
    }

    return true;
}
