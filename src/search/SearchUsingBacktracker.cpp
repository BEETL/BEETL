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

#include "SearchUsingBacktracker.hh"

#include "BwtReader.hh"
#include "IntervalFile.hh"
#include "KmerSearchIntervalHandler.hh"
#include "KmerSearchRange.hh"
#include "OneBwtBackTracker.hh"
#include "RangeStore.hh"
#include "Timer.hh"
#include "Tools.hh"
#include "config.h"
#include "parameters/SearchParameters.hh"
#include "libzoo/util/Logger.hh"

#include <algorithm>
#include <sstream>

#ifdef _OPENMP
# include <omp.h>
#endif //ifdef _OPENMP

using namespace std;


vector<KmerSearchItem> kmerList2;


SearchUsingBacktracker::SearchUsingBacktracker(
    const SearchParameters &searchParams
)
    : searchParams_( searchParams )
{
}

void SearchUsingBacktracker::run()
{
    Timer  timer;
    vector <BwtReaderBase *> inBwt( alphabetSize );

    LetterCountEachPile countsPerPile, countsCumulative;

    RangeStoreExternal r;

    string bwtPrefix = searchParams_["input"];
    string subset_ = "";
    EndPosFile endPosFile( bwtPrefix );

    vector<string> kmerList;
    if ( searchParams_["one kmer string"].isSet() )
    {
        kmerList.push_back( searchParams_["one kmer string"] );
    }
    else
    {
        assert( searchParams_["kmers input file"].isSet() );
        string filename = searchParams_["kmers input file"];
        ifstream ifs( filename );
        string line;
        // get 1st word of each line as kmer
        while ( getline( ifs, line ) )
        {
            istringstream iss( line );
            string kmer;
            iss >> kmer;
            assert( kmer.size() != 1 && "todo: 1-mer search" );
            if ( kmer.size() >= 2 )
                kmerList.push_back( kmer );
        }
    }

    SequenceNumber originalIndex = 0;
    for ( auto kmer : kmerList )
    {
        std::reverse( kmer.begin(), kmer.end() );
        kmerList2.push_back( KmerSearchItem( kmer, 0, 0, originalIndex++ ) );
    }
    std::sort( kmerList2.begin(), kmerList2.end() );
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
    {
        Logger::out() << "kmerList2:" << endl;
        for ( auto kmerItem : kmerList2 )
            Logger::out() << "  " << kmerItem.kmer << endl;
    }

    inBwt = instantiateBwtPileReaders( bwtPrefix, searchParams_.getStringValue( "use shm" ) );

    for ( int i( 0 ); i < alphabetSize; i++ )
    {
        inBwt[i]->readAndCount( countsPerPile[i] );
    }
    countsCumulative = countsPerPile;
    for ( int i( 1 ); i < alphabetSize; i++ )
    {
        countsCumulative[i] += countsCumulative[i - 1];
    }

    Logger_if( LOG_SHOW_IF_VERBOSE )
    {
        cout << "countsPerPile:" << endl;
        countsPerPile.print();
        cout << "countsCumulative:" << endl;
        countsCumulative.print();
    }

    string currentWord_notUsed;
    int startPosInKmerList2 = 0;
    int endPosInKmerList2 = 0;

    for ( int i( 1 ); i < alphabetSize; ++i )
    {
        for ( int j( 1 ); j < alphabetSize; ++j )
        {
            startPosInKmerList2 = endPosInKmerList2;
            while ( endPosInKmerList2 < ( int )kmerList2.size() &&
                    kmerList2[endPosInKmerList2].kmer[0] == alphabet[i] &&
                    kmerList2[endPosInKmerList2].kmer[1] == alphabet[j] )
                ++endPosInKmerList2;

            if ( countsPerPile[i].count_[j] != 0 )
            {
                if ( startPosInKmerList2 != endPosInKmerList2 )
                    r.addRange(
                        KmerSearchRange(
                            currentWord_notUsed,
                            countsCumulative[i - 1].count_[j],
                            countsPerPile[i].count_[j],
                            false
                            , startPosInKmerList2
                            , endPosInKmerList2
                        )
                        , j,
                        i,
                        subset_,
                        1
                    );
            }
            else
            {
                if ( startPosInKmerList2 != endPosInKmerList2 )
                {
                    cout << "Oh nooo, I can't find: " << endl;
                    for ( int k = startPosInKmerList2; k < endPosInKmerList2; ++k )
                    {
                        string rev = kmerList2[k].kmer;
                        std::reverse( rev.begin(), rev.end() );
                        cout << "  " << kmerList2[k].kmer << " (" << rev << ")" << endl;
                    }
                }
            }
        } // ~for j
    } // ~for i

    LetterCount countsSoFar;
    LetterNumber currentPos;
    LetterNumber numRanges, numSingletonRanges;

    r.clear();
    for ( int cycle( 1 ); ; ++cycle )
    {
        Logger_if( LOG_LEVEL_NORMAL )
        Logger::out() << "cycle: " << cycle << endl;

        Logger_if( LOG_SHOW_IF_VERBOSE )
        {
            Logger::out() << "   time now: " << timer.timeNow();
            Logger::out() << "   usage: " << timer << endl;
        }

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

                const bool propagateSequence = false;//( compareParams_ ? ( *compareParams_ )["propagate sequence"] : false );
                OneBwtBackTracker backTracker(
                    inBwt[i],
                    currentPos,
                    r,
                    countsSoFar,
                    subset_,
                    cycle + 1,
                    false, // don't propagate breakpoints until $ sign
                    true, // skip-already-processed-intervals deactivated
                    propagateSequence,
                    endPosFile
                );

                KmerSearchIntervalHandler intervalHandler;
                KmerSearchRange rangeObject;
                backTracker.process( i, currentWord_notUsed, intervalHandler, rangeObject );

                numRanges += backTracker.numRanges_;
                numSingletonRanges += backTracker.numSingletonRanges_;

            } // ~for j
            //cerr << "Done i " << i <<endl;
        }     // ~for i
        r.clear();
        //    return 0; // %%%
        Logger_if( LOG_SHOW_IF_VERBOSE )
        {
            Logger::out() << "Finished cycle " << cycle << ": ranges=" << numRanges << " singletons=" << numSingletonRanges << endl;
        }

        if ( numRanges == 0 ) break;

    } // ~for c

    for ( int i = 0; i < alphabetSize; i++ )
        delete inBwt[i];

    // Output
    {
        ostream *outputStreamPtr = &std::cout;
        string outputFilename = searchParams_["output"];
        ofstream ofs;
        if ( outputFilename != "-" )
        {
            ofs.open( outputFilename );
            if ( ofs.good() )
                outputStreamPtr = &ofs;
            else
                cerr << "Warning: Couldn't open output file " << outputFilename << ". Sending output to stdout." << endl;
        }

        IntervalWriter writer( *outputStreamPtr );
        for ( auto kmerItem : kmerList2 )
        {
            IntervalRecord rec( kmerList[kmerItem.originalIndex], kmerItem.position, kmerItem.count );
            writer.write( rec );
        }
    }
}
