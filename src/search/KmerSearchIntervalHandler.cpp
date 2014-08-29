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

#include "KmerSearchIntervalHandler.hh"

#include "KmerSearchRange.hh"
#include "libzoo/util/Logger.hh"

#include <algorithm>
#include <numeric>

using namespace std;

extern vector<KmerSearchItem> kmerList2;


KmerSearchIntervalHandler::KmerSearchIntervalHandler()
{
}

void KmerSearchIntervalHandler::foundInAOnly
( const int pileNum,
  const LetterCount &countsSoFarA,
  const LetterCount &countsThisRangeA,
  const char *bwtSubstring,
  Range &thisRangeBaseA,
  AlphabetFlag &propagateIntervalA,
  const int cycle
)
{
    KmerSearchRange &thisRangeA = dynamic_cast< KmerSearchRange & >( thisRangeBaseA );
    thisRangeA.clearDataForSubIntervals();

    for ( int i = 0; i < alphabetSize; i++ )
        propagateIntervalA[i] = false;

    assert( thisRangeA.data_.start < thisRangeA.data_.end );
    int lastPile = 0;
    int lastPileEnd = thisRangeA.data_.start;

    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "this interval's kmers:" << endl;
    for ( int k = thisRangeA.data_.start; k < thisRangeA.data_.end; ++k )
    {
        KmerSearchItem &kmerSearchItem = kmerList2[k];
        string &kmer = kmerSearchItem.kmer;
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "  " << kmer << endl;
        assert( kmer[cycle - 1] == alphabet[pileNum] );
        if ( ( int )kmer.size() == cycle )
        {
            LetterNumber total = std::accumulate( &countsThisRangeA.count_[0], &countsThisRangeA.count_[alphabetSize], ( LetterNumber )0 );
            assert( total == thisRangeBaseA.num_ );
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "    -> FOUND " << total << " " << kmer << endl;
            //            LetterNumber totalSoFar = std::accumulate( &countsSoFarA.count_[0], &countsSoFarA.count_[alphabetSize], 0 );
            kmerSearchItem.position = thisRangeBaseA.pos_;
            kmerSearchItem.count = thisRangeBaseA.num_;
            assert( lastPile == 0 );
        }
        else
        {
            int pile = whichPile[( int )kmer[cycle]];
            propagateIntervalA[pile] = true;

            if ( pile != lastPile )
            {
                thisRangeA.getDataForSubInterval( pile ).start = lastPileEnd;
            }
            thisRangeA.getDataForSubInterval( pile ).end = k + 1;
            lastPile = pile;
        }
        lastPileEnd = k + 1;
    }

    for ( int i = 0; i < alphabetSize; i++ )
        if ( propagateIntervalA[i] )
        {
            assert( thisRangeA.getDataForSubInterval( i ).start < thisRangeA.getDataForSubInterval( i ).end );
        }

}
