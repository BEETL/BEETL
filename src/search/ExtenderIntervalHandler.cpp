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

#include "ExtenderIntervalHandler.hh"

#include "IntervalFile.hh"
#include "libzoo/util/ColorText.hh"
#include "libzoo/util/Logger.hh"

#include <algorithm>

using namespace std;

extern vector<string> kmerList2;


ExtenderIntervalHandler::ExtenderIntervalHandler( EndPosFile &endPosFile )
    : endPosFile_( endPosFile )
{
}

void ExtenderIntervalHandler::foundInAOnly
( const int pileNum,
  const LetterCount &countsSoFarA,
  const LetterCount &countsThisRangeA,
  const char *bwtSubstring,
  Range &thisRangeBaseA,
  AlphabetFlag &propagateIntervalA,
  const int cycle
)
{
    if ( countsThisRangeA.count_[0] )
    {
        IntervalRecord *rec = reinterpret_cast< IntervalRecord * >( thisRangeBaseA.userData_ );

        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "$ signs detected for " << *rec << ": " << countsThisRangeA.count_[0] << " items from " << countsSoFarA.count_[0] << endl;
        if ( !thisRangeBaseA.word_.empty() )
        {
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "  Sub-sequence from beginning of read to searched kmer: " << thisRangeBaseA.word_ << endl;
        }

        for ( LetterNumber i = 0; i < countsThisRangeA.count_[0]; ++i )
            rec->dollarSignPositions.push_back( countsSoFarA.count_[0] + i );
    }

    for ( int l( 0 ); l < alphabetSize; l++ )
        propagateIntervalA[l] = ( countsThisRangeA.count_[l] > 0 );

    Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
    {
        IntervalRecord *rec = reinterpret_cast< IntervalRecord * >( thisRangeBaseA.userData_ );
        Logger::out() << "Propagating " << *rec << " to " << countsThisRangeA << endl;
    }

    if ( !thisRangeBaseA.word_.empty() ) // == if --propagate-sequence
    {
        // Full sequence output at last cycle
        const IntervalRecord *rec = reinterpret_cast< IntervalRecord * >( thisRangeBaseA.userData_ );
        const LetterNumber nextBwtPosStart = countsSoFarA.count_[whichPile[( int )rec->kmer[0]]];
        const LetterNumber count = countsThisRangeA.count_[whichPile[( int )rec->kmer[0]]];
        if ( ( count && ( nextBwtPosStart + count > rec->position ) && ( nextBwtPosStart < rec->position + rec->count ) )
             || ( rec->kmer[0] == '$' && countsThisRangeA.count_[0] > 0 ) )
            //    const int cycleCount = 101;
            //    if (cycle >= cycleCount + 1)
        {
            const int cycleCount = cycle - 1;
            const string &rotatedSeq = thisRangeBaseA.word_;
            //        cout << "rec " << *rec << ":\nSeq=" << rotatedSeq << endl;

            const size_t searchedKmerSize = rec->kmer.size() + 1; // TODO: this k-mer size should be fixed in beetl-compare output
            const size_t dollarPos = rotatedSeq.find( '$' );
            size_t kmerPosInRead;
            string part1, part2, part3;
            //string part1_cycle0Letter;
            if ( rec->kmer[0] == '$' )
            {
                kmerPosInRead = 0;
                part1 = rotatedSeq.substr( 0, dollarPos  );
            }
            else if ( dollarPos != string::npos )
            {
                kmerPosInRead = cycleCount - dollarPos;
                part1 = rotatedSeq.substr( dollarPos + 1, kmerPosInRead - 2 );
                //                part1_cycle0Letter = rotatedSeq.substr( dollarPos + kmerPosInRead - 1, 1 );
                part2 = rotatedSeq.substr( 0, searchedKmerSize - 2 );
                if ( dollarPos + 2 >= searchedKmerSize )
                    part3 = rotatedSeq.substr( searchedKmerSize - 2, dollarPos - searchedKmerSize + 2 );
                else
                    part3 = "??? " + rotatedSeq;
            }
            else
            {
                kmerPosInRead = cycleCount + 1;
                part1 = rotatedSeq;
            }

            // We still have 1 character to propagate to get the full sequence, and the possible values of this character are in countsThisRangeA
            for ( int i = 0; i < alphabetSize; i++ )
            {
                for ( uint j = 0; j < countsThisRangeA.count_[i]; ++j ) // Just repeat the same line multiple times in the unlikely event of this count being > 1
                {
                    string lastPropagatedChar;
                    if ( i ) lastPropagatedChar = alphabet[i]; // don't print '$' chars

                    // header
                    static int firstTime = true;
                    if ( firstTime )
                    {
                        cout << "Output:\t//kmer\tposition\tcount\tposInRead\tdollarPos\tseqNum\tseq" << endl;
                        firstTime = false;
                    }

                    SequenceNumber dollarPos = 0, seqNum = 0;
                    if ( !rec->dollarSignPositions.empty() )
                    {
                        dollarPos = rec->dollarSignPositions[0];
                        seqNum = endPosFile_.convertDollarNumToSequenceNum( dollarPos );
                    }

                    cout << "Output:\t" << rec->kmer
                         << '\t' << rec->position
                         << '\t' << rec->count
                         << '\t' << ( ( off_t )kmerPosInRead - 1 )
                         << '\t' << dollarPos
                         << '\t' << seqNum
                         << '\t' << part1 << ColorText::startRed << /*part1_cycle0Letter << */ lastPropagatedChar << part2 << ColorText::endRed << part3 << endl;
                }
            }

            // End propagation
            for ( int i = 0; i < alphabetSize; i++ )
                propagateIntervalA[i] = false;
        }
    }
    else
    {
        // if not --propagate-sequence, then we stop the propagation at $ signs
        propagateIntervalA[0] = false;
    }
}
