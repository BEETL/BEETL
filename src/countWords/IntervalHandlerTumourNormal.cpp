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


#include "IntervalHandlerTumourNormal.hh"

#include "libzoo/util/Logger.hh"

using namespace std;


void IntervalHandlerTumourNormal::foundInBoth
( const int pileNum,
  const LetterCount &countsThisRangeA,
  const LetterCount &countsThisRangeB,
  const Range &thisRangeA,
  const Range &thisRangeB,
  AlphabetFlag &propagateIntervalA,
  AlphabetFlag &propagateIntervalB,
  bool &isBreakpointDetected,
  const int cycle
)
{
    int nonsharedPaths( 0 );
    int sharedPathsA( 0 );
    int sharedPathsB( 0 );
    int zeroPaths( 0 );
    double het_correction_factor( 2.0 );
    LetterNumber meanSignalAOnly( 0 ), meanSignalBOnly( 0 ), current_minOccA_( 0 ), current_minOccB_( 0 ), countsThisRangeA_cimulative( 0 ), countsThisRangeB_cimulative( 0 );

    if ( cycle < 12 )
    {
        current_minOccA_ =  minOcc_;
        current_minOccB_ =  minOcc_;
    }
    else
    {
        meanSignalAOnly = ( countsThisRangeA.count_[1] + countsThisRangeA.count_[2] + countsThisRangeA.count_[3] + countsThisRangeA.count_[5] ) / 10;
        meanSignalBOnly = ( countsThisRangeB.count_[1] + countsThisRangeB.count_[2] + countsThisRangeB.count_[3] + countsThisRangeB.count_[5] ) / 10;

        current_minOccA_ = max( meanSignalAOnly, minOcc_ );
        current_minOccB_ = max( meanSignalBOnly, minOcc_ );
    }

    for ( int l( 1 ); l < alphabetSize; l++ )
    {
        if ( l == 4 ) continue;
        sharedPathsB +=  countsThisRangeB.count_[l] > 1;
        sharedPathsA +=  countsThisRangeA.count_[l] > 1;
        nonsharedPaths += ( ( ( countsThisRangeA.count_[l] > ( LetterNumber )( ( double )current_minOccA_ * fsizeRatio_ ) && countsThisRangeB.count_[l] == 0 ) || ( countsThisRangeB.count_[l] > current_minOccB_ && countsThisRangeA.count_[l] == 0 ) ) );
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << l << " countsThisRangeA.count_ " << countsThisRangeA.count_[l] << "\t" << "countsThisRangeB.count_ "  << countsThisRangeB.count_[l] << endl;
        if ( countsThisRangeA.count_[0] == ( countsThisRangeA.count_[1] + countsThisRangeA.count_[2] + countsThisRangeA.count_[3] + countsThisRangeA.count_[5] ) ) nonsharedPaths = 0;
        if ( countsThisRangeB.count_[0] == ( countsThisRangeB.count_[1] + countsThisRangeB.count_[2] + countsThisRangeB.count_[3] + countsThisRangeB.count_[5] ) ) nonsharedPaths = 0;
    } // ~for l



    if (  nonsharedPaths > 0 && sharedPathsB < 3 && sharedPathsA < 3 )
    {
        for ( int l( 1 ); l < alphabetSize; l++ )
        {
            propagateIntervalA[l] = ( countsThisRangeA.count_[l] >= current_minOccA_ && countsThisRangeB.count_[l] == 0 );
            propagateIntervalB[l] = ( countsThisRangeB.count_[l] >= current_minOccB_ && countsThisRangeA.count_[l] == 0 );
        }

        isBreakpointDetected = true;
        #pragma omp critical (IO)
        Logger::out()
                << "BKPT"
#ifdef PROPAGATE_SEQUENCE
                << ' ' << thisRangeB.word_
#else
                // Print what we know of the sequence
                << ' ' << alphabet[pileNum] << string( cycle - 1, 'x' )
#endif
                << ' ' << countsThisRangeA.count_[0]
                << ':' << countsThisRangeA.count_[1]
                << ':' << countsThisRangeA.count_[2]
                << ':' << countsThisRangeA.count_[3]
                << ':' << countsThisRangeA.count_[4]
                << ':' << countsThisRangeA.count_[5]
                << ' ' << countsThisRangeB.count_[0]
                << ':' << countsThisRangeB.count_[1]
                << ':' << countsThisRangeB.count_[2]
                << ':' << countsThisRangeB.count_[3]
                << ':' << countsThisRangeB.count_[4]
                << ':' << countsThisRangeB.count_[5]
                << ' ' << ( thisRangeA.pos_ & matchMask )
                << ' ' << ( thisRangeB.pos_ & matchMask )
                << ' ' << thisRangeA.num_
                << ' ' << thisRangeB.num_
                << endl;
    }
    else
    {
        for ( int l( 1 ); l < alphabetSize; l++ )
            propagateIntervalB[l] = ( countsThisRangeB.count_[l] >= current_minOccB_ );
        for ( int l( 1 ); l < alphabetSize; l++ )
            propagateIntervalA[l] = ( countsThisRangeA.count_[l] >= current_minOccA_ );
    }

    // don't bother with Ns
    propagateIntervalA[whichPile[( int )dontKnowChar]] = false;
    propagateIntervalB[whichPile[( int )dontKnowChar]] = false;

} // ~foundInBoth


void IntervalHandlerTumourNormal::foundInAOnly
( const int pileNum,
  const LetterCount &countsSoFarA,
  const LetterCount &countsThisRangeA,
  const char *bwtSubstring,
  Range &thisRangeA,
  AlphabetFlag &propagateIntervalA,
  const int cycle
)
{
    if ( countsThisRangeA.count_[0] > 0 )
    {
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
        {
            #pragma omp critical (IO)
            Logger::out()
                    << "READ"
#ifdef PROPAGATE_SEQUENCE
                    << ' ' << thisRangeA.word_
#else
                    << ' ' << alphabet[pileNum]
#endif
                    << ' ' << thisRangeA.pos_
                    << ' ' << countsThisRangeA.count_[0]
                    << ':' << countsThisRangeA.count_[1]
                    << ':' << countsThisRangeA.count_[2]
                    << ':' << countsThisRangeA.count_[3]
                    << ':' << countsThisRangeA.count_[4]
                    << ':' << countsThisRangeA.count_[5]
                    << ' ' << countsSoFarA.count_[0]
                    << endl;
        }
    }
    // TBD print out IDs of discovered reads

    for ( int l( 1 ); l < alphabetSize; l++ )
    {
        propagateIntervalA[l] = ( countsThisRangeA.count_[l] > 0 );
    }

    // don't bother with Ns
    propagateIntervalA[whichPile[( int )dontKnowChar]] = false;
} // ~foundInBoth

void IntervalHandlerTumourNormal::foundInBOnly
( const int pileNum,
  const LetterCount &countsSoFarB,
  const LetterCount &countsThisRangeB,
  const char *bwtSubstring,
  Range &thisRangeB,
  AlphabetFlag &propagateIntervalB,
  const int cycle
)
{
    if ( countsThisRangeB.count_[0] > 0 )
    {
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
        {
            #pragma omp critical (IO)
            Logger::out()
                    << "INBS"
#ifdef PROPAGATE_SEQUENCE
                    << ' ' << thisRangeB.word_
#else
                    << ' ' << alphabet[pileNum]
#endif
                    << ' ' << thisRangeB.pos_
                    << ' ' << countsThisRangeB.count_[0]
                    << ':' << countsThisRangeB.count_[1]
                    << ':' << countsThisRangeB.count_[2]
                    << ':' << countsThisRangeB.count_[3]
                    << ':' << countsThisRangeB.count_[4]
                    << ':' << countsThisRangeB.count_[5]
                    << ' ' << countsSoFarB.count_[0]
                    << endl;
        }
    }
    // TBD print out IDs of discovered reads

    for ( int l( 1 ); l < alphabetSize; l++ )
    {
        propagateIntervalB[l] = ( countsThisRangeB.count_[l] > 0 );
    }

    // don't bother with Ns
    propagateIntervalB[whichPile[( int )dontKnowChar]] = false;
} // ~foundInBoth

