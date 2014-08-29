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

#include "IntervalHandlerReference.hh"

#include "libzoo/util/Logger.hh"

using namespace std;


//
// IntervalHandlerReference member function declarations
//

void IntervalHandlerReference::foundInBoth
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
    bool significantNonRef( false );
    //  LetterNumber maxSignalAOnly(0), maxSignalBOnly(0);

    if ( thisRangeB.num_ > 1 )
    {
        // if k-mer is not unique in B (reference),
        // propagate all B and all A that matches B
        for ( int l( 1 ); l < alphabetSize; l++ )
        {
            if ( countsThisRangeB.count_[l] > 0 )
            {
                propagateIntervalB[l] = true;
                if ( countsThisRangeA.count_[l] > 0 )
                    propagateIntervalA[l] = true;
                else
                    propagateIntervalA[l] = false;
            } // ~if
            else
            {
                propagateIntervalA[l] = false;
                propagateIntervalB[l] = false;
            } // ~else
        } // ~for
    } // ~if
    else
    {
        if ( thisRangeB.num_ != 1 )
        {
            cerr << "thisRandB is no 1. Aborting." << endl;
        }
        for ( int l( 1 ); l < alphabetSize; l++ )
        {
            if ( countsThisRangeB.count_[l] > 0 )
            {
                propagateIntervalB[l] = true;
                if ( countsThisRangeA.count_[l] > 0 )
                    propagateIntervalA[l] = true;
                else
                    propagateIntervalA[l] = false;
            } // ~if
            else
            {
                propagateIntervalB[l] = false;
                if ( countsThisRangeA.count_[l] > minOcc_ )
                {
                    propagateIntervalA[l] = true;
                    significantNonRef = true;
                }
                else
                    propagateIntervalA[l] = false;
            } // ~else
        } // ~for


    } // ~else

    if ( significantNonRef == true )
    {
        isBreakpointDetected = true;
        if ( !thisRangeB.word_.empty() )
        {
            #pragma omp critical (IO)
            Logger::out()
                    << "BKPT"
                    << ' ' << thisRangeB.word_
                    << ' ' << ( thisRangeB.pos_ & matchMask )
                    << ' ' << countsThisRangeA.count_[0]
                    << ':' << countsThisRangeA.count_[1]
                    << ':' << countsThisRangeA.count_[2]
                    << ':' << countsThisRangeA.count_[3]
                    << ':' << countsThisRangeA.count_[4]
                    << ':' << countsThisRangeA.count_[5]
                    << ':' << countsThisRangeB.count_[0]
                    << ':' << countsThisRangeB.count_[1]
                    << ':' << countsThisRangeB.count_[2]
                    << ':' << countsThisRangeB.count_[3]
                    << ':' << countsThisRangeB.count_[4]
                    << ':' << countsThisRangeB.count_[5]
                    << endl;
        }
        else
        {
            #pragma omp critical (IO)
            Logger::out()
                    << "BKPT"
                    << ' ' << alphabet[pileNum]
                    << ' ' << countsThisRangeA.count_[0]
                    << ':' << countsThisRangeA.count_[1]
                    << ':' << countsThisRangeA.count_[2]
                    << ':' << countsThisRangeA.count_[3]
                    << ':' << countsThisRangeA.count_[4]
                    << ':' << countsThisRangeA.count_[5]
                    << ':' << countsThisRangeB.count_[0]
                    << ':' << countsThisRangeB.count_[1]
                    << ':' << countsThisRangeB.count_[2]
                    << ':' << countsThisRangeB.count_[3]
                    << ':' << countsThisRangeB.count_[4]
                    << ':' << countsThisRangeB.count_[5]
                    << ' ' << ( thisRangeA.pos_ & matchMask )
                    << ' ' << ( thisRangeB.pos_ & matchMask )
                    << endl;
        }
    }

    // don't bother with Ns
    propagateIntervalA[whichPile[( int )dontKnowChar]] = false;
    propagateIntervalB[whichPile[( int )dontKnowChar]] = false;

} // ~foundInBoth


void IntervalHandlerReference::foundInAOnly
( const int pileNum,
  const LetterCount &countsSoFarA,
  const LetterCount &countsThisRangeA,
  const char *bwtSubstring,
  Range &thisRangeA,
  AlphabetFlag &propagateIntervalA,
  const int cycle
)
{

    bool significantPath( false );

    for ( int l( 1 ); l < alphabetSize; l++ )
    {
        if ( countsThisRangeA.count_[l] >= minOcc_ )
        {
            significantPath = true;
            propagateIntervalA[l] = true;
        } // ~if
        else
        {
            propagateIntervalA[l] = false;
        } // ~else

    } // ~for l

    if ( significantPath == false )
        #pragma omp critical (IO)
    {
        Logger::out() << "READ ";
        if ( thisRangeA.word_.empty() )
            Logger::out() << alphabet[pileNum]; // No propagated sequence
        else
            Logger::out() << thisRangeA.word_;
        Logger::out() << " " << thisRangeA.pos_;
        for ( int l( 0 ); l < alphabetSize; l++ )
            Logger::out() << ( ( l == 0 ) ? " " : ":" ) << countsThisRangeA.count_[l];
        Logger::out() << endl;
    }


#ifdef OLD
    // For now this is same as for Splice - continue until all reads found
    if ( countsThisRangeA.count_[0] > 0 )
        #pragma omp critical (IO)
    {
        Logger::out() << "READ " << thisRangeA.word_;
        Logger::out() << " " << thisRangeA.pos_;
        for ( int l( 0 ); l < alphabetSize; l++ )
            Logger::out() << ( ( l == 0 ) ? " " : ":" ) << countsThisRangeA.count_[l];
        Logger::out() << endl;
    }
    // TBD print out IDs of discovered reads

    for ( int l( 1 ); l < alphabetSize; l++ )
    {
        propagateIntervalA[l] = ( countsThisRangeA.count_[l] > 0 );
    } // ~for l
#endif

    // don't bother with Ns
    propagateIntervalA[whichPile[( int )dontKnowChar]] = false;
} // ~foundInBoth

void IntervalHandlerReference::foundInBOnly
( const int pileNum,
  const LetterCount &countsSoFarB,
  const LetterCount &countsThisRangeB,
  const char *bwtSubstring,
  Range &thisRangeB,
  AlphabetFlag &propagateIntervalB,
  const int cycle
)
{
    // TBD
}

