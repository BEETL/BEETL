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


#include "IntervalHandlerSplice.hh"

using namespace std;


//
// IntervalHandlerSplice member function declarations
//

void IntervalHandlerSplice::foundInBoth
( const int pileNum,
  const LetterCount &countsThisRangeA, const LetterCount &countsThisRangeB,
  const Range &thisRangeA, const Range &thisRangeB,
  AlphabetFlag &propagateIntervalA, AlphabetFlag &propagateIntervalB )
{
    bool sharedPath( false );
    LetterCountType maxSignalAOnly( 0 ), maxSignalBOnly( 0 );

    for ( int l( 1 ); l<alphabetSize; l++ )
    {
        if ( ( countsThisRangeB.count_[l]==0 )
             &&( countsThisRangeA.count_[l]>maxSignalAOnly ) )
            maxSignalAOnly=countsThisRangeA.count_[l];

        if ( ( countsThisRangeA.count_[l]==0 )
             &&( countsThisRangeB.count_[l]>maxSignalBOnly ) )
            maxSignalBOnly=countsThisRangeB.count_[l];

        sharedPath|=( ( countsThisRangeA.count_[l]>0 )
                      &&( countsThisRangeB.count_[l]>0 ) );
        // for a general shared interval, we only care about backward
        // extensions in B (genome) that are common to A (RNA)
        propagateIntervalA[l]=( countsThisRangeA.count_[l]>=minOcc_ );
        propagateIntervalB[l]=( propagateIntervalA[l]&&
                                ( countsThisRangeB.count_[l]>=0 ) );
    } // ~for l
    if ( ( sharedPath==false )
         &&( maxSignalAOnly>=minOcc_ )
         &&( maxSignalBOnly>=minOcc_ ) )
    {
        // At a breakpoint position, we *do* care about the B-only intervals
        // extensions in B (genome) that are common to A (RNA)
        for ( int l( 1 ); l<alphabetSize; l++ )
            propagateIntervalB[l]=( countsThisRangeB.count_[l]>=minOcc_ );


#ifdef PROPAGATE_PREFIX
        printf(
            "BKPT %s %llu:%llu:%llu:%llu:%llu:%llu %llu:%llu:%llu:%llu:%llu:%llu\n",
            thisRangeB.word_.c_str(),
            countsThisRangeA.count_[0],
            countsThisRangeA.count_[1],
            countsThisRangeA.count_[2],
            countsThisRangeA.count_[3],
            countsThisRangeA.count_[4],
            countsThisRangeA.count_[5],
            countsThisRangeB.count_[0],
            countsThisRangeB.count_[1],
            countsThisRangeB.count_[2],
            countsThisRangeB.count_[3],
            countsThisRangeB.count_[4],
            countsThisRangeB.count_[5] );
#else
        printf(
            "BKPT %c %llu:%llu:%llu:%llu:%llu:%llu %llu:%llu:%llu:%llu:%llu:%llu %llu %llu\n",
            alphabet[pileNum],
            countsThisRangeA.count_[0],
            countsThisRangeA.count_[1],
            countsThisRangeA.count_[2],
            countsThisRangeA.count_[3],
            countsThisRangeA.count_[4],
            countsThisRangeA.count_[5],
            countsThisRangeB.count_[0],
            countsThisRangeB.count_[1],
            countsThisRangeB.count_[2],
            countsThisRangeB.count_[3],
            countsThisRangeB.count_[4],
            countsThisRangeB.count_[5],
            ( thisRangeA.pos_&matchMask ),
            ( thisRangeB.pos_&matchMask ) );
#endif

#ifdef OLD
        cout << "BKPT ";
#ifdef PROPAGATE_PREFIX
        cout << thisRangeA.word_;
#endif
        for ( int l( 0 ); l<alphabetSize; l++ )
            cout << ( ( l==0 )?" ":":" ) << countsThisRangeA.count_[l];
        for ( int l( 0 ); l<alphabetSize; l++ )
            cout << ( ( l==0 )?" ":":" ) << countsThisRangeB.count_[l];
        cout << endl;
#endif
    }

    // don't bother with Ns
    propagateIntervalA[whichPile[( int )dontKnowChar]]=false;
    propagateIntervalB[whichPile[( int )dontKnowChar]]=false;

} // ~foundInBoth

void IntervalHandlerSplice::foundInAOnly
( const int pileNum,
  const LetterCount &countsSoFarA,
  const LetterCount &countsThisRangeA,
  const Range &thisRangeA,
  AlphabetFlag &propagateIntervalA )
{
    if ( countsThisRangeA.count_[0]>0 )
    {
#ifdef PROPAGATE_PREFIX
        printf( "READ %s %llu %llu:%llu:%llu:%llu:%llu:%llu %llu\n",
                thisRangeA.word_.c_str(),
                thisRangeA.pos_,
                countsThisRangeA.count_[0],
                countsThisRangeA.count_[1],
                countsThisRangeA.count_[2],
                countsThisRangeA.count_[3],
                countsThisRangeA.count_[4],
                countsThisRangeA.count_[5],
                countsSoFarA.count_[0] );
#else
        printf( "READ %c %llu %llu:%llu:%llu:%llu:%llu:%llu %llu\n",
                alphabet[pileNum],
                thisRangeA.pos_,
                countsThisRangeA.count_[0],
                countsThisRangeA.count_[1],
                countsThisRangeA.count_[2],
                countsThisRangeA.count_[3],
                countsThisRangeA.count_[4],
                countsThisRangeA.count_[5],
                countsSoFarA.count_[0] );
#endif

#ifdef OLD
#ifdef PROPAGATE_PREFIX
        cout << "READ " << thisRangeA.word_;
#else
        cout << "READ " << thisRangeA.num_;
#endif
        cout << " " << thisRangeA.pos_;
        for ( int l( 0 ); l<alphabetSize; l++ )
            cout << ( ( l==0 )?" ":":" ) << countsThisRangeA.count_[l];
        cout << endl;
#endif
    }
    // TBD print out IDs of discovered reads

    for ( int l( 1 ); l<alphabetSize; l++ )
    {
        propagateIntervalA[l]=( countsThisRangeA.count_[l]>0 );
    } // ~for l

    // don't bother with Ns
    propagateIntervalA[whichPile[( int )dontKnowChar]]=false;
} // ~foundInBoth

void IntervalHandlerSplice::foundInBOnly
( const int pileNum,
  const LetterCount &countsSoFarB,
  const LetterCount &countsThisRangeB,
  const Range &thisRangeB,
  AlphabetFlag &propagateIntervalB )
{
    if ( countsThisRangeB.count_[0]>0 )
    {
#ifdef PROPAGATE_PREFIX
        printf( "INBS %s %llu %llu:%llu:%llu:%llu:%llu:%llu %llu\n",
                thisRangeB.word_.c_str(),
                thisRangeB.pos_,
                countsThisRangeB.count_[0],
                countsThisRangeB.count_[1],
                countsThisRangeB.count_[2],
                countsThisRangeB.count_[3],
                countsThisRangeB.count_[4],
                countsThisRangeB.count_[5],
                countsSoFarB.count_[0] );
#else
        printf( "INBS %c %llu %llu:%llu:%llu:%llu:%llu:%llu %llu\n",
                alphabet[pileNum],
                thisRangeB.pos_,
                countsThisRangeB.count_[0],
                countsThisRangeB.count_[1],
                countsThisRangeB.count_[2],
                countsThisRangeB.count_[3],
                countsThisRangeB.count_[4],
                countsThisRangeB.count_[5],
                countsSoFarB.count_[0] );
#endif

#ifdef OLD
#ifdef PROPAGATE_PREFIX
        cout << "READ " << thisRangeB.word_;
#else
        cout << "READ " << thisRangeB.num_;
#endif
        cout << " " << thisRangeB.pos_;
        for ( int l( 0 ); l<alphabetSize; l++ )
            cout << ( ( l==0 )?" ":":" ) << countsThisRangeB.count_[l];
        cout << endl;
#endif
    }
    // TBD print out IDs of discovered reads

    for ( int l( 1 ); l<alphabetSize; l++ )
    {
        propagateIntervalB[l]=( countsThisRangeB.count_[l]>0 );
    } // ~for l

    // don't bother with Ns
    propagateIntervalB[whichPile[( int )dontKnowChar]]=false;
} // ~foundInBoth

