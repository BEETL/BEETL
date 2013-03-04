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

#include "IntervalHandlerReference.hh"

using namespace std;


//
// IntervalHandlerReference member function declarations
//

void IntervalHandlerReference::foundInBoth
( const int pileNum,
  const LetterCount &countsThisRangeA, const LetterCount &countsThisRangeB,
  const Range &thisRangeA, const Range &thisRangeB,
  AlphabetFlag &propagateIntervalA, AlphabetFlag &propagateIntervalB )
{
    bool significantNonRef( false );
    //  LetterCountType maxSignalAOnly(0), maxSignalBOnly(0);

    if ( thisRangeB.num_>1 )
    {
        // if k-mer is not unique in B (reference),
        // propagate all B and all A that matches B
        for ( int l( 1 ); l<alphabetSize; l++ )
        {
            if ( countsThisRangeB.count_[l]>0 )
            {
                propagateIntervalB[l]=true;
                if ( countsThisRangeA.count_[l]>0 )
                    propagateIntervalA[l]=true;
                else
                    propagateIntervalA[l]=false;
            } // ~if
            else
            {
                propagateIntervalA[l]=false;
                propagateIntervalB[l]=false;
            } // ~else
        } // ~for
    } // ~if
    else
    {
        if ( thisRangeB.num_!=1 )
        {
            cerr << "thisRandB is no 1. Aborting." << endl;
        }
        for ( int l( 1 ); l<alphabetSize; l++ )
        {
            if ( countsThisRangeB.count_[l]>0 )
            {
                propagateIntervalB[l]=true;
                if ( countsThisRangeA.count_[l]>0 )
                    propagateIntervalA[l]=true;
                else
                    propagateIntervalA[l]=false;
            } // ~if
            else
            {
                propagateIntervalB[l]=false;
                if ( countsThisRangeA.count_[l]>minOcc_ )
                {
                    propagateIntervalA[l]=true;
                    significantNonRef=true;
                }
                else
                    propagateIntervalA[l]=false;
            } // ~else
        } // ~for


    } // ~else

    if ( significantNonRef==true )
    {
#ifdef PROPAGATE_PREFIX
        printf(
            "BKPT %s %llu %llu:%llu:%llu:%llu:%llu:%llu %llu:%llu:%llu:%llu:%llu:%llu\n",
            thisRangeB.word_.c_str(),
            thisRangeB.pos_&matchMask,
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
            thisRangeA.pos_,
            thisRangeB.pos_ );
#endif

#ifdef OLD
        cout << "BKPT " << ( thisRangeB.num_&matchMask );
        cout << " " << thisRangeB.pos_;
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

#ifdef OLD
{
    for ( int l( 1 ); l<alphabetSize; l++ )
    {
        propagateIntervalA[l]=( countsThisRangeA.count_[l]>=minOcc_ );
        propagateIntervalB[l]=( countsThisRangeB.count_[l]>0 );
    } // ~for l
} // ~foundInBoth
#endif



void IntervalHandlerReference::foundInAOnly
( const int pileNum,
  const LetterCount &countsSoFarA,
  const LetterCount &countsThisRangeA,
  const Range &thisRangeA,
  AlphabetFlag &propagateIntervalA )
{

    bool significantPath( false );

    for ( int l( 1 ); l<alphabetSize; l++ )
    {
        if ( countsThisRangeA.count_[l]>=minOcc_ )
        {
            significantPath=true;
            propagateIntervalA[l]=true;
        } // ~if
        else
        {
            propagateIntervalA[l]=false;
        } // ~else

    } // ~for l

    if ( significantPath==false )
    {
        cout << "READ ";
#ifdef PROPAGATE_PREFIX
        cout << thisRangeA.word_;
#endif
        cout << " " << thisRangeA.pos_;
        for ( int l( 0 ); l<alphabetSize; l++ )
            cout << ( ( l==0 )?" ":":" ) << countsThisRangeA.count_[l];
        cout << endl;
    }


#ifdef OLD
    // For now this is same as for Splice - continue until all reads found
    if ( countsThisRangeA.count_[0]>0 )
    {
        cout << "READ " << thisRangeA.word_;
        cout << " " << thisRangeA.pos_;
        for ( int l( 0 ); l<alphabetSize; l++ )
            cout << ( ( l==0 )?" ":":" ) << countsThisRangeA.count_[l];
        cout << endl;
    }
    // TBD print out IDs of discovered reads

    for ( int l( 1 ); l<alphabetSize; l++ )
    {
        propagateIntervalA[l]=( countsThisRangeA.count_[l]>0 );
    } // ~for l
#endif

    // don't bother with Ns
    propagateIntervalA[whichPile[( int )dontKnowChar]]=false;
} // ~foundInBoth

void IntervalHandlerReference::foundInBOnly
( const int pileNum,
  const LetterCount &countsSoFarB,
  const LetterCount &countsThisRangeB,
  const Range &thisRangeB,
  AlphabetFlag &propagateIntervalB )
{
    // TBD
}

