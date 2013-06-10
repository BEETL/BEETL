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

#ifndef INCLUDED_BACKTRACKER_HH
#define INCLUDED_BACKTRACKER_HH

#include "BwtReader.hh"
#include "IntervalHandlerBase.hh"
#include "LetterCount.hh"
#include "RangeStore.hh"
#include "Types.hh"

#include <string>


typedef bool AlphabetFlag[alphabetSize];

// BackTracker class implements the backward search and takes its
// cue from IntervalHandler as to whether to continue down a particular
// branch of the search tree
struct BackTracker
{
    BackTracker( BwtReaderBase *inBwtA, BwtReaderBase *inBwtB,
                 LetterCountType &currentPosA, LetterCountType &currentPosB,
                 RangeStoreExternal &rA, RangeStoreExternal &rB,
                 LetterCount &countsSoFarA, LetterCount &countsSoFarB,
                 int minOcc, const int maxLength, const string &subset );

    void skipIfNecessary( const Range &thisRange,
                          LetterCountType &currentPos,
                          BwtReaderBase &inBwt,
                          LetterCount &countsSoFar );

    template<class IntervalHandlerParam>  void operator()
    ( int i, string &thisWord, IntervalHandlerParam &intervalHandler );

    BwtReaderBase *inBwtA_;
    BwtReaderBase *inBwtB_;

    LetterCountType &currentPosA_;
    LetterCountType &currentPosB_;

    RangeStoreExternal &rA_;
    RangeStoreExternal &rB_;
    LetterCount &countsSoFarA_;
    LetterCount &countsSoFarB_;

    int minOcc_;
    int maxLength_;
    const string &subset_;

    LetterCountType numRanges_;
    LetterCountType numSingletonRanges_;

    AlphabetFlag propagateIntervalA_;
    AlphabetFlag propagateIntervalB_;
};

template<class IntervalHandlerParam>
void BackTracker::operator()
( int i, string &thisWord, IntervalHandlerParam &intervalHandler_  )
{
    LetterCount countsThisRangeA, countsThisRangeB;
    Range thisRangeA, thisRangeB;
    //  string thisWord;
    bool notAtLastA( true ), notAtLastB( true );
    bool hasChild;

    while ( 1 )
    {
        while ( notAtLastA )
        {
            notAtLastA = rA_.getRange( thisRangeA );
            if ( ( notAtLastA == false ) || ( ( thisRangeA.pos_ & matchFlag ) != 0 ) ) break;

#ifdef DEBUG
            cout << "RangeA: " << i << " " << j << " "
                 << thisRangeA.word_ << " " << thisRangeA.pos_ << " " << thisRangeA.num_
                 << " -- " << currentPosA_ << " < " << thisRangeB.word_ << endl;
#endif

            skipIfNecessary( thisRangeA, currentPosA_, ( *inBwtA_ ), countsSoFarA_ );
            // count children
            countsThisRangeA.clear();
            inBwtA_->readAndCount( countsThisRangeA, thisRangeA.num_ );
#ifdef DEBUG
            countsThisRangeA.print();
#endif
            intervalHandler_.foundInAOnly
            ( i,
              countsSoFarA_,
              countsThisRangeA,
              thisRangeA,
              propagateIntervalA_ );


            // add ranges for any children


#ifdef PROPAGATE_PREFIX
            hasChild = false;
#endif

            for ( int l( 1 ); l < alphabetSize; l++ )
            {
                //       if (countsThisRangeA.count_[l]>=minOcc)
                if ( propagateIntervalA_[l] == true )
                {
#ifdef PROPAGATE_PREFIX
                    if ( hasChild == false )
                    {
                        // assert(thisWord.size()==thisRangeA.word_.size()+1);
                        thisWord.replace( 1, thisRangeA.word_.size(), thisRangeA.word_ );
                    } // ~if
                    thisWord[0] = alphabet[l];
#endif
                    hasChild = true;

                    rA_.addRange( l, i, thisWord,
                                  countsSoFarA_.count_[l],
                                  countsThisRangeA.count_[l],
                                  subset_ );
                } // ~if
            } // ~for l

            if ( hasChild == false )
            {
                //  if no children, print word itself
#ifdef OLD
                cout << "GOLD ";
#ifdef PROPAGATE_PREFIX
                cout << thisRangeA.word_;
#endif
                cout << " " << thisRangeA.num_ << endl;
#endif
                numSingletonRanges_++;
            } // ~if

            countsSoFarA_ += countsThisRangeA;
            currentPosA_ += thisRangeA.num_;

            numRanges_++;
        } // ~while notAtLastA

        while ( notAtLastB )
        {
            notAtLastB = rB_.getRange( thisRangeB );
            if ( ( notAtLastB == false ) || ( ( thisRangeB.pos_ & matchFlag ) != 0 ) ) break;

#ifdef DEBUG
            cout << "RangeB: " << i << " " << j << " "
                 << thisRangeB.word_ << " " << thisRangeB.pos_ << " " << thisRangeB.num_
                 << " -- " << currentPosB_ << " < " << thisRangeB.word_ << endl;
#endif

            skipIfNecessary( thisRangeB, currentPosB_, ( *inBwtB_ ), countsSoFarB_ );
            // count children
            countsThisRangeB.clear();
            uint charsRead = inBwtB_->readAndCount( countsThisRangeB, thisRangeB.num_ );
            assert( charsRead == thisRangeB.num_ );
#ifdef DEBUG
            countsThisRangeB.print();
#endif

            //#ifdef TBD
            intervalHandler_.foundInBOnly
            ( i,
              countsSoFarB_,
              countsThisRangeB,
              thisRangeB,
              propagateIntervalB_ );
            //#endif

            // add ranges for any children
#ifdef PROPAGATE_PREFIX
            hasChild = false;
#endif
            for ( int l( 1 ); l < alphabetSize; l++ )
            {
                //       if (countsThisRangeB.count_[l]>=minOcc)
                if ( propagateIntervalB_[l] == true )
                {
#ifdef PROPAGATE_PREFIX
                    if ( hasChild == false )
                    {
                        // assert(thisWord.size()==thisRangeB.word_.size()+1);
                        thisWord.replace( 1, thisRangeB.word_.size(), thisRangeB.word_ );
                        hasChild = true;
                    } // ~if
                    thisWord[0] = alphabet[l];
#endif

                    rB_.addRange( l, i, thisWord,
                                  countsSoFarB_.count_[l],
                                  countsThisRangeB.count_[l],
                                  subset_ );
                } // ~if
            } // ~for l
            if ( hasChild == false )
            {
                //  if no children, print word itself
                //   cout << "GOLD " << thisRangeB.word_ << " " << thisRangeB.num_ << endl;
                //  numSingletonRanges_++;
            } // ~if

            countsSoFarB_ += countsThisRangeB;
            currentPosB_ += thisRangeB.num_;

            //      numRanges_++;
        } // ~while

        if ( notAtLastA == false )
        {
            assert ( notAtLastB == false );
            break;
        } // ~if
        else
        {
            assert( ( thisRangeA.pos_ & matchFlag ) != 0 );
            assert( ( thisRangeB.pos_ & matchFlag ) != 0 );

#ifdef DEBUG
            cout << "RangeA: " << i << " "
                 << thisRangeA.word_ << " " << thisRangeA.pos_ << " " << thisRangeA.num_
                 << " -- " << currentPosA_ << " = ";
            cout << "RangeB: " << i << " "
                 << thisRangeB.word_ << " " << thisRangeB.pos_ << " " << thisRangeB.num_
                 << " -- " << currentPosB_ << endl;
            cout << ( thisRangeA.pos_ & matchFlag ) << " " << ( thisRangeB.pos_ & matchFlag ) << endl;
#endif


            skipIfNecessary( thisRangeA, currentPosA_, ( *inBwtA_ ), countsSoFarA_ );
            // count children
            countsThisRangeA.clear();
            inBwtA_->readAndCount( countsThisRangeA, thisRangeA.num_ );
#ifdef DEBUG
            countsThisRangeA.print();
#endif

            skipIfNecessary( thisRangeB, currentPosB_, ( *inBwtB_ ), countsSoFarB_ );
            // count children
            countsThisRangeB.clear();
            inBwtB_->readAndCount( countsThisRangeB, thisRangeB.num_ );
#ifdef DEBUG
            countsThisRangeB.print();
#endif

            assert( thisRangeA.word_.size() == thisRangeB.word_.size() );
            intervalHandler_.foundInBoth
            ( i,
              countsThisRangeA, countsThisRangeB,
              thisRangeA, thisRangeB,
              propagateIntervalA_, propagateIntervalB_ );


#ifdef PROPAGATE_PREFIX
            hasChild = false;
#endif
            for ( int l( 1 ); l < alphabetSize; l++ )
            {
                if ( ( propagateIntervalA_[l] == true )
                     || ( propagateIntervalB_[l] == true ) )
                {
#ifdef PROPAGATE_PREFIX
                    if ( hasChild == false )
                    {
                        // assert(thisWord.size()==thisRangeA.word_.size()+1);
                        thisWord.replace( 1, thisRangeA.word_.size(), thisRangeA.word_ );
                        hasChild = true;
                    }
                    thisWord[0] = alphabet[l];
#endif
                    //  thisWord+=thisRangeA.word_;
                    LetterCountType thisFlag
                    ( ( ( propagateIntervalA_[l] == true )
                        && ( propagateIntervalB_[l] == true ) ) ? matchFlag : 0 );


                    if ( propagateIntervalA_[l] == true )
                        rA_.addRange( l, i, thisWord,
                                      ( countsSoFarA_.count_[l]
                                        | thisFlag ),
                                      countsThisRangeA.count_[l],
                                      subset_ );
                    // if (countsThisRangeB.count_[l]>0)
                    if ( propagateIntervalB_[l] == true )
                        rB_.addRange( l, i, thisWord,
                                      ( countsSoFarB_.count_[l]
                                        | thisFlag ),
                                      countsThisRangeB.count_[l],
                                      subset_ );
                } // ~if
            } // ~for

            countsSoFarA_ += countsThisRangeA;
            currentPosA_ += thisRangeA.num_;
            countsSoFarB_ += countsThisRangeB;
            currentPosB_ += thisRangeB.num_;

            numRanges_++;


        }


    } // ~while


} // ~BackTracker::operator()

#endif

