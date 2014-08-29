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

#include "TwoBwtBackTracker.hh"

#include "libzoo/util/Logger.hh"

using namespace std;


TwoBwtBackTracker::TwoBwtBackTracker( BwtReaderBase *inBwtA, BwtReaderBase *inBwtB,
                                      LetterNumber &currentPosA, LetterNumber &currentPosB,
                                      RangeStoreExternal &rA, RangeStoreExternal &rB,
                                      LetterCount &countsSoFarA, LetterCount &countsSoFarB,
                                      int minOcc, const int maxLength, const string &subset, const int cycle,
                                      const bool doesPropagateBkptToSeqNumInSetA,
                                      const bool doesPropagateBkptToSeqNumInSetB,
                                      const bool noComparisonSkip,
                                      const bool propagateSequence )
    : BackTrackerBase( subset, cycle, noComparisonSkip, propagateSequence )
    , inBwtA_( inBwtA ), inBwtB_( inBwtB )
    , currentPosA_( currentPosA ), currentPosB_( currentPosB )
    , rA_( rA ), rB_( rB ), countsSoFarA_( countsSoFarA ), countsSoFarB_( countsSoFarB )
    , minOcc_( minOcc ), maxLength_( maxLength ), subset_( subset ), cycle_( cycle )
    , numNotSkippedA_( 0 )
    , numNotSkippedB_( 0 )
    , numSkippedA_( 0 )
    , numSkippedB_( 0 )
    , doesPropagateBkptToSeqNumInSetA_( doesPropagateBkptToSeqNumInSetA )
    , doesPropagateBkptToSeqNumInSetB_( doesPropagateBkptToSeqNumInSetB )
{
    for ( int l( 0 ); l < alphabetSize; ++l )
        propagateIntervalA_[l] = propagateIntervalB_[l] = false;
}

/*
void TwoBwtBackTracker::process (
    int pileNum
    , string &thisWord
    , IntervalHandlerBase &intervalHandler
)
{
    LetterCount countsThisRangeA, countsThisRangeB;
    Range thisRangeA, thisRangeB;
    //string thisWord;
    bool notAtLastA( true ), notAtLastB( true );
    bool hasChild;

    while ( 1 )
    {
        processSingletons(
            pileNum
            , notAtLastA
            , rA_
            , thisRangeA
            , currentPosA_
            , inBwtA_
            , countsSoFarA_
            , countsThisRangeA
            , intervalHandler
            , propagateIntervalA_
            , thisWord
            , doesPropagateBkptToSeqNumInSetA_
            , NULL //( IntervalHandler_FoundCallbackPtr )( &IntervalHandlerBase::foundInAOnly )
            , 1
        );

        processSingletons(
            pileNum
            , notAtLastB
            , rB_
            , thisRangeB
            , currentPosB_
            , inBwtB_
            , countsSoFarB_
            , countsThisRangeB
            , intervalHandler
            , propagateIntervalB_
            , thisWord
            , doesPropagateBkptToSeqNumInSetB_
            , NULL //( IntervalHandler_FoundCallbackPtr )( &IntervalHandlerBase::foundInBOnly )
            , 2
        );

        if ( notAtLastA == false )
        {
            assert ( notAtLastB == false );
            break;
        } // ~if
        else
        {
            assert( ( thisRangeA.pos_ & matchFlag ) != 0 );
            assert( ( thisRangeB.pos_ & matchFlag ) != 0 );

            char *bwtSubstringA;
            prepareCallbackArgs(
                thisRangeA
                , currentPosA_
                , inBwtA_
                , countsSoFarA_
                , countsThisRangeA
                , intervalHandler
                , bwtSubstringA
            );

            char *bwtSubstringB;
            prepareCallbackArgs(
                thisRangeB
                , currentPosB_
                , inBwtB_
                , countsSoFarB_
                , countsThisRangeB
                , intervalHandler
                , bwtSubstringB
            );


#ifdef PROPAGATE_SEQUENCE
            assert( thisRangeA.word_.size() == thisRangeB.word_.size() );
#endif
            bool isBreakpointDetected = false;
            intervalHandler.foundInBoth
            ( pileNum,
              countsThisRangeA,
              countsThisRangeB,
              thisRangeA,
              thisRangeB,
              propagateIntervalA_,
              propagateIntervalB_,
              isBreakpointDetected,
              cycle_
            );
            if ( isBreakpointDetected )
            {
                thisRangeA.isBkptExtension_ = true;
                thisRangeB.isBkptExtension_ = true;
            }


#ifdef PROPAGATE_SEQUENCE
            hasChild = false;
#endif
            for ( AlphabetSymbol l( 1 ); l < alphabetSize; l++ )
            {
                if ( ( propagateIntervalA_[l] == true )
                     || ( propagateIntervalB_[l] == true ) )
                {
                    updatePropagatedSuffixWord( hasChild, thisRangeA, thisWord, l );
                    //  thisWord+=thisRangeA.word_;
                    LetterNumber thisFlag
                    ( ( ( propagateIntervalA_[l] == true )
                        && ( propagateIntervalB_[l] == true ) ) ? matchFlag : 0 );

                    Range newRangeA( thisWord,
                                     ( countsSoFarA_.count_[l]
                                       | thisFlag ),
                                     countsThisRangeA.count_[l],
                                     thisRangeA.isBkptExtension_ );
                    Range newRangeB( thisWord,
                                     ( countsSoFarB_.count_[l]
                                       | thisFlag ),
                                     countsThisRangeB.count_[l],
                                     thisRangeB.isBkptExtension_ );

                    bool doAddRangeA = noComparisonSkip_ || !rA_.isRangeKnown( newRangeA, l, pileNum, subset_, cycle_ );
                    bool doAddRangeB = noComparisonSkip_ || !rB_.isRangeKnown( newRangeB, l, pileNum, subset_, cycle_ );
                    if ( thisFlag )
                    {
                        doAddRangeA = doAddRangeB = doAddRangeA || doAddRangeB;
                    }

                    if ( propagateIntervalA_[l] == true )
                    {
                        if ( doAddRangeA )
                        {
                            rA_.addRange( newRangeA, l, pileNum, subset_, cycle_ );
                            ++numNotSkippedA_;
                        }
                        else
                        {
                            ++numSkippedA_;
                        }
                    }
                    // if (countsThisRangeB.count_[l]>0)
                    if ( propagateIntervalB_[l] == true )
                    {
                        if ( doAddRangeB )
                        {
                            rB_.addRange( newRangeB, l, pileNum, subset_, cycle_ );
                            ++numNotSkippedB_;
                        }
                        else
                        {
                            ++numSkippedB_;
                        }
                    }
                } // ~if
            } // ~for

            countsSoFarA_ += countsThisRangeA;
            currentPosA_ += thisRangeA.num_;
            countsSoFarB_ += countsThisRangeB;
            currentPosB_ += thisRangeB.num_;

            numRanges_++;


        }


    } // ~while


}
*/




void TwoBwtBackTracker::process ( int i, string &thisWord, IntervalHandlerBase &intervalHandler_  )
{
    LetterCount countsThisRangeA, countsThisRangeB;
    Range thisRangeA, thisRangeB;
    //  string thisWord;
    bool notAtLastA( true ), notAtLastB( true );
    bool hasChild = false;

    while ( 1 )
    {
        while ( notAtLastA )
        {
            notAtLastA = rA_.getRange( thisRangeA );
            if ( ( notAtLastA == false ) || ( ( thisRangeA.pos_ & matchFlag ) != 0 ) ) break;

            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "RangeA: " << i << " "
                    << thisRangeA.word_ << " "
                    << thisRangeA.pos_ << " " << thisRangeA.num_
                    << " -- " << currentPosA_
                    << " < " << thisRangeB.word_
                    << endl;

            skipIfNecessary( thisRangeA, currentPosA_, ( *inBwtA_ ), countsSoFarA_ );
            // count children
            countsThisRangeA.clear();
            inBwtA_->readAndCount( countsThisRangeA, thisRangeA.num_ );
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                Logger::out() << countsThisRangeA << endl;
            }

            intervalHandler_.foundInAOnly
            ( i,
              countsSoFarA_,
              countsThisRangeA,
              NULL,
              thisRangeA,
              propagateIntervalA_ ,
              cycle_
            );

            // Sequence numbers extraction
            if ( doesPropagateBkptToSeqNumInSetA_ )
            {
                if ( countsThisRangeA.count_[0] > 0 )
                {
                    #pragma omp critical (IO)
                    cout << "READNUM"
                         << "(" << thisRangeA.isBkptExtension_ << ")"
                         << " " << i
                         << " " << ( thisRangeA.pos_ + 1 )
                         << " " << thisRangeA.num_
                         << " " << countsThisRangeA.count_[0]
                         << endl;
                }
                for ( int l( 1 ); l < alphabetSize; l++ )
                    propagateIntervalA_[l] = ( countsThisRangeA.count_[l] > 0 );
            }

            // add ranges for any children
            hasChild = false;
            for ( int l( 1 ); l < alphabetSize; l++ )
            {
                //       if (countsThisRangeA.count_[l]>=minOcc)
                if ( propagateIntervalA_[l] == true )
                {
                    if ( propagateSequence_ )
                    {
                        if ( hasChild == false )
                        {
                            assert( thisWord.size() == thisRangeA.word_.size() + 1 );
                            thisWord.replace( 1, thisRangeA.word_.size(), thisRangeA.word_ );
                        } // ~if
                        thisWord[0] = alphabet[l];
                    }
                    hasChild = true;

                    Range newRange( thisWord,
                                    countsSoFarA_.count_[l],
                                    countsThisRangeA.count_[l],
                                    thisRangeA.isBkptExtension_ );
                    if ( noComparisonSkip_ ||
                         !rA_.isRangeKnown( newRange, l, i, subset_, cycle_ ) )
                        rA_.addRange( newRange, l, i, subset_, cycle_ );
                } // ~if
            } // ~for l

            if ( hasChild == false )
            {
#ifdef OLD
                //  if no children, print word itself
                cout << "GOLD ";
                cout << thisRangeA.word_;
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

            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "RangeB: " << i << " "
                    << thisRangeB.word_ << " "
                    << thisRangeB.pos_ << " " << thisRangeB.num_
                    << " -- " << currentPosB_
                    << " < " << thisRangeB.word_
                    << endl;

            skipIfNecessary( thisRangeB, currentPosB_, ( *inBwtB_ ), countsSoFarB_ );
            // count children
            countsThisRangeB.clear();
            uint charsRead = inBwtB_->readAndCount( countsThisRangeB, thisRangeB.num_ );
            assert( charsRead == thisRangeB.num_ );
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                Logger::out() << countsThisRangeB << endl;
            }

            //#ifdef TBD
            intervalHandler_.foundInBOnly
            ( i,
              countsSoFarB_,
              countsThisRangeB,
              NULL,
              thisRangeB,
              propagateIntervalB_,
              cycle_
            );
            //#endif

            // add ranges for any children
            if ( propagateSequence_ )
                hasChild = false;
            for ( int l( 1 ); l < alphabetSize; l++ )
            {
                //       if (countsThisRangeB.count_[l]>=minOcc)
                if ( propagateIntervalB_[l] == true )
                {
                    if ( propagateSequence_ )
                    {
                        if ( hasChild == false )
                        {
                            assert( thisWord.size() == thisRangeB.word_.size() + 1 );
                            thisWord.replace( 1, thisRangeB.word_.size(), thisRangeB.word_ );
                            hasChild = true;
                        } // ~if
                        thisWord[0] = alphabet[l];
                    }

                    Range newRange( thisWord,
                                    countsSoFarB_.count_[l],
                                    countsThisRangeB.count_[l],
                                    thisRangeB.isBkptExtension_ );
                    if ( noComparisonSkip_ ||
                         !rB_.isRangeKnown( newRange, l, i, subset_, cycle_ ) )
                        rB_.addRange( newRange, l, i, subset_, cycle_ );
                } // ~if
            } // ~for l
            /*
                        if ( hasChild == false )
                        {
                            //  if no children, print word itself
                            //   cout << "GOLD " << thisRangeB.word_ << " " << thisRangeB.num_ << endl;
                            //  numSingletonRanges_++;
                        } // ~if
            */
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

            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "RangeA: " << i << " "
                    << thisRangeA.word_ << " "
                    << thisRangeA.pos_ << " " << thisRangeA.num_
                    << " -- " << currentPosA_ << " = ";
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "RangeB: " << i << " "
                    << thisRangeB.word_ << " "
                    << thisRangeB.pos_ << " " << thisRangeB.num_
                    << " -- " << currentPosB_ << endl;
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << ( thisRangeA.pos_ & matchFlag ) << " " << ( thisRangeB.pos_ & matchFlag ) << endl;


            skipIfNecessary( thisRangeA, currentPosA_, ( *inBwtA_ ), countsSoFarA_ );
            // count children
            countsThisRangeA.clear();
            LetterNumber count = inBwtA_->readAndCount( countsThisRangeA, thisRangeA.num_ );
            assert( count == thisRangeA.num_ );
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                Logger::out() << countsThisRangeA << endl;
            }

            skipIfNecessary( thisRangeB, currentPosB_, ( *inBwtB_ ), countsSoFarB_ );
            // count children
            countsThisRangeB.clear();
            count = inBwtB_->readAndCount( countsThisRangeB, thisRangeB.num_ );
            assert( count == thisRangeB.num_ );
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                Logger::out() << countsThisRangeB << endl;
            }

            if ( propagateSequence_ )
                assert( thisRangeA.word_.size() == thisRangeB.word_.size() );

            bool isBreakpointDetected = false;
            intervalHandler_.foundInBoth
            ( i,
              countsThisRangeA, countsThisRangeB,
              thisRangeA, thisRangeB,
              propagateIntervalA_, propagateIntervalB_,
              isBreakpointDetected, cycle_ );
            if ( isBreakpointDetected )
            {
                thisRangeA.isBkptExtension_ = true;
                thisRangeB.isBkptExtension_ = true;
            }

            if ( propagateSequence_ )
                hasChild = false;

            for ( int l( 1 ); l < alphabetSize; l++ )
            {
                if ( ( propagateIntervalA_[l] == true )
                     || ( propagateIntervalB_[l] == true ) )
                {
                    if ( propagateSequence_ )
                    {
                        if ( hasChild == false )
                        {
                            assert( thisWord.size() == thisRangeA.word_.size() + 1 );
                            thisWord.replace( 1, thisRangeA.word_.size(), thisRangeA.word_ );
                            hasChild = true;
                        }
                        thisWord[0] = alphabet[l];
                    }

                    //  thisWord+=thisRangeA.word_;
                    LetterNumber thisFlag
                    ( ( ( propagateIntervalA_[l] == true )
                        && ( propagateIntervalB_[l] == true ) ) ? matchFlag : 0 );

                    Range newRangeA( thisWord,
                                     ( countsSoFarA_.count_[l]
                                       | thisFlag ),
                                     countsThisRangeA.count_[l],
                                     thisRangeA.isBkptExtension_ );
                    Range newRangeB( thisWord,
                                     ( countsSoFarB_.count_[l]
                                       | thisFlag ),
                                     countsThisRangeB.count_[l],
                                     thisRangeB.isBkptExtension_ );

                    bool doAddRangeA = noComparisonSkip_ || !rA_.isRangeKnown( newRangeA, l, i, subset_, cycle_ );
                    bool doAddRangeB = noComparisonSkip_ || !rB_.isRangeKnown( newRangeB, l, i, subset_, cycle_ );
                    if ( thisFlag )
                    {
                        doAddRangeA = doAddRangeB = doAddRangeA || doAddRangeB;
                    }

                    if ( propagateIntervalA_[l] == true )
                    {
                        if ( doAddRangeA )
                        {
                            rA_.addRange( newRangeA, l, i, subset_, cycle_ );
                            ++numNotSkippedA_;
                        }
                        else
                        {
                            ++numSkippedA_;
                        }
                    }
                    // if (countsThisRangeB.count_[l]>0)
                    if ( propagateIntervalB_[l] == true )
                    {
                        if ( doAddRangeB )
                        {
                            rB_.addRange( newRangeB, l, i, subset_, cycle_ );
                            ++numNotSkippedB_;
                        }
                        else
                        {
                            ++numSkippedB_;
                        }
                    }
                } // ~if
            } // ~for

            countsSoFarA_ += countsThisRangeA;
            currentPosA_ += thisRangeA.num_;
            countsSoFarB_ += countsThisRangeB;
            currentPosB_ += thisRangeB.num_;

            numRanges_++;


        }


    } // ~while
}
