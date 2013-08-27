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

#include "BackTracker.hh"

#include "libzoo/util/Logger.hh"

using namespace std;


BackTracker::BackTracker( BwtReaderBase *inBwtA, BwtReaderBase *inBwtB,
                          LetterNumber &currentPosA, LetterNumber &currentPosB,
                          RangeStoreExternal &rA, RangeStoreExternal &rB,
                          LetterCount &countsSoFarA, LetterCount &countsSoFarB,
                          int minOcc, const int maxLength, const string &subset, const int cycle,
                          const bool doesPropagateBkptToSeqNumInSetA,
                          const bool doesPropagateBkptToSeqNumInSetB,
                          const bool noComparisonSkip ):
    inBwtA_( inBwtA ), inBwtB_( inBwtB ),
    currentPosA_( currentPosA ), currentPosB_( currentPosB ),
    rA_( rA ), rB_( rB ), countsSoFarA_( countsSoFarA ), countsSoFarB_( countsSoFarB ),
    minOcc_( minOcc ), maxLength_( maxLength ), subset_( subset ), cycle_( cycle ), numRanges_( 0 ), numSingletonRanges_( 0 )
    , doesPropagateBkptToSeqNumInSetA_( doesPropagateBkptToSeqNumInSetA )
    , doesPropagateBkptToSeqNumInSetB_( doesPropagateBkptToSeqNumInSetB )
    , noComparisonSkip_( noComparisonSkip )
{}

void BackTracker::skipIfNecessary( const Range &thisRange,
                                   LetterNumber &currentPos,
                                   BwtReaderBase &inBwt,
                                   LetterCount &countsSoFar )
{
    Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Want 2 skip " << thisRange.pos_ - currentPos << ": " << currentPos << " to " << thisRange.pos_ << endl;
    if ( ( thisRange.pos_ & matchMask ) > currentPos )
    {
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Skipping " << thisRange.pos_ - currentPos << ": " << currentPos << " to " << thisRange.pos_ << endl;
        inBwt.readAndCount( countsSoFar, ( thisRange.pos_ & matchMask ) - currentPos );
        currentPos = ( thisRange.pos_ & matchMask );
    } // ~if
    if ( thisRange.pos_ < currentPos )
    {
        cerr << "thisRange is too low. Should be > " << currentPos << "." << endl;
    }
} // ~BackTracker::skipIfNecessary

void BackTracker::operator() ( int i, string &thisWord, IntervalHandlerBase &intervalHandler_  )
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

            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "RangeA: " << i << " "
#ifdef PROPAGATE_PREFIX
                    << thisRangeA.word_ << " "
#endif
                    << thisRangeA.pos_ << " " << thisRangeA.num_
                    << " -- " << currentPosA_
#ifdef PROPAGATE_PREFIX
                    << " < " << thisRangeB.word_
#endif
                    << endl;

            skipIfNecessary( thisRangeA, currentPosA_, ( *inBwtA_ ), countsSoFarA_ );
            // count children
            countsThisRangeA.clear();
            inBwtA_->readAndCount( countsThisRangeA, thisRangeA.num_ );
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                countsThisRangeA.print();
            }

            intervalHandler_.foundInAOnly
            ( i,
              countsSoFarA_,
              countsThisRangeA,
              thisRangeA,
              propagateIntervalA_ );

            // Sequence numbers extraction
            if ( doesPropagateBkptToSeqNumInSetA_ )
            {
                if ( countsThisRangeA.count_[0] > 0 )
                {
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
#ifdef PROPAGATE_PREFIX
                    if ( hasChild == false )
                    {
                        // assert(thisWord.size()==thisRangeA.word_.size()+1);
                        thisWord.replace( 1, thisRangeA.word_.size(), thisRangeA.word_ );
                    } // ~if
                    thisWord[0] = alphabet[l];
#endif
                    hasChild = true;

                    if ( noComparisonSkip_ ||
                         rA_.isRangeKnown( l, i, thisWord,
                                           countsSoFarA_.count_[l],
                                           countsThisRangeA.count_[l],
                                           thisRangeA.isBkptExtension_,
                                           subset_, cycle_ ) )
                        rA_.addRange( l, i, thisWord,
                                      countsSoFarA_.count_[l],
                                      countsThisRangeA.count_[l],
                                      thisRangeA.isBkptExtension_,
                                      subset_, cycle_ );
                } // ~if
            } // ~for l

            if ( hasChild == false )
            {
#ifdef OLD
                //  if no children, print word itself
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

            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "RangeB: " << i << " "
#ifdef PROPAGATE_PREFIX
                    << thisRangeB.word_ << " "
#endif
                    << thisRangeB.pos_ << " " << thisRangeB.num_
                    << " -- " << currentPosB_
#ifdef PROPAGATE_PREFIX
                    << " < " << thisRangeB.word_
#endif
                    << endl;

            skipIfNecessary( thisRangeB, currentPosB_, ( *inBwtB_ ), countsSoFarB_ );
            // count children
            countsThisRangeB.clear();
            uint charsRead = inBwtB_->readAndCount( countsThisRangeB, thisRangeB.num_ );
            assert( charsRead == thisRangeB.num_ );
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                countsThisRangeB.print();
            }

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

                    if ( noComparisonSkip_ ||
                         rB_.isRangeKnown( l, i, thisWord,
                                           countsSoFarB_.count_[l],
                                           countsThisRangeB.count_[l],
                                           thisRangeB.isBkptExtension_,
                                           subset_, cycle_ ) )
                        rB_.addRange( l, i, thisWord,
                                      countsSoFarB_.count_[l],
                                      countsThisRangeB.count_[l],
                                      thisRangeB.isBkptExtension_,
                                      subset_, cycle_ );
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
#ifdef PROPAGATE_PREFIX
                    << thisRangeA.word_ << " "
#endif
                    << thisRangeA.pos_ << " " << thisRangeA.num_
                    << " -- " << currentPosA_ << " = ";
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "RangeB: " << i << " "
#ifdef PROPAGATE_PREFIX
                    << thisRangeB.word_ << " "
#endif
                    << thisRangeB.pos_ << " " << thisRangeB.num_
                    << " -- " << currentPosB_ << endl;
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << ( thisRangeA.pos_ & matchFlag ) << " " << ( thisRangeB.pos_ & matchFlag ) << endl;


            skipIfNecessary( thisRangeA, currentPosA_, ( *inBwtA_ ), countsSoFarA_ );
            // count children
            countsThisRangeA.clear();
            inBwtA_->readAndCount( countsThisRangeA, thisRangeA.num_ );
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                countsThisRangeA.print();
            }

            skipIfNecessary( thisRangeB, currentPosB_, ( *inBwtB_ ), countsSoFarB_ );
            // count children
            countsThisRangeB.clear();
            inBwtB_->readAndCount( countsThisRangeB, thisRangeB.num_ );
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
            {
                countsThisRangeB.print();
            }

#ifdef PROPAGATE_PREFIX
            assert( thisRangeA.word_.size() == thisRangeB.word_.size() );
#endif
            bool isBreakpointDetected = false;
            intervalHandler_.foundInBoth
            ( i,
              countsThisRangeA, countsThisRangeB,
              thisRangeA, thisRangeB,
              propagateIntervalA_, propagateIntervalB_,
              isBreakpointDetected );
            if ( isBreakpointDetected )
            {
                thisRangeA.isBkptExtension_ = true;
                thisRangeB.isBkptExtension_ = true;
            }


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
                    LetterNumber thisFlag
                    ( ( ( propagateIntervalA_[l] == true )
                        && ( propagateIntervalB_[l] == true ) ) ? matchFlag : 0 );

                    bool isRangeKnownA = noComparisonSkip_ ||
                                         rA_.isRangeKnown( l, i, thisWord,
                                                           ( countsSoFarA_.count_[l]
                                                                   | thisFlag ),
                                                           countsThisRangeA.count_[l],
                                                           thisRangeA.isBkptExtension_,
                                                           subset_, cycle_ );
                    bool isRangeKnownB = noComparisonSkip_ ||
                                         rB_.isRangeKnown( l, i, thisWord,
                                                           ( countsSoFarB_.count_[l]
                                                                   | thisFlag ),
                                                           countsThisRangeB.count_[l],
                                                           thisRangeB.isBkptExtension_,
                                                           subset_, cycle_ );
                    if ( thisFlag )
                    {
                        isRangeKnownA = isRangeKnownB = isRangeKnownA || isRangeKnownB;
                    }

                    if ( propagateIntervalA_[l] == true && isRangeKnownA )
                        rA_.addRange( l, i, thisWord,
                                      ( countsSoFarA_.count_[l]
                                        | thisFlag ),
                                      countsThisRangeA.count_[l],
                                      thisRangeA.isBkptExtension_,
                                      subset_, cycle_ );
                    // if (countsThisRangeB.count_[l]>0)
                    if ( propagateIntervalB_[l] == true && isRangeKnownB )
                        rB_.addRange( l, i, thisWord,
                                      ( countsSoFarB_.count_[l]
                                        | thisFlag ),
                                      countsThisRangeB.count_[l],
                                      thisRangeB.isBkptExtension_,
                                      subset_, cycle_ );
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
