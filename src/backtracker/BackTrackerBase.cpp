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

#include "BackTrackerBase.hh"

#include "EndPosFile.hh"
#include "IntervalHandlerBase.hh"

#include "libzoo/util/Logger.hh"

#include <algorithm>

using namespace std;


BackTrackerBase::BackTrackerBase( const string &subset, const int cycle, const bool noComparisonSkip, const bool propagateSequence )
    : subset_( subset )
    , cycle_( cycle )
    , noComparisonSkip_( noComparisonSkip )
    , numRanges_( 0 )
    , numSingletonRanges_( 0 )
    , propagateSequence_( propagateSequence )
{}

void BackTrackerBase::skipIfNecessary( const Range &thisRange,
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
} // ~BackTrackerBase::skipIfNecessary



void BackTrackerBase::processSingletons(
    const int pileNum
    , bool &notAtLast
    , RangeStoreExternal &rA_
    , Range &thisRange
    , LetterNumber &currentPos
    , BwtReaderBase *inBwt
    , LetterCount &countsSoFar
    , LetterCount &countsThisRange
    , IntervalHandlerBase &intervalHandler
    , AlphabetFlag &propagateInterval
    , string &thisWord
    , const bool doesPropagateToEnd
    , IntervalHandler_FoundCallbackPtr foundCallbackPtr
    , EndPosFile &endPosFile
    , int sampleId
)
{
    vector<Range> rAPile0;

    while ( notAtLast )
    {
        notAtLast = rA_.getRange( thisRange );
        Logger_if( LOG_FOR_DEBUGGING )
        {
            Logger::out() << "notAtLast=" << notAtLast << endl;
            Logger::out() << "thisRange.pos_=" << thisRange.pos_ << endl;
        }

        Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
        {
            if ( notAtLast )
            {
                Logger::out() << "RangeA: " << pileNum << " "
                              << thisRange.word_ << " "
                              << thisRange.pos_ << " " << thisRange.num_
                              << " -- " << currentPos
                              << " < " << thisRange.word_
                              << " match=" << ( thisRange.pos_ & matchFlag )
                              << endl;
            }
            else
            {
                Logger::out() << "RangeA: last reached"  << endl;
            }
        }
        if ( ( notAtLast == false ) || ( ( thisRange.pos_ & matchFlag ) != 0 ) ) break;

        char *bwtSubstring = NULL;
        prepareCallbackArgs(
            thisRange
            , currentPos
            , inBwt
            , countsSoFar
            , countsThisRange
            , intervalHandler
            , bwtSubstring
        );

        if ( sampleId == 1 )
            intervalHandler.foundInAOnly
            (
                pileNum,
                countsSoFar,
                countsThisRange,
                bwtSubstring,
                thisRange,
                propagateInterval,
                cycle_
            );
        else
            intervalHandler.foundInBOnly
            (
                pileNum,
                countsSoFar,
                countsThisRange,
                bwtSubstring,
                thisRange,
                propagateInterval,
                cycle_
            );
        /*
                ( intervalHandler.*foundCallbackPtr )
                (
                    pileNum,
                    countsSoFar,
                    countsThisRange,
                    bwtSubstring,
                    thisRange,
                    propagateInterval,
                    cycle_
                );
        */

        // Sequence numbers extraction
        if ( doesPropagateToEnd )
        {
            if ( countsThisRange.count_[0] > 0 && thisRange.isBkptExtension_ )
            {
                #pragma omp critical (IO)
                cout << "READNUM"
                     << "(" << sampleId  << ")"
                     << " " << pileNum
                     << " " << ( thisRange.pos_ + 1 )
                     << " " << thisRange.num_
                     << " " << countsThisRange.count_[0]
                     << endl;
            }
            for ( int l( 0 ); l < alphabetSize; l++ )
                propagateInterval[l] = ( countsThisRange.count_[l] > 0 );
        }

        // add ranges for any children
        bool hasChild = false;
        for ( int l( 0 ); l < alphabetSize; l++ )
        {
            //       if (countsThisRange.count_[l]>=minOcc)
            if ( propagateInterval[l] == true )
            {
                updatePropagatedSuffixWord( hasChild, thisRange, thisWord, l );

#ifdef OLD_BACKTRACKER_COMPAT
                Range newRange( thisWord,
                                countsSoFar.count_[l],
                                countsThisRange.count_[l],
                                thisRange.isBkptExtension_ );
#else
                Range &newRange = intervalHandler.getSubIntervalRange (
                                      thisWord,
                                      countsSoFar.count_[l],
                                      countsThisRange.count_[l],
                                      thisRange.isBkptExtension_,
                                      thisRange,
                                      l
                                  );
#endif
                if ( noComparisonSkip_ ||
                     !rA_.isRangeKnown( newRange, l, pileNum, subset_, cycle_ ) )
                {
                    if ( l != 0 )
                    {
                        rA_.addRange( newRange, l, pileNum, subset_, cycle_ );
                    }
                    else
                    {
                        Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
                        {
                            Logger::out() << "Adding range item to $ pile: ";
                            newRange.prettyPrint( Logger::out() );
                            Logger::out() << ", pileNum=" << pileNum << ", subset=" << subset_ << ", cycle=" << cycle_ << endl;
                        }
                        LetterNumber rangeStartPos = newRange.pos_;
                        LetterNumber rangeLength = newRange.num_;
                        newRange.num_ = 1;
                        for ( LetterNumber i = 0; i < rangeLength; ++i )
                        {
                            newRange.pos_ = rangeStartPos + i;
                            SequenceNumber newPos = endPosFile.convertDollarNumToSequenceNum( newRange.pos_ );
                            newRange.pos_ = newPos;
                            //                        rAPile0.push_back( newRange );

                            rA_.addOutOfOrderRange( newRange, l, pileNum, subset_, cycle_ );
                        }
                    }
                }

#ifdef OLD_BACKTRACKER_COMPAT
#else
                //                delete newRange;
#endif
            } // ~if

            /*
                        std::sort( rAPile0.begin(), rAPile0.end(), compareRangeByPos );
                        for( Range &r: rAPile0 )
                            rA_.addRange( r, l, pileNum, subset_, cycle_ );
            */
            rAPile0.clear();

        } // ~for l

        if ( hasChild == false )
        {
#ifdef OLD
            //  if no children, print word itself
            cout << "GOLD ";
            cout << thisRange.word_;
            cout << " " << thisRange.num_ << endl;
#endif
            numSingletonRanges_++;
        } // ~if

        countsSoFar += countsThisRange;
        currentPos += thisRange.num_;

        numRanges_++;
    } // ~while notAtLast
} /// END OF BLOCK1




void BackTrackerBase::prepareCallbackArgs(
    Range &thisRange
    , LetterNumber &currentPos
    , BwtReaderBase *inBwt
    , LetterCount &countsSoFar
    , LetterCount &countsThisRange
    , IntervalHandlerBase &intervalHandler
    , char *&bwtSubstring
)
{
    skipIfNecessary( thisRange, currentPos, ( *inBwt ), countsSoFar );

    if ( intervalHandler.needSubstring )
    {
        if ( thisRange.num_ > bwtSubstringStore_.size() )
            bwtSubstringStore_.resize( thisRange.num_ );
        bwtSubstring = &bwtSubstringStore_[0];

        LetterNumber charsRead = ( *inBwt )( bwtSubstring, thisRange.num_ );
        assert( charsRead == thisRange.num_ );
        countsThisRange.clear();
        countsThisRange.countString( bwtSubstring, thisRange.num_ );
    }
    else
    {
        // count children
        countsThisRange.clear();
        LetterNumber charsRead = inBwt->readAndCount( countsThisRange, thisRange.num_ );
        assert( charsRead == thisRange.num_ );
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
        {
            Logger::out() << countsThisRange << endl;
        }

    }
} // END_OF_BLOCK2



void BackTrackerBase::updatePropagatedSuffixWord( bool &hasChild, const Range &thisRange, string &thisWord, const AlphabetSymbol l )
{
    if ( propagateSequence_ )
    {
        if ( hasChild == false )
        {
            assert( thisWord.size() == thisRange.word_.size() + 1 );
            thisWord.replace( 1, thisRange.word_.size(), thisRange.word_ );
        } // ~if
        thisWord[0] = alphabet[l];
    }
    hasChild = true;
} // END_OF_BLOCK3
