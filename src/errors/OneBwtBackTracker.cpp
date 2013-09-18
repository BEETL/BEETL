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

#include "OneBwtBackTracker.hh"

#include "BwtCorrectorIntervalHandler.hh"
#include "ErrorCorrectionRange.hh"
#include "libzoo/util/Logger.hh"

using namespace std;

OneBwtBackTracker::OneBwtBackTracker(
    BwtReaderBase *inBwt,
    LetterNumber &currentPos,
    RangeStoreExternal &r,
    LetterCount &countsSoFar,
    int minOcc,
    const int maxLength,
    const string &subset,
    const int cycle,
    const bool doesPropagateBkptToSeqNumInSet,
    const bool noComparisonSkip
):
    inBwt_( inBwt ),
    currentPos_( currentPos ),
    r_( r ),
    countsSoFar_( countsSoFar ),
    minOcc_( minOcc ),
    maxLength_( maxLength ),
    subset_( subset ),
    cycle_( cycle ),
    numRanges_( 0 ),
    numSingletonRanges_( 0 ),
    doesPropagateBkptToSeqNumInSet_( doesPropagateBkptToSeqNumInSet ),
    noComparisonSkip_( noComparisonSkip )
{}

void OneBwtBackTracker::skipIfNecessary(
    const Range &thisRange,
    LetterNumber &currentPos,
    BwtReaderBase &inBwt,
    LetterCount &countsSoFar
)
{
    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "Want 2 skip " << thisRange.pos_ - currentPos << ": " << currentPos << " to " << thisRange.pos_ << endl;
    if ( ( thisRange.pos_ & matchMask ) > currentPos )
    {
        Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "Skipping " << thisRange.pos_ - currentPos << ": " << currentPos << " to " << thisRange.pos_ << endl;
        inBwt.readAndCount( countsSoFar, ( thisRange.pos_ & matchMask ) - currentPos );
        currentPos = ( thisRange.pos_ & matchMask );
    } // ~if
    if ( thisRange.pos_ < currentPos )
    {
        cerr << "thisRange is too low. Should be > " << currentPos << "." << endl;
    }
} // ~OneBwtBackTracker::skipIfNecessary

void OneBwtBackTracker::operator() (
    int pileNum,
    string &thisWord,
    IntervalHandlerBase &intervalHandler_
)
{
    LetterCount countsThisRange;
    ErrorCorrectionRange thisRange;
    bool notAtLast( true );
    bool hasChild;

    while ( notAtLast )
    {
        notAtLast = r_.getRange( thisRange );
        Logger::out( LOG_FOR_DEBUGGING ) << "notAtLast=" << notAtLast << endl;
        Logger::out( LOG_FOR_DEBUGGING ) << "thisRange.pos_=" << thisRange.pos_ << endl;

        if ( ( notAtLast == false ) || ( ( thisRange.pos_ & matchFlag ) != 0 ) ) break;


        skipIfNecessary( thisRange, currentPos_, ( *inBwt_ ), countsSoFar_ );
        // count children
        countsThisRange.clear();

        IntervalType errorIntervalType[alphabetSize];
        vector<LetterNumber> correctionBwtPosns[alphabetSize];
        vector<LetterNumber> errBwtPosns[alphabetSize];

        for ( int i = 0; i < alphabetSize; i++ )
        {
            errorIntervalType[i] = INTERVAL_TYPE_DEFAULT;
            correctionBwtPosns[i] = vector<LetterNumber>();
            errBwtPosns[i] = vector<LetterNumber>();
        }

        char *bwtSubstring = new char[thisRange.num_];
        ( *inBwt_ )( bwtSubstring, thisRange.num_ );

        dynamic_cast< BwtCorrectorIntervalHandler & >( intervalHandler_ ).foundInAOnly(
            pileNum,
            countsSoFar_,
            bwtSubstring,
            countsThisRange,
            thisRange,
            propagateInterval_,
            errorIntervalType,
            correctionBwtPosns,
            errBwtPosns,
            cycle_
        );

        delete[] bwtSubstring;

        // Sequence numbers extraction
        if ( doesPropagateBkptToSeqNumInSet_ )
        {
            if ( countsThisRange.count_[0] > 0 )
            {
                // cout << "READNUM"
                //      << "(" << thisRange.isBkptExtension_ << ")"
                //      << " " << i
                //      << " " << ( thisRange.pos_ + 1 )
                //      << " " << thisRange.num_
                //      << " " << countsThisRange.count_[0]
                //      << " " << countsSoFar_.count_[0]
                //      << endl;

            }
            for ( int l( 1 ); l < alphabetSize; l++ )
                propagateInterval_[l] = ( countsThisRange.count_[l] > 0 );
            // cout << "propagated!" << endl;
        }

        // add ranges for any children


#ifdef PROPAGATE_PREFIX
        hasChild = false;
#endif

        for ( int l( 1 ); l < alphabetSize; l++ )
        {
            //       if (countsThisRangeA.count_[l]>=minOcc)
            if ( propagateInterval_[l] == true )
            {
#ifdef PROPAGATE_PREFIX
                if ( hasChild == false )
                {
                    // assert(thisWord.size()==thisRangeA.word_.size()+1);
                    thisWord.replace( 1, thisRange.word_.size(), thisRange.word_ );
                } // ~if
                thisWord[0] = alphabet[l];
#endif
                hasChild = true;

                ErrorCorrectionRange newRange(
                    thisWord,
                    countsSoFar_.count_[l],
                    countsThisRange.count_[l],
                    thisRange.isBkptExtension_,
                    errorIntervalType[l],
                    correctionBwtPosns[l],
                    errBwtPosns[l]
                );
                if ( noComparisonSkip_ || !r_.isRangeKnown( newRange, l, pileNum, subset_, cycle_ ) )
                    r_.addRange( newRange, l, pileNum, subset_, cycle_ );
            } // ~if
        } // ~for l

        if ( hasChild == false )
        {
            //  if no children, print word itself
#ifdef OLD
            cout << "GOLD ";
#ifdef PROPAGATE_PREFIX
            cout << thisRange.word_;
#endif
            cout << " " << thisRange.num_ << endl;
#endif
            numSingletonRanges_++;
        } // ~if

        countsSoFar_ += countsThisRange;
        currentPos_ += thisRange.num_;

        numRanges_++;
    } // ~while notAtLast
} // ~OneBwtBackTracker::operator()
