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

#include "BwtCorrectorIntervalHandler.hh"

#include "ErrorCorrectionRange.hh"

using namespace std;


bool BwtCorrectorIntervalHandler::defaultDetermineErrors( LetterCount intervalLetterCount, int &correct )
{
    bool hasErrors = false;
    bool correctLetterSet = false;
    for ( int i = 1; i < alphabetSize; i++ )
        if ( intervalLetterCount.count_[i] >= minOccurrences_ )
            if ( correctLetterSet )
                return false;
            else
            {
                correct = i;
                correctLetterSet = true;
            }
        else if ( intervalLetterCount.count_[i] > 0 )
            hasErrors = true;

    if ( ( correctLetterSet == false ) || ( hasErrors == false ) )
        return false;
    else
        return true;
}

void BwtCorrectorIntervalHandler::foundInAOnly(
    const int pileNum,
    const LetterCount &countsSoFarA,
    const LetterCount &countsThisRangeA,
    const char *bwtSubstring,
    Range &thisRangeBaseA,
    AlphabetFlag &propagateIntervalA,
    int cycle
)
{
    ErrorCorrectionRange &thisRangeA = dynamic_cast< ErrorCorrectionRange & >( thisRangeBaseA );
    thisRangeA.clearDataForSubIntervals();

    int correct;
    int intervalWordLength_ = cycle;

    // Set default propagation + we never propagate '$' signs
    propagateIntervalA[0] = false;
    for ( int l( 1 ); l < alphabetSize; l++ )
        propagateIntervalA[l] = ( countsThisRangeA.count_[l] > 0 );

    if ( thisRangeA.data_.errorIntervalType == INTERVAL_TYPE_DEFAULT )
        if ( thisRangeA.num_ <= minOccurrences_ )
        {
            for ( int i = 0; i < alphabetSize; i++ )
                propagateIntervalA[i] = false;
            return;
        }

    if ( intervalWordLength_ < minWitnessLength_ )
        return;

    if ( thisRangeA.data_.errorIntervalType == INTERVAL_TYPE_ERROR )
    {
        //non-dollar backward extensions of the error interval are flagged as error intervals with
        //same bwtpos...

        int totalRangeSize = 0;
        for ( int i = 0; i < alphabetSize; i++ )
            totalRangeSize += countsThisRangeA.count_[i];

        int totalErrors = thisRangeA.data_.errorsForBwtPosns.size();
        assert( totalErrors == totalRangeSize );

        //we may have more than one dollar in this interval, hence the end of more than one read.
        //thus we must keep tally of the dollars as we scan through the BWT substring to find where they are
        //(we need to do this to find the BWT positions of the dollars, so we can look them up in the errorStore_ and
        //update the seqNum field of their error object). this way we make sure each error object gets assigned to
        //its intended read.
        int dollarCount = 0;

        assert( needSubstring );
        for ( int relPos = 0; relPos < totalErrors; relPos++ )
        {
            if ( bwtSubstring[relPos] == alphabet[0] )
            {
                if ( errorStore_[thisRangeA.data_.errorsForBwtPosns[relPos]].seqNum == -1 )
                {
                    errorStore_[thisRangeA.data_.errorsForBwtPosns[relPos]].seqNum = countsSoFarA.count_[0] + dollarCount;
                    errorStore_[thisRangeA.data_.errorsForBwtPosns[relPos]].readEnd = intervalWordLength_;
                }
                dollarCount++;
            }

            //non-dollar letters need to be backward extended so we can get closer to finding their terminating characters
            //hence, we.... tag the 'A' (for example) extension of this error interval with the BWT positions
            for ( int i = 1; i < alphabetSize; i++ )
            {
                // tag the 'A' (for example) extension of this error interval with the BWT positions (of the original
                // error containing range) corresponding to 'A's in this range.
                thisRangeA.getDataForSubInterval( i ).errorIntervalType = INTERVAL_TYPE_ERROR;
                if ( bwtSubstring[relPos] == alphabet[i] )
                    thisRangeA.getDataForSubInterval( i ).errorsForBwtPosns.push_back( thisRangeA.data_.errorsForBwtPosns[relPos] );
            }
        }
        return;
    }

    if ( defaultDetermineErrors( countsThisRangeA, correct ) )
    {
        assert( needSubstring );
        //'correct' is the letter we believe is the correct one for the interval
        //when scanning along the bwtSubstring, anything other than a dollar or this 'correct' letter is treated as an error
        for ( uint relativePos = 0; relativePos < thisRangeA.num_; relativePos++ )
            if ( bwtSubstring[relativePos] != alphabet[correct] && bwtSubstring[relativePos] != '$' )
            {
                //we have identified that bwtSubstring[relativePos] is not the correct letter for this range and not a dollar
                //so call it putative error...
                char putativeError = bwtSubstring[relativePos];

                LetterNumber errBwtPos = relativePos;
                for ( int i = 0; i < alphabetSize; i++ )
                    errBwtPos += countsSoFarA.count_[i];

                if ( errorStore_.find( errBwtPos ) == errorStore_.end() )
                {
                    ErrorInfo newError;
                    newError.firstCycle = intervalWordLength_;
                    newError.lastCycle = intervalWordLength_;
                    newError.corrector += alphabet[correct];
                    errorStore_[errBwtPos] = newError;

                    //finding putative error for the first time, so flag next generation of intervals...

                    //flag the 'extend by putativeError' interval as error type...
                    thisRangeA.getDataForSubInterval( whichPile[( int )putativeError] ).errorIntervalType = INTERVAL_TYPE_ERROR;

                    //flag the 'extend by correct letter' interval as corrector type...
                    thisRangeA.getDataForSubInterval( correct ).errorIntervalType = INTERVAL_TYPE_CORRECTOR;

                    //tag the 'extend by putativeError' interval and 'extend by correct letter' interval
                    //with position of putativeError in the BWT...
                    thisRangeA.getDataForSubInterval( whichPile[( int )putativeError] ).errorsForBwtPosns.push_back( errBwtPos );
                    thisRangeA.getDataForSubInterval( correct ).correctionForBwtPosns.push_back( errBwtPos );
                }
                else
                {
                    //re-finding, so don't flag any intervals... just update 'last cycle we saw this error'
                    errorStore_[errBwtPos].lastCycle = intervalWordLength_;
                }
            }
    }

    //corrector interval...
    if ( thisRangeA.data_.errorIntervalType == INTERVAL_TYPE_CORRECTOR )
    {
        //does one backward extension 'dominate' the others?
        int rangelength = 0;
        for ( int i = 1; i < alphabetSize; i++ )
            rangelength += countsThisRangeA.count_[i];
        int dominator = 0;

        for ( int i = 1; i < alphabetSize; i++ )
            if ( countsThisRangeA.count_[i] >= minOccurrences_ )
                dominator = i;

        if ( ( dominator > 0 ) && ( rangelength > 0 ) )
        {
            //if so, we'll flag the extension by the 'dominator' as corrector type...
            thisRangeA.getDataForSubInterval( dominator ).errorIntervalType = INTERVAL_TYPE_CORRECTOR;

            //now for each BWT position with which this corrector interval is tagged, we'll
            //extend the corrector string of the corresponding error objects in the error store...
            //and add it to BWT positions with which the extension by the 'dominator' is tagged...
            for ( uint errNo = 0; errNo < thisRangeA.data_.correctionForBwtPosns.size(); ++errNo )
            {
                errorStore_[thisRangeA.data_.correctionForBwtPosns[errNo]].corrector += alphabet[dominator];
                thisRangeA.getDataForSubInterval( dominator ).correctionForBwtPosns.push_back( thisRangeA.data_.correctionForBwtPosns[errNo] );
            }
        }
    }
}
