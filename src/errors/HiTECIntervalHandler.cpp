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

#include "HiTECIntervalHandler.hh"
//NOTE! Designed for 'single bwt backtracking' only, in which only 'foundInAOnly' behaviour is implemented
using namespace std;


void ErrorIntervalHandler::foundInBoth
( const int pileNum,
  const LetterCount &countsThisRangeA, const LetterCount &countsThisRangeB,
  const Range &thisRangeA, const Range &thisRangeB,
  AlphabetFlag &propagateIntervalA, AlphabetFlag &propagateIntervalB,
  bool &isBreakpointDetected )
{

} // ~foundInBoth

void ErrorIntervalHandler::foundInAOnly
( const int pileNum,
  const LetterCount &countsSoFarA,
  const LetterCount &countsThisRangeA,
  const Range &thisRangeA,
  AlphabetFlag &propagateIntervalA )
{
    //seems like we didn't need any of the stuf that got copied from 'IntervalHandlerReference' so
    //just deleted it - but be prepared to cutnpaste back!

    //todo: i think when there are two or more identical reads this is only going to tell you
    //there is an error in one of them...
    if ( countsThisRangeA.count_[0] > 0 )
    {
        stringstream ss;
#ifdef PROPAGATE_PREFIX
        ss << thisRangeA.word_.at( thisRangeA.word_.length() - 1 );
#else
        ss << alphabet[pileNum];
#endif

                errorsGathered_.push_back(
                    HiTECErrorLocation(
                        countsSoFarA.count_[0],
                        iterationNumber_,
                        ss.str()[0]
                    )
                );
    }

} // ~foundInAOnly

        void ErrorIntervalHandler::foundInBOnly
        ( const int pileNum,
          const LetterCount &countsSoFarB,
          const LetterCount &countsThisRangeB,
          const Range &thisRangeB,
          AlphabetFlag &propagateIntervalB )
{
    // TBD
}

