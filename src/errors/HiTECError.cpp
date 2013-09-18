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

#include "HiTECError.hh"
//NOTE! Designed for 'single bwt backtracking' only, in which only 'foundInAOnly' behaviour is implemented
using namespace std;


void HiTECError::process
( const int pileNum,
  const LetterCount &countsSoFar,
  const LetterCount &countsThisRange,
  AlphabetFlag &propagateInterval )
{
    if ( countsThisRange.count_[0] > 0 )
    {
        errorsGathered_.push_back(
            ReadErrorLocation(
                countsSoFar.count_[0],
                1,
                word_.at( word_.length() - 1 )
            )
        );
    }

} // ~foundInAOnly
