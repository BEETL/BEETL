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

#include "ExtenderIntervalHandler.hh"

#include "IntervalFile.hh"
#include "libzoo/util/Logger.hh"

#include <algorithm>

using namespace std;

extern vector<string> kmerList2;

ExtenderIntervalHandler::ExtenderIntervalHandler()
{
}

void ExtenderIntervalHandler::foundInAOnly
( const int pileNum,
  const LetterCount &countsSoFarA,
  const LetterCount &countsThisRangeA,
  const char *bwtSubstring,
  Range &thisRangeBaseA,
  AlphabetFlag &propagateIntervalA,
  const int cycle
)
{
    if ( countsThisRangeA.count_[0] )
    {
        IntervalRecord *rec = reinterpret_cast< IntervalRecord * >( thisRangeBaseA.userData_ );

        Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "$ signs detected for " << *rec << ": " << countsThisRangeA.count_[0] << " items from " << countsSoFarA.count_[0] << endl;

        for ( LetterNumber i = 0; i < countsThisRangeA.count_[0]; ++i )
            rec->dollarSignPositions.push_back( countsSoFarA.count_[0] + i );
    }
    //    for ( int i = 1; i < alphabetSize; i++ )
    //        propagateIntervalA[i] = true;

}
