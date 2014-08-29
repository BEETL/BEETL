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

#include "IntervalHandlerBase.hh"

#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;


void IntervalHandlerBase::countString( char *bwtSubstring, int length, LetterCount &countsThisRange )
{
    for ( int i = 0; i < length; i++ )
        countsThisRange += bwtSubstring[i];
}


void IntervalHandlerBase::createOutputFile( const int subsetThreadNum, const int i, const int j, const int cycle, const string &outputDirectory )
{
#define CONCATENATE_J_PILES
    //    if ( cycle >= ( int )minWordLength_ )
    {
        ostringstream filename;
#ifdef CONCATENATE_J_PILES
        filename << outputDirectory << ( outputDirectory.empty() ? "" : "/" ) << "cycle" << cycle << ".subset" << subsetThreadNum << "." << i;
        outFile_.open( filename.str(), ( j == 1 ) ? ios::out : ios::app );
#else
        filename << outputDirectory << ( outputDirectory.empty() ? "" : "/" ) << "cycle" << cycle << ".subset" << subsetThreadNum << "." << i << "." << j;
        outFile_.open( filename.str() );
#endif
    }
}
