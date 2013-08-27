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

#include "PredictiveEncoding.hh"

#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;


PredictionStatistics::PredictionStatistics()
    : removedCount_( 0 )
    , notRemovedCount_( 0 )
{
}

void PredictionStatistics::add( const bool isQualityRemoved, const char qual )
{
    assert( qual >= 33 &&  "Error: Invalid Phred Q-score" );

    if ( isQualityRemoved )
    {
        assert( qual != 33 &&  "Error: Trying to remove a zero Q-score" );
        ++removedCount_;
        ++( removedQualities_[qual - 33] );
    }
    else
    {
        ++notRemovedCount_;
    }
}

void PredictionStatistics::outputToFile( const string &filename )
{
    ofstream os( filename.c_str() );
    os << "Removed Qscores\t" << removedCount_ << endl;
    os << "Kept Qscores\t" << notRemovedCount_ << endl;

    double removedCount = 0;
    double removedSum = 0;
    double associatedErrorRateSum = 0;
    for ( unsigned int i = 1; i < removedQualities_.size(); ++i )
    {
        os << "removed Qscore\t" << i << "\t" << removedQualities_[i] << endl;
        removedCount += removedQualities_[i];
        removedSum += removedQualities_[i] * i;
        associatedErrorRateSum += removedQualities_[i] * pow( 10, -( double )i / 10.0 );
    }

    assert( removedCount == removedCount_ );
    double bestReplacementQuality = -10 * log10( associatedErrorRateSum / removedCount );
    os << "Mean removed Qscore\t" << ( removedSum / removedCount ) << "\t" << char( round( 33 + removedSum / removedCount ) ) << endl;
    os << "Mean associated error rate\t" << ( associatedErrorRateSum / removedCount ) << endl;
    os << "Best replacement QScore\t" << bestReplacementQuality << "\t" << char( round( 33 + bestReplacementQuality ) ) << endl;
}
