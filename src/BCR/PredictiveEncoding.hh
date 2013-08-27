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

#ifndef PREDICTIVE_ENCODING_HH
#define PREDICTIVE_ENCODING_HH

#include "libzoo/util/AutoGrowVector.hh"

#include <string>
#include <vector>

using std::string;
using std::vector;


class PredictionStatistics
{
public:
    PredictionStatistics();

    void add( const bool isQualityRemoved, const char qual );
    void outputToFile( const string &filename );

private:
    int removedCount_;
    int notRemovedCount_;
    AutoGrowVector<int> removedQualities_;
};


#endif // PREDICTIVE_ENCODING_HH
