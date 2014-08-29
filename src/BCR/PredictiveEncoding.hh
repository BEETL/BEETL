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
