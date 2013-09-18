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

#ifndef INCLUDED_HITECSTATS_HH
#define INCLUDED_HITECSTATS_HH

#include "Config.hh"
#include "LetterCount.hh"

#include "Types.hh"

#include <string>
#include <math.h>

using std::string;

#define HITEC_MAX_ITERATIONS 9

class HiTECStats
{

public:
    HiTECStats(
        const double errorRate,
        const double genomeLength,
        const int numberOfReads,
        const int readLength
    );
    int Calculate_wm();
    int Calculate_wm( int maxValue );
    int Calculate_wM( double maxDestructibleRate = 0.0001 );
    int CalculateSupport( int witnessLength );

private:
    const int numberOfReads_;
    const int readLength_;
    const double genomeLength_;
    const double errorRate_;

    double probDestructible( int witnessLength );
    double expectedDestructible( int witnessLength );
    double expectedUncorrectable( int witnessLength );
    double expectedCorrectWitnessNeighbourPairs( int support, int witnessLength );
    double expectedIncorrectWitnessNeighbourPairs( int support, int witnessLength );
    double expectedErroneousReads();
};

#endif
