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

#include "Sorting.hh"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>

//#define USE_OPENMP
#ifdef USE_OPENMP
#include <omp.h>
#include <parallel/algorithm>
#endif //ifdef USE_OPENMP


bool cmpSortEl ( sortElement a, sortElement b )
{
    if ( a.pileN == b.pileN )
        return ( a.posN < b.posN );
    else
        return ( a.pileN < b.pileN );
}

void quickSort( vector< sortElement > &v )
{
#ifdef USE_OPENMP
    __gnu_parallel::sort( v.begin(), v.end(), cmpSortEl );
#else //ifdef USE_OPENMP
    sort( v.begin(), v.end(), cmpSortEl );
#endif //ifdef USE_OPENMP
}
