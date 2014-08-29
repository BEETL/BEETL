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

#include "Sorting.hh"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#include <parallel/algorithm>
#endif //ifdef _OPENMP


bool cmpSortEl ( sortElement a, sortElement b )
{
    if ( a.pileN == b.pileN )
        return ( a.posN < b.posN );
    else
        return ( a.pileN < b.pileN );
}

void quickSort( vector< sortElement > &v )
{
#ifdef _OPENMP
    __gnu_parallel::sort( v.begin(), v.end(), cmpSortEl );
#else //ifdef _OPENMP
    sort( v.begin(), v.end(), cmpSortEl );
#endif //ifdef _OPENMP
}
