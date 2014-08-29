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

#ifndef INCLUDED_EXTENDER_HH
#define INCLUDED_EXTENDER_HH

#include "search/IntervalFile.hh"
#include "search/ExtenderIntervalHandler.hh"
#include "Algorithm.hh"
#include "Config.hh"
#include "LetterCount.hh"
#include "OneBwtBackTracker.hh"
#include "RangeStore.hh"
#include "Types.hh"
#include "parameters/ExtendParameters.hh"

#include <string>
#include <math.h>

using std::string;

class Extender : public Algorithm
{

public:
    Extender(
        const ExtendParameters &extendParams
    );
    virtual ~Extender() {}
    void run( void );

private:
    void fillRangeStore( RangeStoreExternal &rangeStore, const LetterCountEachPile &countsPerPile, const LetterCountEachPile &countsCumulative );

    const ExtendParameters &extendParams_;
};

#endif
