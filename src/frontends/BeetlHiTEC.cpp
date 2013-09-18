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

#include "BeetlCorrect.hh"

#include "BCRext.hh"
#include "BCRexternalBWT.hh"
#include "Common.hh"
#include "DatasetMetadata.hh"
#include "parameters/BwtParameters.hh"
#include "config.h"
#include "libzoo/cli/Common.hh"
#include "libzoo/util/Logger.hh"
#include "countWords/RangeStore.hh"
#include "countWords/IntervalHandlerReference.hh"
#include "errors/OneBwtBackTracker.hh"
#include "errors/WitnessReader.hh"
#include "errors/HiTECStats.hh"
#include "errors/HiTECParameters.hh"
#include "errors/HiTEC.hh"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;
using namespace BeetlBwtParameters;

int main( const int argc, const char **argv )
{
    HiTECParameters params;
    if ( !params.parseArgv( argc, argv ) || params["help"] == 1 || !params.chechRequiredParameters() )
    {
        params.printUsage();
        exit( 1 );
    }

    string readsFile = params.getStringValue( "input filename" );


    double errorRate = ( double )params.getValue( "error rate" ) / ( double )1000000;
    int genomeLength = params.getValue( "genome length" );

    HiTEC *corrector = new HiTEC(
        readsFile,
        errorRate,
        genomeLength
    );

    corrector->showExecutionPlan();

    if ( !params["don't run"].isSet() )
        corrector->run();

    delete corrector;
    return 0;
}
