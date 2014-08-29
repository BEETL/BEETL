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

#include "BeetlCorrect.hh"

#include "BCRext.hh"
#include "BCRexternalBWT.hh"
#include "Common.hh"
#include "DatasetMetadata.hh"
#include "parameters/BwtParameters.hh"
#include "config.h"
#include "libzoo/cli/Common.hh"
#include "libzoo/util/Logger.hh"
#include "errors/WitnessReader.hh"
#include "errors/HiTECStats.hh"
#include "errors/BwtCorrectorParameters.hh"
#include "errors/BwtCorrector.hh"
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

    cout << ",-----.  ,------.,------.,--------.,--.        ,-----. ,-----. ,------. ,------. ,------. ,-----.,--------." << endl;
    cout << "|  |) /_ |  .---'|  .---''--.  .--'|  |       '  .--./'  .-.  '|  .--. '|  .--. '|  .---''  .--./'--.  .--' " << endl;
    cout << "|  .-.  \\|  `--, |  `--,    |  |   |  |       |  |    |  | |  ||  '--'.'|  '--'.'|  `--, |  |       |  |    " << endl;
    cout << "|  '--' /|  `---.|  `---.   |  |   |  '--.    '  '--'\\'  '-'  '|  |\\  \\ |  |\\  \\ |  `---.'  '--'\\   |  |    " << endl;
    cout << "`------' `------'`------'   `--'   `-----'     `-----' `-----' `--' '--'`--' '--'`------' `-----'   `--'    " << endl;
    cout << endl;

    BwtCorrectorParameters params;
    if ( !params.parseArgv( argc, argv ) || params["help"] == 1 || !params.chechRequiredParameters() )
    {
        params.printUsage();
        exit( params["help"] == 0 );
    }

    // Use default parameter values where needed
    params.commitDefaultValues();

    string indexPrefix = params.getStringValue( "input filename" );

    int readLength = params.getValue( "read length" );

    bool compressed;
    vector<string> pileNames;
    string dummyStr;
    BwtReaderBase *dollarPile;
    detectInputBwtProperties( indexPrefix, pileNames, compressed, dummyStr );
    dollarPile = instantiateBwtPileReader( pileNames[0], params.getStringValue( "use shm" ) );

    int numReads = 0;
    LetterCount lc;
    dollarPile->readAndCount( lc );
    for ( int i = 0; i < alphabetSize; i++ )
        numReads += lc.count_[i];
    //must divide number of reads in half, as it (should) include reverse complements...
    numReads /= 2;

    delete dollarPile;

    double errorRate = ( double )params.getValue( "error rate" ) / ( double )1000000;
    int genomeLength = params.getValue( "genome length" );

    int minWitnessLength;

    if ( params["min witness length"].isSet() )
        minWitnessLength = params.getValue( "min witness length" );
    else
    {
        HiTECStats stats(
            errorRate,
            genomeLength,
            numReads,
            readLength
        );
        minWitnessLength = stats.Calculate_wm() - 1;
    }

    int minSupport = 0;

    if ( params["min support"].isSet() )
        minSupport = params.getValue( "min support" );

    BwtCorrector *corrector = new BwtCorrector(
        indexPrefix,
        params.getStringValue( "corrections output filename" ),
        numReads,
        readLength,
        errorRate,
        genomeLength,
        minWitnessLength,
        params.getStringValue( "subset" ),
        &params,
        minSupport
    );

    corrector->showExecutionPlan();

    if ( !params["don't run"].isSet() )
        corrector->run();

    delete corrector;
    return 0;
}
