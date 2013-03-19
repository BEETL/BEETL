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

#include "BeetlUnbwt.hh"

#include "BCRexternalBWT.hh"
#include "Common.hh"
#include "Logger.hh"
#include "config.h"
#include "parameters/UnbwtParameters.hh"

#include <cstdlib>
#include <iostream>
#include <sstream>

using namespace std;


void launchBeetlUnbwt( UnbwtParameters &unbwtParams, const string &inputFilename, const string &outputFilename )
{
    int bcrMode = 1 ; // 1=decode BWT
    CompressionFormatType outputCompression = compressionIncrementalRunLength; // not used
    BCRexternalBWT bwt( ( char * )inputFilename.c_str(), ( char * )outputFilename.c_str(), bcrMode, outputCompression, &unbwtParams );
    // automatically calls: result = unbuildBCR( file1, fileOutBwt, intermediateCycFiles, fileOutput );
    return;


    ostringstream cmdLineParams;
    cmdLineParams << " bcr -m 1 "
                  << " -i " << inputFilename
                  << " -o " << outputFilename;

    launchBeetl( cmdLineParams.str() );
}


void printUsage()
{
    cout << "Parameters:" << endl;
    cout << "    --input (-i)             Input prefix" << endl;
    cout << "    --output (-o)            Output file name" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "    --decode-direction (-d)  = backward [forward|backward]" << endl;
    cout << "    --use-vector             = on [on|off]: Pre-compute a sampling vector of the BWT segments. More RAM, faster." << endl;
    cout << "    --verbose                = quiet [quiet|verbose|very-verbose|debug] or [0|1|2|3]" << endl;
    cout << "    --help (-h)              Help" << endl;
    cout << endl;
}


struct BeetlUnbwtArguments
{
    string argInput;
    string argOutput;
    string argVerbosityLevel;
};

int main( const int argc, const char **argv )
{
    // Generated using: http://patorjk.com/software/taag/#p=display&f=Soft&t=BEETL%20unbwt
    cout << ",-----.  ,------.,------.,--------.,--.                       ,-----.  ,--.   ,--.,--------. " << endl;
    cout << "|  |) /_ |  .---'|  .---''--.  .--'|  |       ,--.,--.,--,--, |  |) /_ |  |   |  |'--.  .--' " << endl;
    cout << "|  .-.  \\|  `--, |  `--,    |  |   |  |       |  ||  ||      \\|  .-.  \\|  |.'.|  |   |  |    " << endl;
    cout << "|  '--' /|  `---.|  `---.   |  |   |  '--.    '  ''  '|  ||  ||  '--' /|   ,'.   |   |  |    " << endl;
    cout << "`------' `------'`------'   `--'   `-----'     `----' `--''--'`------' '--'   '--'   `--'    " << endl;
    cout << "Version " << PACKAGE_VERSION << endl;
    cout << endl;

    cout << "Command called:" << endl << "   ";
    for ( int i = 0; i < argc; ++i )
    {
        cout << " " << argv[i];
    }
    cout << "\n" << endl;

    // Todo: reduce these 2 structures to 1
    BeetlUnbwtArguments args;
    UnbwtParameters params;

    for ( int i = 1; i < argc; ++i )
    {
        if ( isNextArgument( "-h", "--help", argc, argv, i ) )
        {
            printUsage();
            exit( 0 );
        }
        else if ( isNextArgument   ( "-i", "--input"               , argc, argv, i, &args.argInput                           ) ) {}
        else if ( isNextArgument   ( "-o", "--output"              , argc, argv, i, &args.argOutput                          ) ) {}
        else if ( parseNextArgument( "-d", "--decode-direction"    , argc, argv, i, params, UNBWT_OPTION_DECODE_DIRECTION    ) ) {}
        else if ( parseNextArgument( ""  , "--use-vector"          , argc, argv, i, params, UNBWT_OPTION_USE_VECTOR          ) ) {}
        else if ( isNextArgument   ( ""  , "--verbose"             , argc, argv, i, &args.argVerbosityLevel                  ) )
        {
            Logger::setVerbosity( args.argVerbosityLevel );
        }
        else
        {
            cerr << "Error: Invalid parameter: " << argv[i] << "\n" << endl;
            printUsage();
            exit( 1 );
        }
    }

    // Checking for required parameters
    if ( args.argInput.empty() || args.argOutput.empty() )
    {
        cerr << "Error: Missing arguments: --input and --output are required\n" << endl;
        printUsage();
        exit( 1 );
    }

    // Use default parameter values where needed
    if ( params.getValue( UNBWT_OPTION_DECODE_DIRECTION ) == MULTIPLE_OPTIONS )
        params.set( UNBWT_OPTION_DECODE_DIRECTION, DECODE_DIRECTION_BACKWARD );

    if ( params.getValue( UNBWT_OPTION_USE_VECTOR ) == MULTIPLE_OPTIONS )
        params.set( UNBWT_OPTION_USE_VECTOR, USE_VECTOR_ON );

    // Launch
    launchBeetlUnbwt( params, args.argInput, args.argOutput );

    return 0;
}
