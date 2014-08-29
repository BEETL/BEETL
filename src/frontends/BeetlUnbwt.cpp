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

#include "BeetlUnbwt.hh"

#include "BCRexternalBWT.hh"
#include "config.h"
#include "parameters/UnbwtParameters.hh"
#include "libzoo/cli/Common.hh"
#include "libzoo/util/Logger.hh"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string.h>

using namespace std;
using namespace BeetlUnbwtParameters;


UnbwtParameters params;

void launchBeetlUnbwt()
{
    const string &inputFilename = params["input filename prefix"];
    const string &outputFilename = params["output filename"];

    Logger::out() << "\nLaunching the following configuration of Beetl-unbwt:" << endl;
    params.print( Logger::out(), false );
    Logger::out() << endl;

    int bcrMode = 1 ; // 1=decode BWT
    CompressionFormatType outputCompression = compressionIncrementalRunLength; // not used
    BCRexternalBWT bwt( ( char * )inputFilename.c_str(), ( char * )outputFilename.c_str(), bcrMode, outputCompression, &params );
    // automatically calls: result = unbuildBCR( file1, fileOutBwt, intermediateCycFiles, fileOutput );
}


void printUsage()
{
    params.printUsage();

    cout << "Notes:" << endl;
    cout << "    Input must be a set of ASCII-encoded BWT files (not run-length-encoded)" << endl;
    cout << "    Fastq output requires {input}-Q0x quality files to be present" << endl;
    cout << endl;
}


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

    if ( !params.parseArgv( argc, argv ) || params["help"] == 1 || !params.chechRequiredParameters() )
    {
        printUsage();
        exit( params["help"] == 0 );
    }

    // Use default parameter values where needed
    params.commitDefaultValues();

    // Auto-detection of missing arguments
    if ( !params["output format"].isSet() )
    {
        const string &filename = params["output filename"];
        string fileFormat = detectFileFormat( filename );
        if ( fileFormat.empty() )
        {
            cerr << "Error: file format not recognised for " << filename << endl;
            exit( -1 );
        }
        params["output format"] = fileFormat;
    }
    checkFileFormat( params["output filename"], params["output format"] );

    if ( !params["input format"].isSet() )
    {
        const string &bwtPrefix = params["input filename prefix"];
        vector<string> filenames;
        bool isBwtCompressed;
        string availableFileLetters;
        detectInputBwtProperties( bwtPrefix, filenames, isBwtCompressed, availableFileLetters );

        if ( filenames.size() < 2 )
        {
            cerr << "Error: too few input files detected (run with -vv for more details)" << endl;
            exit( -1 );
        }

        if ( isBwtCompressed )
        {
            cerr << "Error: BWT files don't seem to be in ASCII format (they probably got created as run-length-encoded)" << endl;
            exit( -1 );
        }
        else
        {
            params["input format"] = "BWT_ASCII";
        }

        // Check that {prefix}-Q0* files are present if output is fastq
        if ( params["output format"] == "fastq" )
        {
            assert( strchr( availableFileLetters.c_str(), 'Q' ) && "{input}-Q0x files cannot be found (those qualities are required for fastq output)" );
        }
    }

    // Launch
    launchBeetlUnbwt();

    return 0;
}
