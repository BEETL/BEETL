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

#include "BeetlSearch.hh"

#include "BCRexternalBWT.hh"
#include "Tools.hh"
#include "config.h"
#include "parameters/SearchParameters.hh"
#include "libzoo/cli/Common.hh"
#include "libzoo/util/Logger.hh"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <unistd.h>

using namespace std;


void printUsage()
{
    cout << "Parameters:" << endl;
    cout << "    --input  (-i)            Input filename prefix (i.e. BWT files are \"prefix-B0[0-6]\")" << endl;
    cout << "    --output (-o)            = \"searchedKmers_positions\"" << endl;
    cout << "    --kmer   (-k)            Single k-mer string to be searched for" << endl;
    cout << " or --kmers  (-j)            File containing one k-mer per line" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "    --verbosity              = normal [quiet|normal|verbose|very-verbose|debug] or [0|1|2|3|4]" << endl;
    cout << "    -v / -vv                 Shortcuts to --verbosity = verbose / very-verbose" << endl;
    cout << "    --help (-h)              Help" << endl;
    cout << endl;
    cout << "Notes:" << endl;
    cout << "    BWT files must be in ASCII format" << endl;
    cout << "    An intermediate file \"searchedKmers\" will be created, so this file name must be available" << endl;
    cout << endl;
}

struct BeetlSearchArguments
{
    string argInput;
    string argOutput;
    string argKmer;
    string argKmers;
    string argVerbosityLevel;

    BeetlSearchArguments()
        : argOutput( "searchedKmers_positions" )
    {}
};

void launchBeetlSearch( SearchParameters &searchParams, const BeetlSearchArguments &args )
{
    // Try to guess if run-length-encoded based on first few bytes
    bool inputCompressed = false;
    string B01Filename = args.argInput + "-B01";
    ifstream is( B01Filename.c_str() );

    for ( int i = 0; i < 10; ++i )
    {
        char c = 'A';
        is.get( c );
        switch ( toupper( c ) )
        {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'N':
                break;
            default:
                inputCompressed = true;
        }
    }
    if ( inputCompressed )
    {
        cerr << "Error: Set A detected as non-ASCII" << endl;
        exit ( -1 );
    }

    // Use the file "searchedKmers" as an intermediate buffer
    const string intermediateFilename = "searchedKmers";
    if ( readWriteCheck( intermediateFilename.c_str(), false, false ) )
    {
        cerr << "Error: intermediate file \"" << intermediateFilename << "\" already exists. Aborting." << endl;
        exit( -1 );
    }
    if ( !args.argKmer.empty() )
    {
        ofstream os( intermediateFilename.c_str() );
        os << args.argKmer << endl;
    }
    else
    {
        if ( link( args.argKmers.c_str(), intermediateFilename.c_str() ) != 0 )
        {
            cerr << "Error: Cannot create link " << intermediateFilename << " -> " << args.argKmers << endl;
        }
    }

    // Launch
    int bcrMode = 2 ; // 2=search BWT
    CompressionFormatType outputCompression = compressionASCII; // not used
    BCRexternalBWT bwt( ( char * )args.argInput.c_str(), ( char * )args.argOutput.c_str(), bcrMode, outputCompression, &searchParams );

    // Cleanup
    unlink( intermediateFilename.c_str() );
}

int main( const int argc, const char **argv )
{
    // Generated using: http://patorjk.com/software/taag/#p=display&f=Soft&t=BEETL%20search
    cout << ",-----.  ,------.,------.,--------.,--.                                          ,--.      " << endl;
    cout << "|  |) /_ |  .---'|  .---''--.  .--'|  |        ,---.  ,---.  ,--,--.,--.--. ,---.|  ,---.  " << endl;
    cout << "|  .-.  \\|  `--, |  `--,    |  |   |  |       (  .-' | .-. :' ,-.  ||  .--'| .--'|  .-.  | " << endl;
    cout << "|  '--' /|  `---.|  `---.   |  |   |  '--.    .-'  `)\\   --.\\ '-'  ||  |   \\ `--.|  | |  | " << endl;
    cout << "`------' `------'`------'   `--'   `-----'    `----'  `----' `--`--'`--'    `---'`--' `--' " << endl;
    cout << "Version " << PACKAGE_VERSION << endl;
    cout << endl;

    cout << "Command called:" << endl << "   ";
    for ( int i = 0; i < argc; ++i )
    {
        cout << " " << argv[i];
    }
    cout << "\n" << endl;

    BeetlSearchArguments args;
    SearchParameters params;

    for ( int i = 1; i < argc; ++i )
    {
        if ( isNextArgument( "-h", "--help", argc, argv, i ) )
        {
            printUsage();
            exit( 0 );
        }
        else if ( isNextArgument   ( "-i", "--input"               , argc, argv, i, &args.argInput                           ) ) {}
        else if ( isNextArgument   ( "-o", "--output"              , argc, argv, i, &args.argOutput                          ) ) {}
        else if ( isNextArgument   ( "-k", "--kmer"                , argc, argv, i, &args.argKmer                            ) ) {}
        else if ( isNextArgument   ( "-j", "--kmers"               , argc, argv, i, &args.argKmers                           ) ) {}
        else if ( isNextArgument   ( ""  , "--verbosity"           , argc, argv, i, &args.argVerbosityLevel                  ) )
        {
            Logger::setVerbosity( args.argVerbosityLevel );
        }
        else if ( isNextArgument   ( "-v"  , ""                    , argc, argv, i ) )
        {
            Logger::setVerbosity   ( "verbose" );
        }
        else if ( isNextArgument   ( "-vv" , ""                    , argc, argv, i ) )
        {
            Logger::setVerbosity( "very-verbose" );
        }
        else
        {
            cerr << "Error: Invalid parameter: " << argv[i] << "\n" << endl;
            printUsage();
            exit( 1 );
        }
    }

    // Checking for required parameters
    if ( args.argInput.empty() || !( args.argKmer.empty() ^ args.argKmers.empty() ) )
    {
        cerr << "Error: Missing or incorrect arguments: -i is required; -j and -k are mutually exclusive, one of them being required\n" << endl;
        printUsage();
        exit( 1 );
    }

    // Launch
    launchBeetlSearch( params, args );

    return 0;
}
