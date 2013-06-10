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

#include "BeetlCompare.hh"

#include "Common.hh"
#include "Logger.hh"
#include "TemporaryFilesManager.hh"
#include "config.h"
#include "countWords/CountWords.hh"
#include "parameters/CompareParameters.hh"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <unistd.h>

using namespace std;
using namespace BeetlCompareParameters;


void printUsage()
{
    cout << "Parameters:" << endl;
    cout << "    --inputA (-a)            Input filename prefix for Set A (such as \"prefix-B0[0-6]\" are Set A's BWT files)" << endl;
    cout << "    --inputB (-b)            Input filename prefix for Set B (such as \"prefix-B0[0-6]\" are Set B's BWT files)" << endl;
    cout << "    --mode (-m)              [split|reference|metagenomics] (See \"Mode\" note below)" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "    --output (-o)            =\"\"          Output filename prefix" << endl;
    cout << "    --max-length (-k)        = 100        Maximal length (number of analysis cycles)" << endl;
    cout << "    --min-occ (-n)           = 1          Minimum number of occurrences (coverage)" << endl;
    cout << "    --inputA-format          = autodetect [bwt_ascii|bwt_rle]" << endl;
    cout << "    --inputB-format          = autodetect [bwt_ascii|bwt_rle]" << endl;
    cout << "    --subset                 = \"\"       Restrict computation to this suffix - Used for distributed computing" << endl;
    cout << "    --temp-directory (-T)    = \".\"      Path for temporary files (hint: choose a fast drive)" << endl;
    cout << "    --verbosity              = normal [quiet|normal|verbose|very-verbose|debug] or [0|1|2|3|4]" << endl;
    cout << "    -v / -vv                 Shortcuts to --verbosity = verbose / very-verbose" << endl;
    cout << "    --help (-h)              Help" << endl;
    cout << endl;
    cout << "Metagenomics mode parameters:" << endl;
    cout << "    --genome-metadata (-c)   (= \"${inputB}-C0\") Input filename \"extended\" prefix for Set B's metadata (for files \"prefix[0-6]\")" << endl;
    cout << "    --taxonomy (-t)          Input filename for Set B's taxonomy information" << endl;
    cout << endl;
    cout << "Metagenomics mode options:" << endl;
    cout << "    --min-kmer-length (-w)   = 50   Minimum length k-mer" << endl;
    cout << "    --report-min-length (-d) = off  [on|off] Report the minimal needed word length for the different taxa in the database" << endl;
    cout << endl;
    cout << "Notes:" << endl;
    cout << "    Mode = split       : Set B are reads of a reference" << endl;
    cout << "    Mode = reference   : Set B is a reference genome" << endl;
    cout << "    Mode = metagenomics: Set B are merged reference genomes. Switches metagenome classifier on." << endl;
    cout << endl;
}

struct BeetlCompareArguments
{
    string argInputA;
    string argInputB;
    string argOutput;
    int argMaxLength;
    int argMinOcc;
    string argSubset;
    string argTempPath;
    string argVerbosityLevel;

    // metagenomics mode
    string argInputC;
    string argTaxonomy;
    int argMinKmerLength;

    BeetlCompareArguments()
        : argMaxLength( 100 )
        , argMinOcc( 1 )
        , argMinKmerLength( 50 )
    {}
};

void launchBeetlCompare( CompareParameters &compareParams, const BeetlCompareArguments &args )
{
    char whichHandler = 'X';
    switch ( compareParams.getValue( COMPARE_OPTION_MODE ) )
    {
        case MODE_SPLIT:
            whichHandler = 's';
            break;
        case MODE_REFERENCE:
            whichHandler = 'r';
            break;
        case MODE_METAGENOMICS:
            whichHandler = 'm';
            break;
    }

    vector<string> setA;
    vector<string> setB;
    vector<string> setC;
    for ( unsigned i = 0; i < 10; ++i )
    {
        stringstream fileA;
        fileA << args.argInputA << "-B0" << i;
        clog << "Probing " << fileA.str() << "..." << endl;
        if ( access( fileA.str().c_str(), R_OK ) == -1 )
            break;
        setA.push_back( fileA.str() );
        clog << "Discovered " << fileA.str() << endl;

        stringstream fileB;
        fileB << args.argInputB << "-B0" << i;
        if ( access( fileB.str().c_str(), R_OK ) == -1 )
        {
            cerr << "Error opening " << fileB.str() << endl;
            exit( -1 );
        }
        setB.push_back( fileB.str() );
        clog << "Discovered " << fileB.str() << endl;

        if ( compareParams.getValue( COMPARE_OPTION_MODE ) == MODE_METAGENOMICS )
        {
            stringstream fileC;
            if ( args.argInputC.empty() )
                fileC << args.argInputB << "-C0" << i;
            else
                fileC << args.argInputC << i;
            if ( access( fileC.str().c_str(), R_OK ) == -1 )
            {
                cerr << "Error opening " << fileC.str() << endl;
                exit( -1 );
            }
            setC.push_back( fileC.str() );
            clog << "Discovered " << fileC.str() << endl;
        }
    }
    assert( setA.size() == setB.size() );
    if ( compareParams.getValue( COMPARE_OPTION_MODE ) == MODE_METAGENOMICS )
        assert( setA.size() == setC.size() );
    if ( setA.size() < 2 )
    {
        cerr << "Error: too few input files detected (all those detected have been reported above)" << endl;
        exit( -1 );
    }

    // Try to guess if run-length-encoded based on first few bytes
    bool inputACompressed = false;
    bool inputBCompressed = false;
    if ( setA.size() > 1 )
    {
        ifstream setAFile2( setA[1].c_str() );
        ifstream setBFile2( setB[1].c_str() );

        for ( int i = 0; i < 10; ++i )
        {
            char c = 'A';
            setAFile2.get( c );
            switch ( toupper( c ) )
            {
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                case 'N':
                    break;
                default:
                    inputACompressed = true;
            }

            c = 'A';
            setBFile2.get( c );
            switch ( toupper( c ) )
            {
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                case 'N':
                    break;
                default:
                    inputBCompressed = true;
            }
        }
    }
    if ( inputACompressed )
        clog << "Set A detected as RLE compressed" << endl;
    else
        clog << "Set A detected as ASCII" << endl;
    if ( inputBCompressed )
        clog << "Set B detected as RLE compressed" << endl;
    else
        clog << "Set B detected as ASCII" << endl;

    bool reportMinLength = ( compareParams.getValue( COMPARE_OPTION_REPORT_MINLENGTH ) == REPORT_MINLENGTH_ON );

    // Initialise temporary directory
    TemporaryFilesManager::get().setTempPath( args.argTempPath );

    clog << "Launching CountWords with "
         << "handler=" << whichHandler
         << ", minOcc=" << args.argMinOcc
         << ", maxLength=" << args.argMaxLength
         << ", #files=" << setA.size()
         << ", reportMinLength=" << reportMinLength
         << ", minKmerLength=" << args.argMinKmerLength
         << ", subset=" << args.argSubset
         << endl;

    Algorithm *pcountWords = new CountWords( inputACompressed,
            inputBCompressed, whichHandler, args.argMinOcc,
            args.argMaxLength, setA, setB, setC, args.argTaxonomy, reportMinLength, args.argMinKmerLength, args.argOutput, args.argSubset );

    // run the "main" method
    pcountWords->run();

    // clean up
    delete pcountWords;
    TemporaryFilesManager::get().cleanup();
}

int main( const int argc, const char **argv )
{
    // Generated using: http://patorjk.com/software/taag/#p=display&f=Soft&t=BEETL%20compare
    cout << ",-----.  ,------.,------.,--------.,--.                                                            " << endl;
    cout << "|  |) /_ |  .---'|  .---''--.  .--'|  |        ,---. ,---. ,--,--,--. ,---.  ,--,--.,--.--. ,---.  " << endl;
    cout << "|  .-.  \\|  `--, |  `--,    |  |   |  |       | .--'| .-. ||        || .-. |' ,-.  ||  .--'| .-. : " << endl;
    cout << "|  '--' /|  `---.|  `---.   |  |   |  '--.    \\ `--.' '-' '|  |  |  || '-' '\\ '-'  ||  |   \\   --. " << endl;
    cout << "`------' `------'`------'   `--'   `-----'     `---' `---' `--`--`--'|  |-'  `--`--'`--'    `----' " << endl;
    cout << "                                                                     `--'                          " << endl;
    cout << "Version " << PACKAGE_VERSION << endl;
    cout << endl;

    cout << "Command called:" << endl << "   ";
    for ( int i = 0; i < argc; ++i )
    {
        cout << " " << argv[i];
    }
    cout << "\n" << endl;

    // Todo: reduce these 2 structures to 1
    BeetlCompareArguments args;
    CompareParameters params;

    for ( int i = 1; i < argc; ++i )
    {
        if ( isNextArgument( "-h", "--help", argc, argv, i ) )
        {
            printUsage();
            exit( 0 );
        }
        else if ( isNextArgument   ( "-a", "--inputA"              , argc, argv, i, &args.argInputA                          ) ) {}
        else if ( isNextArgument   ( "-b", "--inputB"              , argc, argv, i, &args.argInputB                          ) ) {}
        else if ( isNextArgument   ( "-o", "--output"              , argc, argv, i, &args.argOutput                          ) ) {}
        else if ( parseNextArgument( "-m", "--mode"                , argc, argv, i, params, COMPARE_OPTION_MODE              ) ) {}
        else if ( isNextArgumentInt( "-k", "--max-length"          , argc, argv, i, &args.argMaxLength                          ) ) {}
        else if ( isNextArgumentInt( "-n", "--min-occ"             , argc, argv, i, &args.argMinOcc                          ) ) {}

        else if ( isNextArgument   ( "-c", "--genome-metadata"     , argc, argv, i, &args.argInputC                          ) ) {}
        else if ( isNextArgument   ( "-t", "--taxonomy"            , argc, argv, i, &args.argTaxonomy                        ) ) {}
        else if ( isNextArgumentInt( "-w", "--min-kmer-length"     , argc, argv, i, &args.argMinKmerLength                   ) ) {}
        else if ( parseNextArgument( "-d", "--report-min-length"   , argc, argv, i, params, COMPARE_OPTION_REPORT_MINLENGTH  ) ) {}
        else if ( isNextArgument   ( ""  , "--subset"              , argc, argv, i, &args.argSubset                          ) ) {}
        else if ( isNextArgument   ( "-T", "--temp-directory"      , argc, argv, i, &args.argTempPath                        ) ) {}
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
    if ( args.argInputA.empty()
         || args.argInputB.empty()
         || params.getValue( COMPARE_OPTION_MODE ) == MULTIPLE_OPTIONS
       )
    {
        cerr << "Error: Missing arguments: All parameters below are required\n" << endl;
        printUsage();
        exit( 1 );
    }

    if ( params.getValue( COMPARE_OPTION_MODE ) == MODE_METAGENOMICS
         && (
             args.argTaxonomy.empty()
         )
       )
    {
        cerr << "Error: Missing Metagenomics-specific parameters:\n" << endl;
        printUsage();
        exit( 1 );
    }

    // Launch
    launchBeetlCompare( params, args );

    return 0;
}
