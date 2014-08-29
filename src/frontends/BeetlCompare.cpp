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

#include "BeetlCompare.hh"

#include "Common.hh"
#include "config.h"
#include "countWords/CountWords.hh"
#include "parameters/CompareParameters.hh"
#include "libzoo/cli/Common.hh"
#include "libzoo/util/Logger.hh"
#include "libzoo/util/TemporaryFilesManager.hh"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string.h>
#include <unistd.h>

using namespace std;
using namespace BeetlCompareParameters;


CompareParameters params;

void printUsage()
{
    params.printUsage();

    cout << "Notes:" << endl;
    cout << "    Mode = tumour-normal: Set A = tumour, set B = normal" << endl;
    cout << "    Mode = splice       : Set A = DNA; set B = RNA" << endl;
    cout << "    Mode = reference    : Set B is a reference genome" << endl;
    cout << "    Mode = metagenomics : Set B are merged reference genomes. Switches metagenome classifier on." << endl;
    cout << endl;
}

void launchBeetlCompare()
{
    vector<string> setA_filenames;
    vector<string> setB_filenames;
    vector<string> setB_C_filenames;
    bool setA_isBwtCompressed;
    bool setB_isBwtCompressed;
    string setA_availableFileLetters;
    string setB_availableFileLetters;
    detectInputBwtProperties( params["input setA"], setA_filenames, setA_isBwtCompressed, setA_availableFileLetters );
    detectInputBwtProperties( params["input setB"], setB_filenames, setB_isBwtCompressed, setB_availableFileLetters );

    if ( setA_filenames.size() < 2 )
    {
        cerr << "Error: too few input files detected (run with -vv for more details)" << endl;
        exit( -1 );
    }

    assert( setA_filenames.size() == setB_filenames.size() );
    if ( params["mode"] == MODE_METAGENOMICS )
    {
        assert( strchr( setB_availableFileLetters.c_str(), 'C' ) && "{inputB}-C0x files cannot be found" );
        setB_C_filenames = setB_filenames;
        for ( unsigned int i = 0; i < setB_C_filenames.size(); ++i )
        {
            int pos = setB_C_filenames[i].size() - 3;
            assert( setB_C_filenames[i][pos] == 'B' );
            setB_C_filenames[i][pos] = 'C';
        }

        // Check that the 'C' files are encoded as ints and not anymore as shorts,
        // by checking that the first few 32-bit non-zero values are <65535 (unsigned short encoding would lead to large ints)
        ifstream fileC0( setB_C_filenames[0].c_str() );
        uint32_t intVal;
        for ( unsigned int i = 0; i < 10; ++i )
        {
            if ( !fileC0.read( ( char * )&intVal, 4 ) )
                break;
            if ( intVal > 65535 )
            {
                cerr << "Error: inputSetB-C0x files seem to be encoded as 2-bytes values. This version of BEETL expects 4 bytes per value.\nYou can use the tool scripts/misc/shortToInt.pl provided with BEETL source code to convert your files." << endl;
                exit( -1 );
            }
        }
    }

    if ( !params["inputA format"].isSet() )
    {
        if ( setA_isBwtCompressed )
        {
            clog << "Set A detected as RLE compressed" << endl;
            params["inputA format"] = "BWT_RLE";
        }
        else
        {
            clog << "Set A detected as ASCII" << endl;
            params["inputA format"] = "BWT_ASCII";
        }
    }
    if ( !params["inputB format"].isSet() )
    {
        if ( setB_isBwtCompressed )
        {
            clog << "Set B detected as RLE compressed" << endl;
            params["inputB format"] = "BWT_RLE";
        }
        else
        {
            clog << "Set B detected as ASCII" << endl;
            params["inputB format"] = "BWT_ASCII";
        }
    }

    // Use default parameter values where needed
    params.commitDefaultValues();

    // Update variable with optinal user values
    setA_isBwtCompressed = ( params["inputA format"] == INPUT_FORMAT_BWT_RLE );
    setB_isBwtCompressed = ( params["inputB format"] == INPUT_FORMAT_BWT_RLE );

    bool reportMinLength = ( params["report min length"] == REPORT_MINLENGTH_ON );

    Logger::out() << "\nLaunching the following configuration of Beetl-compare:" << endl;
    params.print( Logger::out(), false );
    Logger::out() << endl;

    Algorithm *pcountWords = new CountWords( setA_isBwtCompressed, setB_isBwtCompressed
            , 'X'
            , params["min occ"]
            , params["max length"]
            , setA_filenames, setB_filenames
            , setB_C_filenames
            , params["taxonomy"]
            , reportMinLength
            , params["min kmer length"]
            , params["subset"]
            , &params
                                           );

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

    if ( !params.parseArgv( argc, argv ) || params["help"] == 1 || !params.chechRequiredParameters() )
    {
        printUsage();
        exit( params["help"] == 0 );
    }

    // Checking for extra parameters required in metagenomics mode
    if ( params["mode"] == MODE_METAGENOMICS && !params["taxonomy"].isSet() )
    {
        cerr << "Error: Missing Metagenomics-specific parameter: --taxonomy\n" << endl;
        printUsage();
        exit( 1 );
    }

    // Auto-detection of missing arguments
    if ( !params["memory limit MB"].isSet() )
    {
        params["memory limit MB"] = detectMemoryLimitInMB();
    }
    TemporaryFilesManager::get().setRamLimit( params["memory limit MB"] );

    // Launch
    launchBeetlCompare();

    return 0;
}
