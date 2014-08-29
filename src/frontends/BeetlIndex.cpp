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

#include "BeetlIndex.hh"

#include "BwtIndex.hh"
#include "BwtReader.hh"
#include "Tools.hh"
#include "config.h"
//#include "search/SearchUsingBacktracker.hh"
#include "parameters/IndexParameters.hh"
#include "libzoo/cli/Common.hh"
#include "libzoo/util/Logger.hh"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <unistd.h>

using namespace std;


IndexParameters params;

void printUsage()
{
    params.printUsage();

    cout << "Notes:" << endl;
    cout << "    -j and -k are mutually exclusive, one of them being required.\n" << endl;
    cout << endl;
}

void launchBeetlIndex()
{
    //    SearchUsingBacktracker search( params );
    //    search.run();


    string indexPrefix = params.getStringValue( "input" );
    bool compressed, forceOverwrite;
    vector<string> pileNames;
    string dummyStr;

    int blockSize( 2048 );
    if ( params["block size"].isSet() )
    {
        blockSize = params["block size"];
        if ( blockSize <= 0 )
        {
            cerr << "Index block size must be greater than zero!" << endl;
            exit( EXIT_FAILURE );
        } // ~if
    } // ~if
    else
    {
        cerr << "No index block size specified, using default" << endl;
    } // ~else
    cerr << "Will index every " << blockSize << " bytes " << endl;

    if ( params["force"].isSet() )
    {
        forceOverwrite = true;
    } // ~else


    detectInputBwtProperties( indexPrefix, pileNames, compressed, dummyStr );

    if ( pileNames.empty() )
    {
        cerr << "Did not find any BWT files matching prefix " << indexPrefix
             << "." << endl
             << "If BWT files are named 'myName-B0?' specify '-i myname'."
             << endl;
        exit ( EXIT_FAILURE );
    }

    //    if(!compressed)
    //    {
    //        cerr<<"BWT files seems to be in ASCII format." << endl
    //            <<"Only run-length compressed BWT files can currently be indexed."
    //            << endl <<"Use beetl-convert to change format." << endl;
    //        exit( EXIT_FAILURE );
    //    }


#pragma omp parallel for
    for ( unsigned int thisPile = 0; thisPile < pileNames.size(); ++thisPile )
    {
        const string pileName = pileNames[thisPile];
#pragma omp critical (IO)
        cerr << "Indexing file " << pileName << endl;
//        BwtReaderIndex<BwtReaderRunLength> reader( pileName.c_str(), params.getStringValue( "use shm" ) );
        unique_ptr<BwtReaderBase> reader( instantiateBwtPileReader( pileName.c_str(), params.getStringValue( "use shm" ), false, true ) );
        string indexFileName = pileName + ".idx";


        FILE *pFile;
        if ( !forceOverwrite && readWriteCheck( indexFileName.c_str(), false, false ) )
        {
#pragma omp critical (IO)
            cerr << "File " << indexFileName << " already exists! Rerun with --force to remove." << endl;
            exit( EXIT_FAILURE );
            fclose( pFile );
        }

        pFile = fopen( indexFileName.c_str() , "w" );

        if ( pFile == NULL )
        {
#pragma omp critical (IO)
            cerr << "Problem opening file " << indexFileName << " for writing" << endl;
            exit( EXIT_FAILURE );
        }
        buildIndex( reader.get(), pFile, blockSize );
        fclose ( pFile );
    }
}

int main( const int argc, const char **argv )
{
    // Generated using: http://patorjk.com/software/taag/#p=display&f=Soft&t=BEETL%20index
    cout << ",-----.  ,------.,------.,--------.,--.       ,--.           ,--.                  " << endl;
    cout << "|  |) /_ |  .---'|  .---''--.  .--'|  |       `--',--,--,  ,-|  | ,---. ,--.  ,--. " << endl;
    cout << "|  .-.  \\|  `--, |  `--,    |  |   |  |       ,--.|      \\' .-. || .-. : \\  `'  /  " << endl;
    cout << "|  '--' /|  `---.|  `---.   |  |   |  '--.    |  ||  ||  |\\ `-' |\\   --. /  /.  \\  " << endl;
    cout << "`------' `------'`------'   `--'   `-----'    `--'`--''--' `---'  `----''--'  '--' " << endl;
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

    // Checking for required parameters - TBD set to default (prob 2048) if not specified
    //    if ( ! ( params["kmers input file"].isSet() ) )
    //    {
    //        cerr << "Error: Missing argument: -i is required\n" << endl;
    //        printUsage();
    //        exit( 1 );
    //    }

    // Launch
    launchBeetlIndex();

    return 0;
}
