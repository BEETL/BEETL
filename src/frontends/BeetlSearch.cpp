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

#include "config.h"
#include "search/SearchUsingBacktracker.hh"
#include "parameters/SearchParameters.hh"

SearchParameters params;


void printUsage()
{
    params.printUsage();

    cout << "Notes:" << endl;
    cout << "    -j and -k are mutually exclusive, one of them being required.\n" << endl;
    cout << endl;
}

void launchBeetlSearch()
{
    SearchUsingBacktracker search( params );
    search.run();
}

int main( const int argc, const char **argv )
{
    // Generated using: http://patorjk.com/software/taag/#p=display&f=Soft&t=BEETL%20search
    clog << ",-----.  ,------.,------.,--------.,--.                                          ,--.      " << endl;
    clog << "|  |) /_ |  .---'|  .---''--.  .--'|  |        ,---.  ,---.  ,--,--.,--.--. ,---.|  ,---.  " << endl;
    clog << "|  .-.  \\|  `--, |  `--,    |  |   |  |       (  .-' | .-. :' ,-.  ||  .--'| .--'|  .-.  | " << endl;
    clog << "|  '--' /|  `---.|  `---.   |  |   |  '--.    .-'  `)\\   --.\\ '-'  ||  |   \\ `--.|  | |  | " << endl;
    clog << "`------' `------'`------'   `--'   `-----'    `----'  `----' `--`--'`--'    `---'`--' `--' " << endl;
    clog << "Version " << PACKAGE_VERSION << endl;
    clog << endl;

    clog << "Command called:" << endl << "   ";
    for ( int i = 0; i < argc; ++i )
    {
        clog << " " << argv[i];
    }
    clog << "\n" << endl;

    if ( !params.parseArgv( argc, argv ) || params["help"] == 1 || !params.chechRequiredParameters() )
    {
        printUsage();
        exit( params["help"] == 0 );
    }

    // Use default parameter values where needed
    params.commitDefaultValues();

    // Checking for required parameters
    if ( ! ( params["kmers input file"].isSet() ^ params["one kmer string"].isSet() ) )
    {
        cerr << "Error: Missing or incorrect arguments: -i is required; -j and -k are mutually exclusive, one of them being required\n" << endl;
        printUsage();
        exit( 1 );
    }

    // Launch
    launchBeetlSearch();

    return 0;
}
