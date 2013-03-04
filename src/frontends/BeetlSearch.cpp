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

#include "Common.hh"
#include "config.h"

#include <iostream>
#include <cstdlib>

using namespace std;


void printUsage()
{
    cout << "Parameters:" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "    --help (-h)              Help" << endl;
    cout << endl;
    cout << "Notes:" << endl;
    cout << "    This interface is not yet available. Please use:" << endl;
    cout << "      OldBeetl bcr -m 2" << endl;
    cout << endl;
}

struct BeetlSearchArguments
{
};

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
    for ( int i=0; i < argc; ++i )
    {
        cout << " " << argv[i];
    }
    cout << "\n" << endl;

    BeetlSearchArguments args;

    for ( int i = 1; i < argc; ++i )
    {
        if ( isNextArgument( "-h", "--help", argc, argv, i ) )
        {
            printUsage();
            exit( 0 );
        }
        else
        {
            cerr << "Error: Invalid parameter: " << argv[i] << "\n" << endl;
            printUsage();
            exit( 1 );
        }
    }

    printUsage();
    return 0;
}
