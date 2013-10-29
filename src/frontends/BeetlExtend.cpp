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

#include "config.h"
#include "parameters/ExtendParameters.hh"
#include "search/Extender.hh"

using namespace std;

ExtendParameters params;


void printUsage()
{
    params.printUsage();
}

void launchBeetlExtend()
{
    Extender extender( params );
    extender.run();
}

int main( const int argc, const char **argv )
{
    // Generated using: http://patorjk.com/software/taag/#p=display&f=Soft&t=BEETL%20extend
    cout << ",-----.  ,------.,------.,--------.,--.                          ,--.                    ,--. " << endl;
    cout << "|  |) /_ |  .---'|  .---''--.  .--'|  |        ,---. ,--.  ,--.,-'  '-. ,---. ,--,--,  ,-|  | " << endl;
    cout << "|  .-.  \\|  `--, |  `--,    |  |   |  |       | .-. : \\  `'  / '-.  .-'| .-. :|      \\' .-. | " << endl;
    cout << "|  '--' /|  `---.|  `---.   |  |   |  '--.    \\   --. /  /.  \\   |  |  \\   --.|  ||  |\\ `-' | " << endl;
    cout << "`------' `------'`------'   `--'   `-----'     `----''--'  '--'  `--'   `----'`--''--' `---'" << endl;
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
        exit( 1 );
    }

    // Use default parameter values where needed
    params.commitDefaultValues();

    // Checking for required parameters
    if ( ! ( params["sequence numbers output filename"].isSet() || params["dollar positions output filename"].isSet() ) )
    {
        cerr << "Error: Missing or incorrect arguments: at least one output (-o or -p) is required\n" << endl;
        printUsage();
        exit( 1 );
    }

    // Launch
    launchBeetlExtend();

    return 0;
}
