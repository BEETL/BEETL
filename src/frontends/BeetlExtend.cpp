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

#include "libzoo/util/ColorText.hh"

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
        exit( params["help"] == 0 );
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
