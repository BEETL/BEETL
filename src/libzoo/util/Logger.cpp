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

#include "Logger.hh"

using namespace std;


NullStream nullStream;
int Logger::currentVerbosity = 1;


void Logger::setVerbosity( const string &verbosityString )
{
    if ( verbosityString == "quiet" ||
         verbosityString == "0" )
        currentVerbosity = 0;

    else if ( verbosityString == "normal" ||
              verbosityString == "1" )
        currentVerbosity = 1;

    else if ( verbosityString == "verbose" ||
              verbosityString == "2" )
        currentVerbosity = 2;

    else if ( verbosityString == "very-verbose" ||
              verbosityString == "3" )
        currentVerbosity = 3;

    else if ( verbosityString == "debug" ||
              verbosityString == "4" )
        currentVerbosity = 4;

    else
    {
        clog << "Warning: Invalid verbosity value. Setting to \"normal\"" << endl;
        currentVerbosity = 1;
    }

    Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Setting logging to level " << currentVerbosity << endl;
}
