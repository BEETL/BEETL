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

#include "Logger.hh"

using namespace std;


NullStream nullStream;
int Logger::currentVerbosity = 0;


void Logger::setVerbosity( const string &verbosityString )
{
    if ( verbosityString == "quiet" ||
         verbosityString == "0" )
        currentVerbosity = 0;

    else if ( verbosityString == "verbose" ||
              verbosityString == "1" )
        currentVerbosity = 1;

    else if ( verbosityString == "very-verbose" ||
              verbosityString == "2" )
        currentVerbosity = 2;

    else if ( verbosityString == "debug" ||
              verbosityString == "3" )
        currentVerbosity = 3;

    else
    {
        clog << "Warning: Invalid verbosity value. Setting to \"verbose\"" << endl;
        currentVerbosity = 1;
    }

    Logger::out( LOG_SHOW_IF_VERBOSE ) << "Setting logging to level " << currentVerbosity << endl;
}
