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

#ifndef LOGGER_HH
#define LOGGER_HH

#include <iostream>
#include <string>

using std::string;


enum LOG_LEVEL { LOG_LEVEL_QUIET = 0, LOG_ALWAYS_SHOW, LOG_SHOW_IF_VERBOSE, LOG_SHOW_IF_VERY_VERBOSE, LOG_FOR_DEBUGGING };

class NullStream : public std::ostream
{
public:
    NullStream() : std::ios( 0 ), std::ostream( 0 ) {}
};

extern NullStream nullStream;

#define Logger_if(verbosity) if (Logger::currentVerbosity >= verbosity)

class Logger
{
public:
    static int currentVerbosity;

    static inline std::ostream &out()
    {
        return std::cout;
    }

    /*
        static std::ostream &out( const int verbosity )
        {
            if ( currentVerbosity >= verbosity )
            {
                return std::cout;
            }
            else
            {
                return nullStream;
            }
        }
    */

    static std::ostream &error( const int verbosity = 0 )
    {
        return std::cerr;
    }

    static void setVerbosity( const string &verbosityString );
};

#endif // LOGGER_HH
