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

#ifndef LOGGER_HH
#define LOGGER_HH

#include <iostream>
#include <string>

using std::string;


enum LOG_LEVEL { LOG_LEVEL_QUIET = 0, LOG_LEVEL_NORMAL, LOG_SHOW_IF_VERBOSE, LOG_SHOW_IF_VERY_VERBOSE, LOG_FOR_DEBUGGING };

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
