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

#ifndef INCLUDED_TIMER
#define INCLUDED_TIMER

#include <ctime>
#include <fstream>
#include <sys/resource.h>
#include <sys/time.h>


// Class Name : Timer
// Description: Maintains info on actual and processing time
class Timer
{

public:
    Timer( void );
    std::ostream &print( std::ostream &os );
    // timeNow: returns current date and time as an ASCII string
    const char *timeNow( void ) const;

private:
    rusage thisUsage_;
    rusage lastUsage_;
    timeval thisTime_;
    timeval lastTime_;
    //  timeb thisTime_;
    //  timeb lastTime_;

}; // Timer

std::ostream &operator<<( std::ostream &os, Timer &timer );

#endif
// end of Timer.hh
