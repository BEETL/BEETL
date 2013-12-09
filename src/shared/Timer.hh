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
