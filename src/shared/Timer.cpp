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

#include "Timer.hh"

#include <cstdlib>


Timer::Timer( void )
{
    if ( gettimeofday( &lastTime_, NULL ) != 0 )
        exit( EXIT_FAILURE );
    if ( getrusage( RUSAGE_SELF, &lastUsage_ ) != 0 )
        exit( EXIT_FAILURE );
} // ~Timer::Timer( void )


std::ostream &Timer::print( std::ostream &os )
{
    if ( gettimeofday( &thisTime_, NULL ) != 0 )
        exit( EXIT_FAILURE );
    if ( getrusage( RUSAGE_SELF, &thisUsage_ ) != 0 )
        exit( EXIT_FAILURE );

    static double elapsedActual, elapsedUser, elapsedSystem;

    elapsedUser = thisUsage_.ru_utime.tv_sec - lastUsage_.ru_utime.tv_sec
                  + ( thisUsage_.ru_utime.tv_usec
                      - ( double ) lastUsage_.ru_utime.tv_usec ) / 1000000;

    elapsedSystem = thisUsage_.ru_stime.tv_sec - lastUsage_.ru_stime.tv_sec
                    + ( thisUsage_.ru_stime.tv_usec
                        - ( double ) lastUsage_.ru_stime.tv_usec ) / 1000000;

    elapsedActual = thisTime_.tv_sec - lastTime_.tv_sec + ( thisTime_.tv_usec
                    - ( double ) lastTime_.tv_usec ) / 1000000;

    os << "User: " << elapsedUser << "s System: " << elapsedSystem
       << "s Actual: " << elapsedActual << "s Efficiency: "
       << ( ( elapsedActual == 0 ) ? 0 : ( ( elapsedUser + elapsedSystem )
                                           * 100.0 / elapsedActual ) ) << '\045'; // only way to print % liked by both Intel and gcc!
    lastTime_ = thisTime_;
    lastUsage_ = thisUsage_;
    return os;
} // ~std::ostream& Timer::print( std::ostream& os )

std::ostream &operator<<( std::ostream &os, Timer &timer )
{
    return timer.print( os );
} // ~std::ostream& operator<<( std::ostream& os, Timer& timer )

const char *Timer::timeNow( void ) const
{
    time_t tt( time( NULL ) );
    return ctime( &tt );
} // ~const char* Timer::timeNow( void ) const
