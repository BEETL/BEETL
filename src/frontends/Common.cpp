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

#include "Common.hh"

#include "config.h"
#include "libzoo/util/Logger.hh"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifdef HAVE_SYSINFO
#include <sys/sysinfo.h>
#endif //ifdef HAVE_SYSINFO

#ifdef HAVE_SYSCTLBYNAME
#include <sys/types.h>
#include <sys/sysctl.h>
#endif //ifdef HAVE_SYSCTLBYNAME

#ifdef HAVE_PROC_PIDPATH
#include "libproc.h"
#endif //ifdef HAVE_PROC_PIDPATH

using namespace std;


int detectMemoryLimitInMB()
{
    // Strategy: Retrieve the smallest of ulimit -v and /proc/meminfo
    int result = 0;

    // ulimit -v
    struct rlimit rl;
    if ( getrlimit( RLIMIT_AS, &rl ) != -1 )
    {
        if ( RLIM_INFINITY != rl.rlim_cur )
        {
            result = ( rl.rlim_cur >> 20 );
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "ulimit -v: " << result << " MB" << endl;
        }
    }

    // /proc/meminfo
#ifdef HAVE_SYSINFO
    struct sysinfo info;
    if ( sysinfo( &info ) != -1 )
    {
        const int totalRamInMB = info.totalram >> 20;
        if ( result == 0 )
            result = totalRamInMB;
        else
            result = min( result, totalRamInMB );
    }
#endif //ifdef HAVE_SYSINFO

#ifdef HAVE_SYSCTLBYNAME
    const char *name = "hw.memsize";
    unsigned long long value = 0;
    size_t valueLength = sizeof( value );
    if ( sysctlbyname( name, &value, &valueLength, NULL, 0 ) != -1 )
    {
        int totalRamInMB = 0;
        if ( valueLength == sizeof( unsigned long long ) )
            totalRamInMB = value >> 20;
        else if ( valueLength == sizeof( unsigned int ) )
            totalRamInMB = ( *( ( unsigned int * )&value ) ) >> 20;

        if ( totalRamInMB )
        {
            if ( result == 0 )
                result = totalRamInMB;
            else
                result = min( result, totalRamInMB );
        }
    }
#endif //ifdef HAVE_SYSCTLBYNAME

    // Use a default minimum in case all detections failed
    if ( result == 0 )
    {
        Logger::out() << "Warning: Couldn't determine RAM constraint. Using default 1024 MB, you can specify it with --memory-limit." << endl;
        result = 1024; // 1 GB
    }

    return result;
}
