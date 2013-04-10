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

#include "Common.hh"

#include "Logger.hh"
#include "config.h"

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


bool beginsWith( const string &s, const string &prefix )
{
    return s.compare( 0, prefix.size(), prefix ) == 0;
}

bool endsWith( const string &s, const string &suffix )
{
    return s.size() >= suffix.size()
           && s.compare( s.size() - suffix.size(), suffix.size(), suffix ) == 0;
}

string detectFileFormat( const string &inputFilename )
{
    if ( endsWith( inputFilename, ".fa" ) || endsWith( inputFilename, ".fasta" ) )
    {
        return fileFormatLabels[FILE_FORMAT_FASTA];
    }
    else if ( endsWith( inputFilename, ".fastq" ) )
    {
        return fileFormatLabels[FILE_FORMAT_FASTQ];
    }
    else if ( endsWith( inputFilename, ".cyc" ) )
    {
        return fileFormatLabels[FILE_FORMAT_CYC];
    }
    else if ( endsWith( inputFilename, ".seq" ) )
    {
        return fileFormatLabels[FILE_FORMAT_SEQ];
    }
    return "";
}

void checkFileFormat( const string &inputFilename, const string &fileFormat )
{
    if ( fileFormat.empty() )
    {
        cerr << "Error: file format not recognised for " << inputFilename << endl;
        exit( EXIT_FAILURE );
    }
}

bool isNextArgument( const string shortPrefix, const string longPrefix, const int argc, const char **argv, int &i, string *argValue )
{
    string arg( argv[i] );
    if ( arg == shortPrefix || arg == longPrefix )
    {
        if ( argValue )
        {
            if ( argc <= i + 1 )
            {
                cerr << "Error: Too few arguments after " << arg << endl;
                exit( 1 );
            }
            *argValue = string( argv[++i] );
        }
        return true;
    }
    if ( argValue
         && beginsWith( arg, longPrefix )
         && arg.size() > longPrefix.size()
         && arg[longPrefix.size()] == '='
         && !longPrefix.empty()
       )
    {
        *argValue = arg.substr( longPrefix.size() + 1 );
        return true;
    }
    return false;
}

bool isNextArgumentInt( const string shortPrefix, const string longPrefix, const int argc, const char **argv, int &i, int *argValue )
{
    string s;
    if ( isNextArgument( shortPrefix, longPrefix, argc, argv, i, &s ) )
    {
        stringstream ss;
        ss << s;
        ss >> *argValue;
        return true;
    }
    return false;
}

bool parseNextArgument( const string shortPrefix, const string longPrefix, const int argc, const char **argv, int &i, ToolParameters &toolParams, const unsigned toolParamKey )
{
    string argValue;
    if ( isNextArgument( shortPrefix, longPrefix, argc, argv, i, &argValue ) )
    {
        toolParams.set( toolParamKey, argValue );
        return true;
    }
    return false;
}

/*
bool string2bool( const string &str )
{
    if (str == "true") return true;
    if (str == "false") return false;
    if (str == "1") return true;
    if (str == "0") return false;

    cerr << "Error: Wrong boolean value \"" << str << "\"" << endl;
    exit( 1 );
}
*/

void launchBeetl( const string &params )
{
    ostringstream path;

#ifdef HAVE_PROC_PIDPATH

    // This should work on mac
    char pathBuf[PROC_PIDPATHINFO_MAXSIZE] = "\0";
    pid_t pid = getpid();
    ssize_t len = proc_pidpath ( pid, pathBuf, sizeof( pathBuf ) );

#else //ifdef HAVE_PROC_PIDPATH

    // Use the same path as the current executable if we manage to extract it
    // Note: "/proc/self/exe" only exists on linux, but launching executables is a temporary solution anyway
    char pathBuf[1024] = "\0";
    ssize_t len = ::readlink( "/proc/self/exe", pathBuf, sizeof( pathBuf ) - 1 );

#endif //ifdef HAVE_PROC_PIDPATH

    if ( len > 0 )
    {
        pathBuf[len] = '\0';
        char *lastSlash = strrchr( pathBuf, '/' );
        if ( lastSlash )
        {
            *( lastSlash + 1 ) = 0;
            path << pathBuf;
        }
    }

    string oldBeetl = path.str() + "OldBeetl";

    FILE *f = fopen( oldBeetl.c_str(), "rb" );
    if ( !f )
    {
        oldBeetl = path.str() + "../OldBeetl";
        f = fopen( oldBeetl.c_str(), "rb" );
    }
    fclose( f );

    stringstream command;
    command << oldBeetl << " " << params;

    cout << "Launching command:" << endl;
    cout << "  " << command.str() << endl;

    int ret = system( command.str().c_str() );
    if ( ret == -1 )
    {
        cerr << "Error in fork while launching Beetl" << endl;
        exit( 1 );
    }
    else
    {
        if ( WIFSIGNALED( ret ) || !WIFEXITED( ret ) )
        {
            exit( -1 );
        }
        int retChild = WEXITSTATUS( ret );
        if ( retChild != 0 )
        {
            exit( retChild );
        }
    }
}

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
            Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "ulimit -v: " << result << " MB" << endl;
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
        Logger::out( LOG_ALWAYS_SHOW ) << "Warning: Couldn't determine RAM constraint. Using default 1024 MB, you can specify it with --memory-limit." << endl;
        result = 1024; // 1 GB
    }

    return result;
}
