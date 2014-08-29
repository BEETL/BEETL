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

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unistd.h>

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
    else if ( endsWith( inputFilename, ".fq" ) )
    {
        return fileFormatLabels[FILE_FORMAT_FASTQ];
    }
    else if ( endsWith( inputFilename, ".seq" ) )
    {
        return fileFormatLabels[FILE_FORMAT_SEQ];
    }
    else if ( endsWith( inputFilename, ".bcl" ) )
    {
        return fileFormatLabels[FILE_FORMAT_BCL];
    }
    else if ( endsWith( inputFilename, ".cyc" ) || endsWith( inputFilename, "cyc." ) || beginsWith( inputFilename, "cyc." ) )
    {
        return fileFormatLabels[FILE_FORMAT_CYC];
    }
    return "";
}

void checkFileFormat( const string &inputFilename, const string &fileFormat )
{
    if ( fileFormat.empty() )
    {
        cerr << "Error: filename suffix not recognised for " << inputFilename << ". (You may try to run with --help, as we usually provide options such as --input-format to specify the file format)" << endl;
        exit( EXIT_FAILURE );
    }
}

bool isNextArgument( const string &shortPrefix, const string &longPrefix, const int argc, const char **argv, int &i, string *argValue )
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

bool isNextArgumentInt( const string &shortPrefix, const string &longPrefix, const int argc, const char **argv, int &i, int *argValue )
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

bool parseNextArgument( const string &shortPrefix, const string &longPrefix, const int argc, const char **argv, int &i, ToolParameters &toolParams, const unsigned toolParamKey )
{
    string argValue;
    if ( isNextArgument( shortPrefix, longPrefix, argc, argv, i, &argValue ) )
    {
        toolParams[toolParamKey] = argValue;
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

    // Use the same path as the current executable if we manage to extract it
    // Note: "/proc/self/exe" only exists on linux, but launching executables is a temporary solution anyway
    char selfPath[1024];
    ssize_t len = ::readlink( "/proc/self/exe", selfPath, sizeof( selfPath ) - 1 );
    if ( len != -1 )
    {
        selfPath[len] = '\0';
        char *lastSlash = strrchr( selfPath, '/' );
        if ( lastSlash )
        {
            *( lastSlash + 1 ) = 0;
            path << selfPath;
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

bool doesFileExist( const string &filename )
{
    return access( filename.c_str(), F_OK ) == 0;
}
