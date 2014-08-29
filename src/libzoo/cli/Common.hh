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

#ifndef LIBZOO_CLI_COMMON_HH
#define LIBZOO_CLI_COMMON_HH

#include "ToolParameters.hh"

#include <string>

using std::string;


// options: file format

enum FileFormat
{
    FILE_FORMAT_FASTA,
    FILE_FORMAT_FASTQ,
    FILE_FORMAT_SEQ,
    FILE_FORMAT_CYC,
    FILE_FORMAT_BCL,
    FILE_FORMAT_COUNT
};

const string fileFormatLabels[] =
{
    "fasta",
    "fastq",
    "seq",
    "cyc",
    "bcl",
    "" // end marker
};


bool beginsWith( const string &s, const string &prefix );
bool endsWith( const string &s, const string &suffix );
string detectFileFormat( const string &inputFilename );
void checkFileFormat( const string &inputFilename, const string &fileFormat );
bool isNextArgument( const string &shortPrefix, const string &longPrefix, const int argc, const char **argv, int &i, string *argValue );
bool isNextArgumentInt( const string &shortPrefix, const string &longPrefix, const int argc, const char **argv, int &i, int *argValue );
bool parseNextArgument( const string &shortPrefix, const string &longPrefix, const int argc, const char **argv, int &i, ToolParameters &, const unsigned toolParamKey );
void launchBeetl( const string &params );
bool doesFileExist( const string &filename );

#endif // LIBZOO_CLI_COMMON_HH
