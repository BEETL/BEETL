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
 **
 ** Collection of basic tools and defines
 **/

#ifndef _TOOLS_H_
#define _TOOLS_H_

#include "Types.hh"

#include <cassert>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using std::string;
using std::vector;


#define verboseEncode 0
#define verboseDecode 0
#define deletePartialBWT 0 //If it is set to 1, it deletes the BWT-segments files and keeps the entire BWT, otherwise renames them.
#define deletePartialLCP 0	//If it is set to 1, it deletes the LCP-segments files and keeps the entire LCP, otherwise renames them.
//#define deletePartialSA 0 //If it is set to 1, it deletes the SA-segments files and keeps the entire BWT, otherwise renames them.
#define deleteCycFile 1  //If it is set to 1, it deletes the cycs files.
#define BUILD_SA 0   //If it is set to 1, it computes the GSA (seqID, position) and the SA (position of the concatenated sequences without a further end-marker).
#define BUILD_LCP 1	    	//If it is set to 1, it uses larger structures to allow --generate-lcp

class Tools
{
private:
    static double startTime;
public:
    static void StartTimer();
    static double GetTime();
    static uchar *GetRandomString( unsigned, unsigned, unsigned & );
    static uchar *GetFileContents( char *, size_t = 0 );
    static uint8_t FloorLog2( uint64_t );
    static uint8_t CeilLog2( uint64_t );
    static uint8_t *MakeTable();
    static uint8_t FastFloorLog2( uint32_t );
};

void getFileName( const string &stem, const char code, const int pile,
                  string &fileName );

bool isValidFastaFile( const char *filename );

bool isValidReadFile( const char *filename );

bool readWriteCheck( const char *fileName, const bool readWrite, const bool failIfError = true );

void checkIfEqual( const int arg1, const int arg2 );

void checkIfNotEqual( const int arg1, const int arg2 );

bool hasPrefix( const string &fullString, const string &prefix );
bool hasSuffix( const string &fullString, const string &suffix );
vector<string> splitString ( string s, const string &token );

void detectInputBwtProperties( const string &prefix, vector<string> &filenames, bool &isBwtCompressed, string &availableFileLetters );

int safeRename( const string &from, const string &to );
void pauseBetweenCycles();
void readProcSelfStat( int &out_pid, int &out_num_threads, int &out_processor );

#endif
