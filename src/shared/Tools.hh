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


#define TERMINATE_CHAR '$'

#define dataTypedimAlpha uchar  //size of the alphabet (in biologic case 6 ($,A,C,G,N,T))
#define dataTypelenSeq uint //length of the sequences (in biologic case 100)
#define dataTypeNumSeq 0  //number of sequences in the input file
#define dataTypeNumChar 1  //numer of characters in the input file (length of the BWT)

// USE_ATTRIBUTE_PACKED: if set, set __attribute__((packed)) for the
// struct sortElement in Sorting.hh; reduces the memory consumption from
// 24 to 17 bytes per input base
#define USE_ATTRIBUTE_PACKED 1


typedef LetterCountType dataTypeNChar;
typedef SequenceNumberType dataTypeNSeq;

#ifdef BCR_AND_BCREXT_NOW_USE_SAME_TYPES

#if dataTypeNumSeq == 1
#   define dataTypeNSeq ulong
#else
#   define dataTypeNSeq uint
#endif

#if dataTypeNumChar == 1
#   define dataTypeNChar ulong
#else
#   define dataTypeNChar uint
#endif

#endif

//It is the definition of each element of the generalized suffix array (GSA)
struct ElementType
{
    dataTypelenSeq sa;          //It is the position in the sequence, so it goes from 0 a length read
    dataTypeNSeq numSeq;  //It is the number of the sequence.
};

#define verboseEncode 0
#define verboseDecode 0
#define deletePartialBWT 0 //If it is set to 1, it deletes the BWT-segments files and keeps the entire BWT, otherwise renames them.
#define deletePartialLCP 0	//If it is set to 1, it deletes the LCP-segments files and keeps the entire LCP, otherwise renames them.
//#define deletePartialSA 0 //If it is set to 1, it deletes the SA-segments files and keeps the entire BWT, otherwise renames them.
#define deleteCycFile 1  //If it is set to 1, it deletes the cycs files.
#define BUILD_SA 0   //If it is set to 1, it computes the GSA (seqID, position) and the SA (position of the concatenated sequences without a further end-marker).
#define BUILD_LCP 1	    	//If it is set to 1, it uses larger structures to allow --generate-lcp
#define BackByVector 1  //if it is set to 1, it uses the sampling of the BWT segments for inverse BWT. More memory, less time.

class Tools
{
private:
    static double startTime;
public:
    static void StartTimer();
    static double GetTime();
    static uchar *GetRandomString( unsigned, unsigned, unsigned & );
    static uchar *GetFileContents( char *, ulong = 0 );
    static unsigned FloorLog2( ulong );
    static unsigned CeilLog2( ulong );
    static unsigned *MakeTable();
    static unsigned FastFloorLog2( unsigned );
};

void getFileName( const string &stem, const char code, const int pile,
                  string &fileName );

bool isValidFastaFile( const char *filename );

bool isValidReadFile( const char *filename );

void readWriteCheck( const char *fileName, const bool readWrite );

void checkIfEqual( const int arg1, const int arg2 );

void checkIfNotEqual( const int arg1, const int arg2 );

bool hasPrefix( const string &fullString, const string &prefix );
bool hasSuffix( const string &fullString, const string &suffix );

#endif
