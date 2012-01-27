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

#include <iostream>
#include <fstream>
#include <string>
//#include <cstdlib>
//#include <stdlib.h>
#include <assert.h>
#include <ctime>
#include "Types.hh"

using namespace std;

#define TERMINATE_CHAR '$'

#ifdef THESE_DEFS_REPEATED_IN_TYPES_DOT_HH
#ifndef uchar
#define uchar unsigned char
#endif
#ifndef uint
#define uint unsigned int
#endif
#ifndef ulong
#define ulong unsigned long
#endif
#endif

#define dataTypedimAlpha uchar  //size of the alphabet (in biologic case 6 ($,A,C,G,N,T))
#define dataTypelenSeq uint	//length of the sequences (in biologic case 100)
#define dataTypeNumSeq 0		//number of sequences in the input file
#define dataTypeNumChar 1		//numer of characters in the input file (length of the BWT)

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

#define verboseEncode 0
#define verboseDecode 0
#define deletePartialBWT 0
#define BackByVector 1

class Tools
{
private:
    static double startTime;
public:
    static void StartTimer();
    static double GetTime();
    static uchar * GetRandomString(unsigned, unsigned, unsigned &);
 	static uchar * GetFileContents(char *, ulong =0);
    static unsigned FloorLog2(ulong);
    static unsigned CeilLog2(ulong);
    static unsigned* MakeTable();
    static unsigned FastFloorLog2(unsigned);
};

void getFileName( const string& stem, const char code, const int pile,
		   string& fileName );

bool isValidFastaFile(const char *filename);

bool isValidReadFile(const char *filename);

void readWriteCheck(const char *fileName, const bool readWrite);

void checkIfEqual(const int arg1, const int arg2);

void checkIfNotEqual(const int arg1, const int arg2);
#endif

