/*
* Collection of basic tools and defines
*/

#ifndef _TOOLS_H_
#define _TOOLS_H_

#include <iostream>
#include <fstream>

#define TERMINATE_CHAR '$'

#ifndef uchar
#define uchar unsigned char
#endif
#ifndef uint
#define uint unsigned int
#endif
#ifndef ulong
#define ulong unsigned long
#endif

#define dataTypedimAlpha uchar //size of the alphabet (in biologic case 6 ($,A,C,G,N,T))
#define dataTypelenSeq uint //length of the sequences (in biologic case 100)

#if dataTypeNumSeq == 1
# define dataTypeNSeq ulong
#else
# define dataTypeNSeq uint
#endif

#if dataTypeNumChar == 1
# define dataTypeNChar ulong
#else
# define dataTypeNChar uint
#endif

//They go in Config.h
//#define dataTypeNumSeq 2 //number of sequences in the input file
//#define dataTypeNumChar 1 //numer of characters in the input file (length of the BWT)
//#define verboseEncode 0
//#define verboseDecode 1
//#define deletePartialBWT 0
//#define BackByVector 1

class Tools
{
public:
static uchar * GetRandomString(unsigned, unsigned, unsigned &);
static uchar * GetFileContents(char *, ulong =0);
static unsigned FloorLog2(ulong);
static unsigned CeilLog2(ulong);
static unsigned* MakeTable();
static unsigned FastFloorLog2(unsigned);
};

#endif
