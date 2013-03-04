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

#ifndef TRANPOSEFASTA_INCLUDED
#define TRANPOSEFASTA_INCLUDED

#include "Tools.hh"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using std::string;
using std::vector;


#define BUFFERSIZE 1024 // 2^20
//#define CYCLENUM 100


class SeqReaderFile;

class TransposeFasta
{
public:
    TransposeFasta();
    void init( SeqReaderFile *pReader, const bool processQualities = true );
    ~TransposeFasta();

    bool convert( const string &input, const string &output );   //Input from Fasta file (converts Fasta File into cyc Files)
    bool inputCycFile( const string &cycPrefix );                                    //Input from cyc files
    bool convertFromCycFileToFastaOrFastq( const string &fileInputPrefix, const string &fileOutput );      //Convert cyc files into Fasta or Fastq File
    bool hasProcessedQualities()
    {
        return processQualities_;
    }

    dataTypelenSeq lengthRead;    //Lenght of each text
    dataTypeNChar lengthTexts;   //Total length of all texts without $-symbols

    dataTypeNSeq nSeq;   //number total of texts in filename1
    dataTypeNChar freq[256];  //contains the distribution of the symbols. It is useful only for testing. It depends on the #characters

private:
    SeqReaderFile *pReader_;
    uint cycleNum_;
    vector<FILE *> outputFiles_;
    vector<vector<uchar> > buf_;
    //    FILE* outputFiles_[CYCLENUM];

    //    uchar buf_[CYCLENUM][BUFFERSIZE];
    bool processQualities_;
};

#endif
