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

#include "Algorithm.hh"
#include "OneBwtBackTracker.hh"
#include "Config.hh"
#include "LetterCount.hh"

#include "countWords/RangeStore.hh"
#include "Types.hh"
#include "WitnessReader.hh"
#include "BCRext.hh"
#include "BCRexternalBWT.hh"

#include "Timer.hh"
#include "config.h"

#include "libzoo/util/Logger.hh"
#include "shared/SeqReader.hh"

#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <math.h>

#ifndef INCLUDED_ERRORINFO_HH
#define INCLUDED_ERRORINFO_HH

using namespace std;

struct ErrorInfo
{
    ErrorInfo( void ):
        corrector( "" ),
        firstCycle( 0 ),
        lastCycle( 0 ),
        readEnd( 0 ),
        seqNum( -1 )
    {}

    ErrorInfo(
        int inSeqNum,
        int inReadEnd,
        int inLastCycle,
        int inFirstCycle,
        string inCorrector
    ):
        seqNum( inSeqNum ),
        readEnd( inReadEnd ),
        positionInRead( inReadEnd - inFirstCycle - 1 ),
        lastCycle( inLastCycle ),
        firstCycle( inFirstCycle ),
        corrector( inCorrector )
    {}

    void print();
    static bool SortByRead( ErrorInfo a, ErrorInfo b );
    static bool SortByReadThenLastCycle( ErrorInfo a, ErrorInfo b );

    int firstCycle; //first cycle the error is noticed
    int lastCycle; // last cycle where error is noticed
    int readEnd; //cycle at which $ for read is reached
    string corrector; // what read should be corrected to
    int seqNum; // position of read in original list (zero indexed)
    int positionInRead; //zero indexed position of error which corrector noticed within read
    int correctorStart; //zero indexed position of first character in correction string relative to the original read before alignment
    bool reverseStrand; //flag depicting which strand the correction was noticed on

    //use the -end-pos file to map the alphabetical read positions back to the original ordering
    static void SetReadNumbersToOriginal( char *endPosFileName, vector<ErrorInfo> &errorsInSortedReads );

    //make sure seqNum field for ErrorInfo objects refers to one of the original reads (as oppposed to their reverse complements)
    //and  make sure that the corrector strings are running in the correct direction (L2R for reverse strand, R2L otherwise)
    static void ConvertRCCorrectionsToOriginal( vector<ErrorInfo> &errors, int numberOfReads, int readLength );

    //write a vector of ErrorInfo objects to a csv file with the fields in the order that is understood by 'ReadCorrectionsFromCsv'
    static void CorrectionsToCsv( const string &fileName, vector<ErrorInfo> &corrections );
    //read a vector of ErrorInfo objects from a file written by CorrectionsToCsv
    static vector<ErrorInfo> ReadCorrectionsFromCsv( const string &fileName );

};

typedef map<LetterNumber, ErrorInfo> ErrorStore;

string strreverse( const string &inStr );

#endif
