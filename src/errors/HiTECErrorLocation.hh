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

#ifndef INCLUDED_HITECERRORLOCATION_HH
#define INCLUDED_HITECERRORLOCATION_HH

using namespace std;

string strreverse( const string &inStr );

struct HiTECErrorLocation
{
    HiTECErrorLocation() :
        readNum( 0 ),
        positionInRead( 0 ),
        corrector( ' ' )
    {}

    HiTECErrorLocation(
        int inReadNum,
        int inPositionInRead,
        char inCorrector
    ):
        readNum( inReadNum ),
        positionInRead( inPositionInRead ),
        corrector( inCorrector )
    {}

    void print()
    {
        cout << "Error located in read " << readNum;
        cout << " at position " << positionInRead;
        cout << " - should be " << corrector << endl;
    }

    int readNum;
    int positionInRead;
    char corrector;

    static void SetReadNumbersToOriginal( char *endPosFileName, vector<HiTECErrorLocation> &errorsInSortedReads );
    static void ConvertRCCorrectionsToOriginal( vector<HiTECErrorLocation> &errors, int numberOfReads, int readLength );
    static bool ErrorLocationSorter( HiTECErrorLocation a, HiTECErrorLocation b );
    static void CorrectionsToCsv( const string &fileName, vector<HiTECErrorLocation> &errorLocations );
    static vector<HiTECErrorLocation> ReadCorrectionsFromCsv( const string &fileName );

};


#endif