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

#ifndef INCLUDED_BWTCORRECTOR_HH
#define INCLUDED_BWTCORRECTOR_HH

#include "Algorithm.hh"
#include "OneBwtBackTracker.hh"
#include "Config.hh"
#include "LetterCount.hh"
#include "HiTECStats.hh"
#include "countWords/RangeStore.hh"
#include "Types.hh"
#include "parameters/BwtParameters.hh"
#include "WitnessReader.hh"
#include "ErrorInfo.hh"
#include "BwtCorrectorIntervalHandler.hh"
#include <string>
#include <math.h>

using std::string;

class BwtCorrector : public Algorithm
{

public:
    BwtCorrector(
        const string &inputFile,
        const string &outputFile,
        int numberOfReads,
        int readLength,
        double errorRate,
        double genomeLength,
        int minWitnessLength,
        const string &subset
    )
        : genomeLength_( genomeLength )
        , indexPrefix_( inputFile )
        , errorRate_( errorRate )
        , subset_( subset )
        , numberOfReads_( numberOfReads )
        , readLength_( readLength )
        , minWitnessLength_( minWitnessLength )
        , outputFile_( outputFile )
    {}

    virtual ~BwtCorrector() {}
    void showExecutionPlan();
    ErrorStore findErrors();
    void run( void );

private:
    double genomeLength_;
    double errorRate_;

    int numberOfReads_;
    int readLength_;

    int minWitnessLength_;

    string subset_;

    string readsFile_;
    string indexPrefix_;
    string outputFile_;
};

#endif
