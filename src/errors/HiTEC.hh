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

#ifndef INCLUDED_HITEC_HH
#define INCLUDED_HITEC_HH

#include "Algorithm.hh"
#include "OneBwtBackTracker.hh"
#include "Config.hh"
#include "LetterCount.hh"

#include "countWords/RangeStore.hh"
#include "Types.hh"
#include "parameters/BwtParameters.hh"
#include "WitnessReader.hh"
#include "HiTECErrorLocation.hh"
#include "HiTECIntervalHandler.hh"
#include "HiTECStats.hh"
#include <string>
#include <math.h>

using std::string;

#define HITEC_MAX_ITERATIONS 9

class HiTEC : public Algorithm
{

public:
    HiTEC(
        const string &readsFile,
        double errorRate,
        double genomeLength
    );

    virtual ~HiTEC() {}
    void showExecutionPlan();
    void run( void );
    void setWitnessLength( int witnessLength );
    void setCorrectThreshold( int minOccurrences );
private:
    bool compressIntermediateBwts_;

    int numberOfReads_;
    int readLength_;
    double genomeLength_;
    double errorRate_;

    //the minimum number of times a letter in an interval needs to appear to be considered the 'correct' letter
    //the HiTEC paper calls this 'T'
    int minSupport_;

    //the 'w little m' statistic as described in the HiTEC paper - calculated from the error rate, genome length, number of reads, read length.
    //represents a short witness, designed at catching as many errors as possible, at the risk of 'falsely accusing' correct bases of being errors
    int wm_;

    //the 'w big m' statistic - chosen to achieve a certain level of specificity (and is greater than 'w little m')
    int wM_;

    //the sequence of witness lengths to use throughout the iterations of HiTEC...
    int witnessLengthSequence_[9];

    //the maximum false positive rate to allow in calculating wM (however, it doesnt translate to a true max false positive rate, as HiTEC uses a combination
    //of witness lengths anyway). The original HiTEC sets this to 0.0001, so we use this as our default.
    double proportionCorrectedThreshold_;

    bool doesPropagateBkptToSeqNumInSet_;

    //no parallel version, as yet
    string subset_;

    string readsFile_;
    string tempFilesPrefix_;
    string workingReadsFile_;

    RangeStoreExternal *getPutativeErrors( int witnessLength );
    vector<HiTECErrorLocation> findErrorsInSortedReads( int witnessLength, RangeStoreExternal *errorIntervals );

    //HiTEC runs in iterations, and at each one, corrections are calculated and applied. So a new BWT must be calculated each iteration.
    void makeBwtAndLcp( void );

    string lcpFileName( int letterNumber );
    string bwtFileName( int letterNumber );

    LetterCountEachPile cumulativePileCounts();

    void setRequiredAttributes(
        const string &readsFile,
        double errorRate,
        double genomeLength
    );

    //set the readLength_ and numberOfReads_ variables...
    void setImpliedAttributes();

    //set the witnessLengthSequence_ and minimumSupport_ variables, based on HiTECStats methods, as described in the HiTEC paper.
    void setRecommendedAttributes();

    //a (not very elegant) cleaning up method
    void removeUsedFiles();

};

#endif
