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

#ifndef INCLUDED_COUNTWORDS_HH
#define INCLUDED_COUNTWORDS_HH

#include "Algorithm.hh"
#include "TwoBwtBackTracker.hh"
#include "Config.hh"
#include "LetterCount.hh"
#include "RangeStore.hh"
#include "Types.hh"
#include "IntervalHandlerBase.hh"
#include "IntervalHandlerReference.hh"
#include "IntervalHandlerMetagenome.hh"
#include "IntervalHandlerSplice.hh"
#include "IntervalHandlerTumourNormal.hh"
#include "parameters/CompareParameters.hh"

#include <string>

using std::string;


class CountWords : public Algorithm
{

public:
    CountWords( bool, bool, char, int, int, const vector<string> &, const vector<string> &
                , const vector<string> &, const string &, bool testDB, uint minWord, string subset
                , const CompareParameters *compareParams = NULL
              );

    virtual ~CountWords() {}
    virtual void run( void );
    void loadFileNumToTaxIds( const string &taxIdNames );
    bool isDistributedProcessResponsibleForPile( const int pile );

private:
    void initialiseMetagomeMode();
    void releaseMetagomeMode();
    void CountWords_parallelSubsetThread(
        const int subsetThreadNum
        , const int cycle
        , RangeStoreExternal &rangeStoreA
        , RangeStoreExternal &rangeStoreB
    );

    bool inputACompressed_;
    bool inputBCompressed_;
    enum BeetlCompareParameters::Mode mode_;

    const int numCycles_;
    const int minOcc_;

    vector <string> setA_;
    vector <string> setB_;
    //Christina: this is needed for the metagenomics variation of the count word algorithm
    // set C of the merging algorithm, it contains for each BWT positions the fileNumber from which the suffix came
    vector <string> setC_;
    vector<char *> mmappedCFiles_;
    //FileNumberToTaxTree. this should contain the taxonomic information of the database.
    //for each file number (corresponding to the Merge C set) there should be the taxonomic Ids ordered this way:
    //Filenumber superkingdom phylum class order family genus species strain
    //if one of those numbers is not available it should be replaced with a 0
    string ncbiInfo_;
    //Christina: added a minimal word length for the output else the parsing of the results takes to long for my liking
    uint minWordLen_;
    //for each fileNumber there should be the same amount of taxIds, this can stop at any level it will be filled up with zeros
    vector< vector< int> > fileNumToTaxIds_;
    //
    bool testDB_;
    //minimal word Lenght. Not exactly needed but speeds the algorithm up,
    //because the fileNumbers and taxa do not have to be checked before the minimal wordLength is reached
    //uint minWordLength_;
    const string subset_; // Compute only the results that have this string as a suffix' suffix
    const CompareParameters *compareParams_;

    // Switches and parameters
    const bool doPauseBetweenCycles_;
    const bool doesPropagateBkptToSeqNumInSetA_;
    const bool doesPropagateBkptToSeqNumInSetB_;
    bool noComparisonSkip_;
    const bool bwtInRam_;
    const bool propagateSequence_;
    const string outputDirectory_;

    // Used for computations
    vector <BwtReaderBase *> inBwtA_;
    vector <BwtReaderBase *> inBwtB_;
    LetterNumber numRanges_;
    LetterNumber numSingletonRanges_;
    LetterNumber numNotSkippedA_;
    LetterNumber numNotSkippedB_;
    LetterNumber numSkippedA_;
    LetterNumber numSkippedB_;
    LetterCountEachPile countsCumulativeA_;
    LetterCountEachPile countsCumulativeB_;
    double fsizeRatio_;
};
#endif
