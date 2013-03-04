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

#ifndef INCLUDED_COUNTWORDS_HH
#define INCLUDED_COUNTWORDS_HH

#include "Algorithm.hh"
#include "BackTracker.hh"
#include "Config.hh"
#include "LetterCount.hh"
#include "RangeStore.hh"
#include "Types.hh"
#include "IntervalHandlerBase.hh"
#include "IntervalHandlerReference.hh"
#include "IntervalHandlerMetagenome.hh"
#include "IntervalHandlerSplice.hh"

#include <string>

using std::string;


//#define PROPAGATE_PREFIX 1

class CountWords : public Algorithm
{

public:
    CountWords( bool, bool, bool, char, int, int, const vector<string> &, const vector<string> &,
                const vector<string> &, const string &, bool testDB, uint minWord, string prefix );

    virtual ~CountWords() {}
    void run( void );
    void loadFileNumToTaxIds( string taxIdNames );

private:
    bool bothSetsCompressed_;
    bool inputACompressed_;
    bool inputBCompressed_;
    char whichHandler_;

    int paramN_;
    int paramK_;

    vector <string> setA_;
    vector <string> setB_;
    //Christina: this is needed for the metagenomics variation of the count word algorithm
    // set C of the merging algorithm, it contains for each BWT positions the fileNumber from which the suffix came
    vector <string> setC_;
    //FileNumberToTaxTree. this should contain the taxonomic information of the database.
    //for each file number (corresponding to the Merge C set) there should be the taxonomic Ids ordered this way:
    //Filenumber superkingdom phylum class order family genus species strain
    //if one of those numbers is not available it should be replaced with a 0
    string ncbiInfo_;
    // flag to indicate if the database should be tested
    bool testDatabase_;
    //Christina: added a minimal word length for the output else the parsing of the results takes to long for my liking
    uint minWordLen_;
    // Prefix for temp output files
    string tmpPrefix_;
    //for each fileNumber there should be the same amount of taxIds, this can stop at any level it will be filled up with zeros
    vector< vector< int> > fileNumToTaxIds_;
    //boolean to indicate if the database should be tested or not
    bool testDB_;
    //minimal word Lenght. Not exactly needed but speeds the algorithm up,
    //because the fileNumbers and taxa do not have to be checked before the minimal wordLength is reached
    //uint minWordLength_;
};
#endif
