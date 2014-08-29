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

#ifndef INCLUDED_INTERVALHANDLER_METAGENOME_HH
#define INCLUDED_INTERVALHANDLER_METAGENOME_HH

#include "Alphabet.hh"
#include "IntervalHandlerBase.hh"
#include "LetterCount.hh"
#include "RangeStore.hh"
#include "Tools.hh"
#include "Types.hh"

#include <cstdlib>
#include <fstream>
#include <unistd.h>

using std::vector;


//
// IntervalHandler
//
// Idea here is that different algorithms can be implemented by defining
// new subclasses of IntervalHandler
struct IntervalHandlerMetagenome : public IntervalHandlerBase
{
    IntervalHandlerMetagenome( unsigned int minOcc,
                               vector<string> &filenamesCSet,
                               const vector<char *> &mmappedCFiles,
                               vector< vector< int> > &fileNumToTaxIds,
                               bool testDB,
                               uint minWordLength,
                               uint maxWordLength );

    virtual ~IntervalHandlerMetagenome();

    virtual void foundInBoth
    ( const int pileNum,
      const LetterCount &countsThisRangeA,
      const LetterCount &countsThisRangeB,
      const Range &thisRangeA,
      const Range &thisRangeB,
      AlphabetFlag &propagateIntervalA,
      AlphabetFlag &propagateIntervalB,
      bool &isBreakpointDetected,
      const int cycle
    );

    virtual void foundInAOnly
    ( const int pileNum,
      const LetterCount &countsSoFarA,
      const LetterCount &countsThisRangeA,
      const char *bwtSubstring,
      Range &thisRangeA,
      AlphabetFlag &propagateIntervalA,
      const int cycle
    );

    virtual void foundInBOnly
    ( const int pileNum,
      const LetterCount &countsSoFarB,
      const LetterCount &countsThisRangeB,
      const char *bwtSubstring,
      Range &thisRangeB,
      AlphabetFlag &propagateIntervalB,
      const int cycle
    );

    vector<bool> intervalInSameTaxa( vector<uint> &sharedTaxIds, vector<MetagFileNumRefType> &fileNumbers );
    void getFileNumbersForRange( const unsigned int &pileNum, const LetterNumber &bwtPosition, const uint &num, vector<MetagFileNumRefType> &fileNumbers );

    const LetterNumber minOcc_;
    //setC of the merging algorithm from Tony,
    //for each bwt positions there should be a (unsigned short) fileNumber indicating from which file the suffix came from
    vector<int> cSetFileDescs_;
    vector<off_t> posInFile_;
    const vector<char *> mmappedCFiles_;
    //for each fileNumber there should be the same amount of taxIds, this can stop at any level it will be filled up with zeros
    vector< vector< int> > &fileNumToTaxIds_;

    bool testDB_;
    //minimal word length. Not exactly needed but speeds the algorithm up,
    //because the fileNumbers and taxa do not have to be checked before the minimal wordLength is reached
    uint minWordLength_;
    uint maxWordLength_;

    //    void createOutputFile( const int subsetThreadNum, const int i, const int j, const int cycle );
    //    std::ofstream outFile_;
};


#endif
