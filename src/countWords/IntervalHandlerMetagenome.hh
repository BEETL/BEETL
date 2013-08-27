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

#ifndef INCLUDED_INTERVALHANDLER_METAGENOME_HH
#define INCLUDED_INTERVALHANDLER_METAGENOME_HH

#include "Alphabet.hh"
#include "IntervalHandlerBase.hh"
#include "LetterCount.hh"
#include "RangeStore.hh"
#include "Tools.hh"
#include "Types.hh"

#include <cstdlib>
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
                               vector<FILE *> mergeCSet,
                               vector< vector< int> > fileNumToTaxIds,
                               bool testDB,
                               uint minWordLength,
                               uint maxWordLength )
        : minOcc_( minOcc ), mergeCSet_( mergeCSet ),
          fileNumToTaxIds_( fileNumToTaxIds ), testDB_( testDB ),
          minWordLength_( minWordLength ),
          maxWordLength_( maxWordLength ) {}

    virtual ~IntervalHandlerMetagenome() {}

    virtual void foundInBoth
    ( const int pileNum,
      const LetterCount &countsThisRangeA, const LetterCount &countsThisRangeB,
      const Range &thisRangeA, const Range &thisRangeB,
      AlphabetFlag &propagateIntervalA, AlphabetFlag &propagateIntervalB,
      bool &isBreakpointDetected );

    virtual void foundInAOnly
    ( const int pileNum,
      const LetterCount &countsSoFarA,
      const LetterCount &countsThisRangeA,
      const Range &thisRangeA,
      AlphabetFlag &propagateIntervalA );

    virtual void foundInBOnly
    ( const int pileNum,
      const LetterCount &countsSoFarB,
      const LetterCount &countsThisRangeB,
      const Range &thisRangeB,
      AlphabetFlag &propagateIntervalB );

    vector<bool> intervalInSameTaxa( vector<uint> &sharedTaxIds, vector<MetagFileNumRefType> &fileNumbers );
    void getFileNumbersForRange( const int &pileNum, const LetterNumber &bwtPosition, const uint &num, vector<MetagFileNumRefType> &fileNumbers );

    const LetterNumber minOcc_;
    //setC of the mering algorithm from Tony,
    //for each bwt positions there should be a (unsigned short) fileNumber indicating from which file the suffix came from
    vector<FILE *> mergeCSet_;
    //for each fileNumber there should be the same amount of taxIds, this can stop at any level it will be filled up with zeros
    vector< vector< int> > fileNumToTaxIds_;

    bool testDB_;
    //minimal word length. Not exactly needed but speeds the algorithm up,
    //because the fileNumbers and taxa do not have to be checked before the minimal wordLength is reached
    uint minWordLength_;
    uint maxWordLength_;
};


#endif
