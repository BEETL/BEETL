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

#ifndef KMER_SEARCH_INTERVAL_HANDLER_HH
#define KMER_SEARCH_INTERVAL_HANDLER_HH

#include "KmerSearchRange.hh"
#include "IntervalHandlerBase.hh"
#include "RangeStore.hh"


struct KmerSearchItem
{
    KmerSearchItem(
        string inKmer = "",
        LetterNumber inPosition = 0,
        LetterNumber inCount = 0,
        SequenceNumber index = 0
    ):
        kmer( inKmer ),
        position( inPosition ),
        count( inCount ),
        originalIndex( index )
    {}

    bool operator<( const KmerSearchItem &rhs ) const
    {
        return ( kmer.compare( rhs.kmer ) < 0 );
    }
    //	void print();

    string kmer;
    LetterNumber position;
    LetterNumber count;
    SequenceNumber originalIndex; // index in list of input kmers
};


struct KmerSearchIntervalHandler : public IntervalHandlerBase
{
    KmerSearchIntervalHandler();

    virtual ~KmerSearchIntervalHandler() {}

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
    )
    {
        assert( false );
    }

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
    )
    {
        assert( false );
    }

    virtual Range &getSubIntervalRange (
        const string &word,
        const LetterNumber pos,
        const LetterNumber num,
        const bool isBkptExtension,
        Range &parentRange,
        const int subIntervalNum
    )
    {
        static KmerSearchRange r;
        r = KmerSearchRange( word, pos, num, isBkptExtension, parentRange, subIntervalNum );
        return r;
    }

};

#endif // KMER_SEARCH_INTERVAL_HANDLER_HH
