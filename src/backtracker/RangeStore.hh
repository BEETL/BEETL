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

#ifndef INCLUDED_RANGESTORE_HH
#define INCLUDED_RANGESTORE_HH

#include "Config.hh"
#include "Alphabet.hh"
#include "LetterCount.hh"
#include "Range.hh"
#include "Tools.hh"
#include "libzoo/util/TemporaryFilesManager.hh"

#include <string>


//
// abstract RangeStore - manage sets of intervals in BWT
//
struct RangeStore
{
    virtual ~RangeStore() {}
    //    virtual void swap( void ) = 0;
    virtual void setCycleNum( const int cycleNum ) = 0;
    virtual void setPortion( int pileNum, int portionNum ) = 0;
    virtual bool getRange( Range &thisRange ) = 0;
    virtual void addRange( const Range &, const int pileNum, const int portionNum, const string &subset, const int cycle ) = 0;
    virtual void addOutOfOrderRange( const Range &, const int pileNum, const AlphabetSymbol portionNum, const string &subset, const int cycle ) = 0;
    virtual bool isRangeKnown( const Range &, const int pileNum, const int portionNum, const string &subset, const int cycle ) = 0;

protected:
    bool isSubsetValid( const string &subset, const int cycle, const int pileNum, const int portionNum );
}; // ~struct RangeStore


//
// RangeStoreExternal: Store to disk
//
struct RangeStoreExternal : public RangeStore
{
    RangeStoreExternal( const bool propagateSequence = false, const string fileStem = "Intervals" );

    virtual ~RangeStoreExternal();

    //    virtual void swap( void );
    virtual void setCycleNum( const int cycleNum );

    virtual void setPortion( int pileNum, int portionNum );
    virtual void deleteInputPortion( int pileNum, int portionNum );

    virtual bool getRange( Range &thisRange );

    virtual void addRange( const Range &, const int pileNum, const int portionNum, const string &subset, const int cycle );
    virtual void addOutOfOrderRange( const Range &, const int pileNum, const AlphabetSymbol portionNum, const string &subset, const int cycle );
    virtual bool isRangeKnown( const Range &, const int pileNum, const int portionNum, const string &subset, const int cycle );

    virtual void clear( bool doDeleteFiles = true );

    void getFileName( const string &stem, const int pile, const int portion,
                      string &fileName );

    string fileStem_;
    string fileStemIn_;
    string fileStemOut_;

    RangeState stateIn_;
    vector< vector< RangeState > > stateOut_; //[alphabetSize][alphabetSize];

private:
    bool getRange( RangeState &stateFile, Range &thisRange );

    vector< vector< RangeState > > stateInForComparison_; //[alphabetSize][alphabetSize];
    Range lastRangeReadForComparison_[alphabetSize][alphabetSize];
    vector< std::pair< Range, AlphabetSymbol > > outOfOrderRangesForPile0_;
}; // ~struct RangeStoreExternal

#endif
