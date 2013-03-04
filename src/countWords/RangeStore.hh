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

#ifndef INCLUDED_RANGESTORE_HH
#define INCLUDED_RANGESTORE_HH

#include "Config.hh"
#include "LetterCount.hh"
#include "Tools.hh"

#include <string>


typedef unsigned short NumberFrag;
const int bitsPerFrag( ( 8*sizeof( NumberFrag ) )-1 );
const NumberFrag needAnotherFrag( 1<<bitsPerFrag );
const NumberFrag fragMask( ( unsigned short ) ~needAnotherFrag );
const LetterCountType letterCountMask( ( LetterCountType )fragMask );
const int fragBufSize( 8 );

struct Range
{
    Range( const string &word, LetterCountType pos, LetterCountType num ) :
#ifdef PROPAGATE_PREFIX
        word_( word ),
#endif
        pos_( pos ), num_( num ) {}
    Range( void )  {}
#ifdef PROPAGATE_PREFIX
    string word_;
#endif
    LetterCountType pos_;
    LetterCountType num_;
};

struct AllRanges: public vector<vector<vector<Range> > >
{
    AllRanges( void ) : vector<vector<vector<Range > > >
        ( alphabetSize, vector<vector<Range> >( alphabetSize ) ) {}
};

//
// abstract RangeStore - manage sets of intervals in BWT
//

struct RangeStore
{
    virtual ~RangeStore() {}
    virtual void swap( void )=0;
    virtual void setPortion( int pileNum, int portionNum )=0;
    virtual bool getRange( Range &thisRange )=0;
    virtual void addRange( int pileNum, int portionNum, const string &seq,
                           LetterCountType pos, LetterCountType num )=0;


}; // ~struct RangeStore

//
// RangeStoreRAM - hold BWT intervals in guess where
//

struct RangeStoreRAM : public RangeStore
{
    RangeStoreRAM();
    virtual ~RangeStoreRAM();

    virtual void swap( void );

    virtual void setPortion( int pileNum, int portionNum );

    virtual bool getRange( Range &thisRange );
    virtual void addRange( int pileNum, int portionNum, const string &seq,
                           LetterCountType pos, LetterCountType num );

    // clear range store ready for next iter
    virtual void clear( void );

    AllRanges r1, r2;
    AllRanges *pThis,*pNext,*pTemp;
    vector<Range>::iterator i_;
    vector<Range>::iterator end_;

}; // ~struct RangeStoreRAM

//
// RangeState: this is a proxy class that sits in RangeStoreExternal
// and handles the conversion between compressed representations on file
// and the structs the code expects
//


struct RangeState
{
    RangeState();
    void clear( void );

    void addSeq( const string &seq );
    void getSeq( string &word );

    void addNum( LetterCountType num );
    bool getNum( LetterCountType &num );

    NumberFrag fragBuf_[fragBufSize];
    char wordLast_[256];
    LetterCountType posLast_;
    FILE *pFile_;
}; // ~struct RangeState


//
// RangeStoreExternal
//

struct RangeStoreExternal : public RangeStore
{
    RangeStoreExternal( const string fileStemIn="countA",
                        const string fileStemOut="countB" );

    virtual ~RangeStoreExternal();

    virtual void swap( void );

    virtual void setPortion( int pileNum, int portionNum );

    virtual bool getRange( Range &thisRange );
    virtual void addRange( int pileNum, int portionNum, const string &seq,
                           LetterCountType pos, LetterCountType num );

    virtual void clear( void );

    void getFileName( const string &stem, const int pile, const int portion,
                      string &fileName );

    string fileStemIn_;
    string fileStemOut_;

    RangeState stateIn_;
    RangeState stateOut_[alphabetSize][alphabetSize];
    char buf_[256]; // get rid, not nice

}; // ~struct RangeStoreExternal

#endif
