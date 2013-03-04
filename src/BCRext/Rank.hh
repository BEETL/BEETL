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

// Encapsulate the functionality for rank and inverse-rank queries;
// operates on a single BWT file (agnostic to whether it's a partial
// or complete BWT)

#ifndef RANK_HH
#define RANK_HH

// data type definitions
#include "Tools.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>


/* class OccComparison  */
/* { */
/* public: */
/*   static int idx; */

/*   static struct _CompareFloatField */
/*   { */
/*     bool operator() (const vector<dataTypeNSeq>& left, int right) */
/*     { */
/*       return left[1]<right; */
/*     } */

/*   } CompareOcc; */
/* }; */


class StaticOccComparison
{
public:
    int idx;

    StaticOccComparison( const int &i )
    {
        idx=i;
    }

    bool operator() ( const vector<dataTypeNSeq> &left, dataTypeNSeq right )
    {
        return left[idx]<right;
    }
};



class Rank
{
public:
    Rank( void );
    ~Rank();

    // set methods
    bool setAlpha( dataTypedimAlpha *alpha,const int &s );
    bool setAlphaInverse( dataTypedimAlpha *alphaInverse,const int &s );
    void setCharCount( const int &cnt )
    {
        charCount_=cnt;
    }

    // read BWT and initialize all the counts
    bool initialize( const string &fn, const int &blockSize );

    // Given character x and a number n, get the blockID that contains
    // the n-th occurrence of x
    dataTypeNChar getInverseRank( const char &c,const dataTypeNChar &n );

private:
    FILE *ifsBWT; // input file streamn of BWT

    vector<dataTypedimAlpha> alpha_;
    vector<dataTypedimAlpha> alphaInverse_;
    bool alphaInitialized_;
    bool alphaInverseInitialized_;
    int charCount_;
    int blockSize_;

    // holding the actual counts
    vector< vector<dataTypeNSeq> > occ_; // addressing: occ_[blockID][char]
    vector< map<dataTypedimAlpha,dataTypedimAlpha> > mapCountToBlock_; // adressing mapCountToBlock_[char][countToFind]==blockID
    vector< StaticOccComparison > compareObjects_;

};


#endif
