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

#ifndef INCLUDED_BCREXT_HH
#define INCLUDED_BCREXT_HH

#include "Algorithm.hh"

#include <string>

using std::string;


// by Tobias, small class interface to call from beetl executable

class BCRext : public Algorithm
{
    const bool useHuffmanEncoder_;
    const bool useRunlengthEncoder_;
    const bool useAsciiEncoder_;
    const bool useImplicitSort_;
    const bool useSeqFile_;
    const string inFile_;
    const string prefix_;
public:
    BCRext( bool, bool, bool, bool, bool, const string &, const string & );
    ~BCRext()
    {
    }
    void run( void );
};
// end


#define BCREXT_HH_ID "@(#) $Id: BCRext.hh,v 1.2 2011/11/28 11:32:09 tjakobi Exp $"




#endif


