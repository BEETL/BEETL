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

#ifndef INCLUDED_BCREXT_HH
#define INCLUDED_BCREXT_HH

#include <string>
#include "Algorithm.hh"
using namespace std;

// by Tobias, small class interface to call from beetl executable 

class BCRext : public Algorithm {
    bool useHuffmanEncoder_;
    bool useAsciiEncoder_;
    bool useRunlengthEncoder_;
    char* inFile_;
    char* prefix_;
public:
    BCRext(bool, bool, bool, string, string);
    ~BCRext() { delete[] inFile_; delete prefix_; }
    void run(void);
};
// end


#define BCREXT_HH_ID "@(#) $Id: BCRext.hh,v 1.2 2011/11/28 11:32:09 tjakobi Exp $"




#endif


