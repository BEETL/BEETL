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

#include "BWTCollection.hh"

#include <cassert>
#include <iostream>

using namespace std;
using SXSI::BWTCollection;


#define BCR_ID "$Id: BCR_BWTCollection.cpp,v 1.6 2011/11/28 16:38:32 acox Exp $"


BCR::BCR( const int mode, const string &in, const string &out,
          const CompressionFormatType outputCompression ) :
    mode_( mode ),
    outputCompression_( outputCompression )
{
    inFile_ = in;
    outFile_ = out;
}

void BCR::run( void )
{

    BWTCollection *BCRexternalBWT = BWTCollection::InitBWTCollection( inFile_, outFile_, mode_, outputCompression_ );

    //cout << "finished iteration, usage: " << timer << endl;

    delete BCRexternalBWT;
}
