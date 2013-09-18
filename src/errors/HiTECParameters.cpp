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

#include <iostream>
#include <fstream>
#include "HiTECParameters.hh"
#include "../shared/SeqReader.hh"

void HiTECParameters::ReadLengthAndCount( const string &readsFile, int &numberOfReads, int &readLength )
{
    FILE *reads = fopen( readsFile.c_str(), "r" );
    if ( !reads )
    {
        cerr << "Error opening file " << readsFile << endl;
        exit( -1 );
    }

    SeqReaderFile *seqReader = SeqReaderFile::getReader( reads );

    readLength = seqReader->length();

    numberOfReads = 0;
    while ( !seqReader->allRead() )
    {
        seqReader->readNext();
        numberOfReads++;
    }

    fclose( reads );
    delete seqReader;

}