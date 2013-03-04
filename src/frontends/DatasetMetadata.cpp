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

#include "DatasetMetadata.hh"

#include "SeqReader.hh"
#include "TransposeFasta.hh"

using namespace std;


DatasetMetadata datasetMetadata;


void DatasetMetadata::init( const string &input, const string &inputFormat )
{
    if ( inputFormat == "cyc" )
    {
        // Use TransposeFasta's ability to extract cyc files metadata, even though we don't need to transpose it
        TransposeFasta transp;
        transp.inputCycFile( input );
        nReads = transp.nSeq;
        nCycles = transp.lengthRead;
    }
    else
    {
        FILE *f = fopen( input.c_str(),"rb" );
        SeqReaderFile *pReader( SeqReaderFile::getReader( f ) );
        nCycles = pReader->length();

        // At this point SeqReader will have read the first dataset entry
        // We estimate the total number of entries by dividing the total file size
        long entrySize = ftell( f );
        fseek( f, 0, SEEK_END );
        long fileSize = ftell( f );
        nReads = fileSize / entrySize;
        delete pReader;
    }

    nBases = nCycles * nReads;
    rleCompressibility = 2; // final bits per base - todo
}
