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

#ifndef DATASET_METADATA_HH
#define DATASET_METADATA_HH

#include <string>

using std::string;


class DatasetMetadata
{
public:
    unsigned int nCycles;
    unsigned long nReads;
    unsigned long nBases;
    float rleCompressibility;

    void init( const string &input, const string &inputFormat );
};
extern DatasetMetadata datasetMetadata;


#endif //ifndef DATASET_METADATA_HH
