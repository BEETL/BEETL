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

#ifndef DATASET_METADATA_HH
#define DATASET_METADATA_HH

#include "Types.hh"

#include <string>

using std::string;


class DatasetMetadata
{
public:
    SequenceLength nCycles;
    SequenceNumber nReads;
    LetterNumber nBases;
    float rleCompressibility;

    void init( const string &input, const string &inputFormat );
};
extern DatasetMetadata datasetMetadata;


#endif //ifndef DATASET_METADATA_HH
