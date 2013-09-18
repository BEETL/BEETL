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

#ifndef ERROR_HITECPARAMETERS_HH
#define ERROR_HITECPARAMETERS_HH

#include "libzoo/cli/ToolParameters.hh"

#include <string>

using std::string;


// options: input format

enum InputFormat
{
    INPUT_FORMAT_BWT_RLE,
    INPUT_FORMAT_FASTA,
    INPUT_FORMAT_FASTQ,
    INPUT_FORMAT_CYC,
    INPUT_FORMAT_SEQ,
    INPUT_FORMAT_BCL,
    INPUT_FORMAT_COUNT
};

static const string inputFormatLabels[] =
{
    "bwt-rle",
    "fasta",
    "fastq",
    "cyc",
    "seq",
    "bcl",
    "" // end marker
};

// options: algorithm

enum AlgorithmOption
{
    ALGORITHM_HITEC,
    ALGORITHM_BWT
};

static const string algorithmLabels[] =
{
    "HiTEC",
    "BWT",
    "" // end marker
};


// options: intermediate format

enum IntermediateFormat
{
    INTERMEDIATE_FORMAT_RLE,
    INTERMEDIATE_FORMAT_ASCII
};

static const string intermediateFormatLabels[] =
{
    "RLE",
    "ASCII",
    "multiRLE",
    "Huffman",
    "" // end marker
};


// options: intermediate storage medium

enum IntermediateStorageMedium
{
    INTERMEDIATE_STORAGE_MEDIUM_DISK,
    INTERMEDIATE_STORAGE_MEDIUM_RAM,
    INTERMEDIATE_STORAGE_MEDIUM_COUNT
};

static const string intermediateStorageMediumLabels[] =
{
    "disk",
    "RAM",
    "" // end marker
};


// Option container

enum HiTECParameterIds
{
    PARAMETER_UNDEFINED = -1,
    PARAMETER_INPUT_FILENAME = 0,
    PARAMETER_INTERMEDIATE_FORMAT,
    PARAMETER_GENOME_LENGTH,
    PARAMETER_ERROR_RATE,
    PARAMETER_READ_LENGTH,
    PARAMETER_MIN_WITNESS_LENGTH,
    PARAMETER_DONT_RUN
};

class HiTECParameters : public ToolParameters
{

public:
    HiTECParameters()
    {
        addEntry( PARAMETER_INPUT_FILENAME, "input filename", "--input", "-i", "Input reads file you want to correct - warning! It gets overwritten - copy your reads first", "", TYPE_STRING | REQUIRED );
        addEntry( PARAMETER_GENOME_LENGTH, "genome length", "--genome-length", "-L", "Approximate length of genome from which your reads came", "", TYPE_INT | REQUIRED );
        addEntry( PARAMETER_ERROR_RATE, "error rate", "--error-rate", "-e", "Rate of errors on sequencing platform generating reads", "", TYPE_INT | REQUIRED );
        addEntry( PARAMETER_MIN_WITNESS_LENGTH, "min witness length", "--min-witness-length", "-w", "Minimum witness length", "", TYPE_INT | AUTOMATED );
        addEntry( PARAMETER_DONT_RUN, "don't run", "--dont-run", "-X", "Don't run the algorithm - just show execution plan", "", TYPE_SWITCH );
        addEntry( -1, "subset", "--subset", "", "Restrict computation to this suffix - Used for distributed computing", "", TYPE_STRING );

        addDefaultVerbosityAndHelpEntries();
    }

    static void ReadLengthAndCount( const string &readsFile, int &numberOfReads, int &readLength );
};


#endif