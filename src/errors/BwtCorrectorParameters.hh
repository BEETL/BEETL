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

#ifndef INCLUDED_BWTCORRECTORPARAMETERS_HH
#define INCLUDED_BWTCORRECTORPARAMETERS_HH

#include "libzoo/cli/ToolParameters.hh"

#include <string>

using std::string;


namespace BeetlCorrectParameters
{

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

enum BwtCorrectorParameterIds
{
    PARAMETER_UNDEFINED = -1,
    PARAMETER_INPUT_FILENAME = 0,
    PARAMETER_INTERMEDIATE_FORMAT,
    PARAMETER_GENOME_LENGTH,
    PARAMETER_ERROR_RATE,
    PARAMETER_READ_LENGTH,
    PARAMETER_MIN_WITNESS_LENGTH,
    PARAMETER_DONT_RUN,
    PARAMETER_WRITE_CORRECTED_READS_ONLY,
    PARAMETER_CORRECTIONS_FILE,
    PARAMETER_CORRECTED_READS_FILE,
    PARAMETER_MIN_SUPPORT
};

} // namespace BeetlCorrectParameters



class BwtCorrectorParameters : public ToolParameters
{

public:
    BwtCorrectorParameters()
    {
        using namespace BeetlCorrectParameters;
        addEntry( PARAMETER_INPUT_FILENAME, "input filename", "--input", "-i", "Prefix of BWT of reads you want to correct", "", TYPE_STRING | REQUIRED );
        addEntry( PARAMETER_GENOME_LENGTH, "genome length", "--genome-length", "-L", "Approximate length of genome from which your reads came", "", TYPE_INT | REQUIRED );
        addEntry( PARAMETER_ERROR_RATE, "error rate", "--error-rate", "-e", "Rate of errors on sequencing platform generating reads", "", TYPE_INT | REQUIRED );
        addEntry( PARAMETER_READ_LENGTH, "read length", "--read-length", "-k", "Length of reads", "", TYPE_INT | REQUIRED );
        addEntry( PARAMETER_MIN_WITNESS_LENGTH, "min witness length", "--min-witness-length", "-w", "Minimum witness length", "", TYPE_INT | AUTOMATED );
        addEntry( PARAMETER_DONT_RUN, "don't run", "--dont-run", "-X", "Don't run the algorithm - just show execution plan", "", TYPE_SWITCH );
        addEntry( -1, "subset", "--subset", "", "Restrict computation to this suffix - Used for distributed computing", "", TYPE_STRING );
        addEntry( PARAMETER_CORRECTIONS_FILE, "corrections output filename", "--corrections-file", "-o", "File to which corrections are written", "", TYPE_STRING | REQUIRED );
        addEntry( PARAMETER_MIN_SUPPORT, "min support", "--minimum-support", "", "Fixed minimum occurrences for a base in an interval to be 'correct'", "", TYPE_INT );
        addDefaultVerbosityAndHelpEntries();
    }

    static void ReadLengthAndCount( const string &readsFile, int &numberOfReads, int &readLength );
};


#endif //ifndef INCLUDED_BWTCORRECTORPARAMETERS_HH
