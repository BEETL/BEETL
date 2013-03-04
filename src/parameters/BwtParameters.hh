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

#ifndef BEETL_BWT_PARAMETERS_HH
#define BEETL_BWT_PARAMETERS_HH

#include <string>


static const int MULTIPLE_OPTIONS = 99;

// options: input format

enum InputFormat
{
    INPUT_FORMAT_FASTA,
    INPUT_FORMAT_FASTQ,
    INPUT_FORMAT_CYC,
    INPUT_FORMAT_SEQ,
    INPUT_FORMAT_BCL,
    INPUT_FORMAT_COUNT
};

static const string inputFormatLabels[] =
{
    "fasta",
    "fastq",
    "cyc",
    "seq",
    "bcl",
    "" // end marker
};


// options: output format

enum OutputFormat
{
    OUTPUT_FORMAT_ASCII,
    OUTPUT_FORMAT_RLE,
    OUTPUT_FORMAT_HUFFMAN,
    OUTPUT_FORMAT_COUNT
};

static const string outputFormatLabels[] =
{
    "ASCII",
    "RLE",
    "Huffman",
    "" // end marker
};


// options: algorithm

enum AlgorithmOption
{
    ALGORITHM_BCR,
    ALGORITHM_EXT,
    ALGORITHM_COUNT
};

static const string algorithmOptionLabels[] =
{
    "bcr",
    "ext",
    "" // end marker
};


// options: intermediate format

enum IntermediateFormat
{
    INTERMEDIATE_FORMAT_ASCII,
    INTERMEDIATE_FORMAT_RLE,
    INTERMEDIATE_FORMAT_MULTIRLE,
    INTERMEDIATE_FORMAT_HUFFMAN,
    INTERMEDIATE_FORMAT_COUNT
};

static const string intermediateFormatLabels[] =
{
    "ASCII",
    "RLE",
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


// options: parallel prefetch on/off

enum ParallelPrefetch
{
    PARALLEL_PREFETCH_OFF,
    PARALLEL_PREFETCH_ON,
    PARALLEL_PREFETCH_COUNT,
};

static const string parallelPrefetchLabels[] =
{
    "prefetch off",
    "prefetch on",
    "" // end marker
};


// options: parallel processing

enum ParallelProcessing // todo: this should only be on/off; the number of cores should be part of the hardware resources
{
    PARALLEL_PROCESSING_OFF,
    PARALLEL_PROCESSING_2CORES,
    PARALLEL_PROCESSING_4CORES,
    PARALLEL_PROCESSING_COUNT,
};

static const string parallelProcessingLabels[] =
{
    "1 core",
    "2 cores",
    "4 cores",
    "" // end marker
};


// options: permute qualities off/on

enum PermuteQualities
{
    PERMUTE_QUALITIES_OFF,
    PERMUTE_QUALITIES_ON,
    PERMUTE_QUALITIES_COUNT,
};

static const string permuteQualitiesLabels[] =
{
    "off",
    "on",
    "" // end marker
};


// options: generate LCP on/off

enum GenerateLcp
{
    GENERATE_LCP_ON,
    GENERATE_LCP_OFF,
    GENERATE_LCP_COUNT,
};

static const string generateLcpLabels[] =
{
    "lcp on",
    "lcp off",
    "" // end marker
};


// options: single cycle on/off

enum SingleCycle
{
    SINGLE_CYCLE_ON,
    SINGLE_CYCLE_OFF,
    SINGLE_CYCLE_COUNT,
};

static const string singleCycleLabels[] =
{
    "single cycle on",
    "single cycle off",
    "" // end marker
};


// options: concatenate output off/on

enum ConcatenateOutput
{
    CONCATENATE_OUTPUT_OFF,
    CONCATENATE_OUTPUT_ON,
    CONCATENATE_OUTPUT_COUNT,
};

static const string concatenateOutputLabels[] =
{
    "off",
    "on",
    "" // end marker
};


// options: sap ordering off/on

enum SapOrdering
{
    SAP_ORDERING_OFF,
    SAP_ORDERING_ON,
    SAP_ORDERING_COUNT,
};

static const string sapOrderingLabels[] =
{
    "off",
    "on",
    "" // end marker
};


// options: generate endPosFile off/on

enum GenerateEndPosFile
{
    GENERATE_ENDPOSFILE_OFF,
    GENERATE_ENDPOSFILE_ON,
    GENERATE_ENDPOSFILE_COUNT,
};

static const string generateEndPosFileLabels[] =
{
    "off",
    "on",
    "" // end marker
};


// Option container

enum Options
{
    OPTION_INPUT_FORMAT,
    OPTION_OUTPUT_FORMAT,
    OPTION_ALGORITHM,
    OPTION_INTERMEDIATE_FORMAT,
    OPTION_INTERMEDIATE_STORAGE_MEDIUM,
    OPTION_PARALLEL_PREFETCH,
    OPTION_PARALLEL_PROCESSING,
    OPTION_PERMUTE_QUALITIES,
    OPTION_GENERATE_LCP,
    OPTION_SINGLE_CYCLE,
    OPTION_CONCATENATE_OUTPUT,
    OPTION_SAP_ORDERING,
    OPTION_GENERATE_ENDPOSFILE,
    OPTION_COUNT // end marker
};

static const string optionNames[] =
{
    "input format",
    "output format",
    "algorithm",
    "intermediate format",
    "intermediate storage medium",
    "parallel prefetch",
    "parallel processing",
    "permute qualities",
    "generate LCP",
    "single cycle",
    "concatenate output",
    "SAP ordering",
    "generate endPosFile",
    "" // end marker
};

static const string *optionPossibleValues[] =
{
    inputFormatLabels,
    outputFormatLabels,
    algorithmOptionLabels,
    intermediateFormatLabels,
    intermediateStorageMediumLabels,
    parallelPrefetchLabels,
    parallelProcessingLabels,
    permuteQualitiesLabels,
    generateLcpLabels,
    singleCycleLabels,
    concatenateOutputLabels,
    sapOrderingLabels,
    generateEndPosFileLabels,
    NULL // end marker
};


class BwtParameters : public vector< int >
{
public:
    int getValue( const enum Options key ) const
    {
        int valNum = this->operator[]( key );
        return valNum;
    }

    string getValueAsString( const enum Options key ) const
    {
        int valNum = getValue( key );
        return optionPossibleValues[key][valNum];
    }

    void set( const enum Options key, const int val ) const
    {
        assert( false && "todo" );
    }
};

#endif //ifndef BEETL_BWT_PARAMETERS_HH
