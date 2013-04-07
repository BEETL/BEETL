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

#ifndef BEETL_COMPARE_PARAMETERS_HH
#define BEETL_COMPARE_PARAMETERS_HH

#include "ToolParameters.hh"

#include <string>


// options: input format
/*
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


// options: process qualities off/on

enum ProcessQualities
{
    PROCESS_QUALITIES_OFF,
    PROCESS_QUALITIES_ON,
    PROCESS_QUALITIES_COUNT,
};

static const string processQualitiesLabels[] =
{
    "off",
    "on",
    "" // end marker
};

*/

// options: mode split/reference/metagenomics

enum Mode
{
    MODE_SPLIT,
    MODE_REFERENCE,
    MODE_METAGENOMICS,
    MODE_COUNT,
};

static const string modeLabels[] =
{
    "split",
    "reference",
    "metagenomics",
    "" // end marker
};


// options: report minlength off/on

enum ReportMinlength
{
    REPORT_MINLENGTH_OFF,
    REPORT_MINLENGTH_ON,
    REPORT_MINLENGTH_COUNT,
};

static const string reportMinlengthLabels[] =
{
    "off",
    "on",
    "" // end marker
};



// Option container

enum CompareOptions
{
    /*
        OPTION_INPUT_FORMAT,
        OPTION_OUTPUT_FORMAT,
        OPTION_ALGORITHM,
        OPTION_INTERMEDIATE_FORMAT,
        OPTION_INTERMEDIATE_STORAGE_MEDIUM,
        OPTION_PARALLEL_PREFETCH,
        OPTION_PARALLEL_PROCESSING,
    COMPARE_OPTION_PROCESS_QUALITIES,
    */
    COMPARE_OPTION_MODE,
    COMPARE_OPTION_REPORT_MINLENGTH,
    COMPARE_OPTION_COUNT // end marker
};

static const string compareOptionNames[] =
{
    /*
        "input format",
        "output format",
        "algorithm",
        "intermediate format",
        "intermediate storage medium",
        "parallel prefetch",
        "parallel processing",
    "process qualities",
    */
    "mode",
    "report minlength",
    "" // end marker
};

static const string *compareOptionPossibleValues[] =
{
    /*
        inputFormatLabels,
        outputFormatLabels,
        algorithmOptionLabels,
        intermediateFormatLabels,
        intermediateStorageMediumLabels,
        parallelPrefetchLabels,
        parallelProcessingLabels,
    processQualitiesLabels,
    */
    modeLabels,
    reportMinlengthLabels,
    NULL // end marker
};


class CompareParameters : public ToolParameters
{
    virtual const string **getOptionPossibleValues() const
    {
        return compareOptionPossibleValues;
    }
public:
    virtual const string getOptionName( const unsigned i ) const
    {
        return compareOptionNames[i];
    }
    //    virtual const string getOptionPossibleValue( const unsigned i, const unsigned j ) const { return compareOptionPossibleValues[i][j]; }
};

#endif //ifndef BEETL_COMPARE_PARAMETERS_HH
