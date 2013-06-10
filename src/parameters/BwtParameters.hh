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

#include "ToolParameters.hh"

#include <string>


namespace BeetlBwtParameters
{

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
    INTERMEDIATE_FORMAT_RLE,
    INTERMEDIATE_FORMAT_ASCII,
    INTERMEDIATE_FORMAT_MULTIRLE,
    INTERMEDIATE_FORMAT_HUFFMAN,
    INTERMEDIATE_FORMAT_COUNT
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


// options: process qualities

enum ProcessQualities
{
    PROCESS_QUALITIES_IGNORE,
    PROCESS_QUALITIES_PERMUTE,
    //    PROCESS_QUALITIES_PRESERVE,
    //    PROCESS_QUALITIES_SMOOTH,
    PROCESS_QUALITIES_COUNT,
};

static const string processQualitiesLabels[] =
{
    "ignore",
    "permute",
    //    "preserve",
    //    "smooth",
    "" // end marker
};


// options: generate LCP on/off

enum GenerateLcp
{
    GENERATE_LCP_OFF,
    GENERATE_LCP_ON,
    GENERATE_LCP_COUNT,
};

static const string generateLcpLabels[] =
{
    "off",
    "on",
    "" // end marker
};


// options: reverse off/on

enum Reverse
{
    REVERSE_OFF,
    REVERSE_ON,
    REVERSE_COUNT,
};

static const string reverseLabels[] =
{
    "off",
    "on",
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


// options: pause between cycles off/on

enum PauseBetweenCycles
{
    PAUSE_BETWEEN_CYCLES_OFF,
    PAUSE_BETWEEN_CYCLES_ON,
    PAUSE_BETWEEN_CYCLES_COUNT,
};

static const string pauseBetweenCyclesLabels[] =
{
    "off",
    "on",
    "" // end marker
};


// Option container

enum BwtOptions
{
    BWT_OPTION_INPUT_FORMAT,
    BWT_OPTION_OUTPUT_FORMAT,
    BWT_OPTION_ALGORITHM,
    BWT_OPTION_INTERMEDIATE_FORMAT,
    BWT_OPTION_INTERMEDIATE_STORAGE_MEDIUM,
    BWT_OPTION_PARALLEL_PREFETCH,
    BWT_OPTION_PARALLEL_PROCESSING,
    BWT_OPTION_PROCESS_QUALITIES,
    BWT_OPTION_GENERATE_LCP,
    BWT_OPTION_REVERSE,
    BWT_OPTION_SINGLE_CYCLE,
    BWT_OPTION_CONCATENATE_OUTPUT,
    BWT_OPTION_SAP_ORDERING,
    BWT_OPTION_GENERATE_ENDPOSFILE,
    BWT_OPTION_PAUSE_BETWEEN_CYCLES,
    BWT_OPTION_COUNT // end marker
};

static const string bwtOptionNames[] =
{
    "input format",
    "output format",
    "algorithm",
    "intermediate format",
    "intermediate storage medium",
    "parallel prefetch",
    "parallel processing",
    "process qualities",
    "generate LCP",
    "reverse",
    "single cycle",
    "concatenate output",
    "SAP ordering",
    "generate endPosFile",
    "pause between cycles",
    "" // end marker
};

static const string *bwtOptionPossibleValues[] =
{
    inputFormatLabels,
    outputFormatLabels,
    algorithmOptionLabels,
    intermediateFormatLabels,
    intermediateStorageMediumLabels,
    parallelPrefetchLabels,
    parallelProcessingLabels,
    processQualitiesLabels,
    generateLcpLabels,
    reverseLabels,
    singleCycleLabels,
    concatenateOutputLabels,
    sapOrderingLabels,
    generateEndPosFileLabels,
    pauseBetweenCyclesLabels,
    NULL // end marker
};

} // namespace BeetlBwtParameters


class BwtParameters : public ToolParameters
{
    virtual const string **getOptionPossibleValues() const
    {
        return BeetlBwtParameters::bwtOptionPossibleValues;
    }
public:
    virtual const string getOptionName( const unsigned i ) const
    {
        return BeetlBwtParameters::bwtOptionNames[i];
    }
    //    virtual const string getOptionPossibleValue( const unsigned i, const unsigned j ) const { return bwtOptionPossibleValues[i][j]; }
};


#endif //ifndef BEETL_BWT_PARAMETERS_HH
