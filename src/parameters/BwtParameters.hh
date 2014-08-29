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

#ifndef BEETL_BWT_PARAMETERS_HH
#define BEETL_BWT_PARAMETERS_HH

#include "libzoo/cli/ToolParameters.hh"

#include <string>

using std::string;


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
    OUTPUT_FORMAT_COUNT

    // de-activated
    ,OUTPUT_FORMAT_HUFFMAN
};

static const string outputFormatLabels[] =
{
    "ASCII",
    "RLE",
//    "Huffman",
    "" // end marker
};


// options: algorithm

enum AlgorithmOption
{
    ALGORITHM_BCR,
    ALGORITHM_EXT,
    ALGORITHM_COUNT
};

static const string algorithmLabels[] =
{
    "BCR",
    "ext",
    "" // end marker
};


// options: intermediate format

enum IntermediateFormat
{
    INTERMEDIATE_FORMAT_RLE,
    INTERMEDIATE_FORMAT_ASCII,
    INTERMEDIATE_FORMAT_COUNT

    // de-activated
    ,INTERMEDIATE_FORMAT_MULTIRLE
    ,INTERMEDIATE_FORMAT_HUFFMAN
};

static const string intermediateFormatLabels[] =
{
    "RLE",
    "ASCII",
//    "multiRLE",
//    "Huffman",
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


// options: generate cycle BWT off/pbe/ascii

enum GenerateCycleBwt
{
    GENERATE_CYCLE_BWT_OFF,
    GENERATE_CYCLE_BWT_PBE,
    GENERATE_CYCLE_BWT_ASCII,
    GENERATE_CYCLE_BWT_COUNT,
};

static const string generateCycleBwtLabels[] =
{
    "off",
    "PBE",
    "ASCII",
    "" // end marker
};


// options: paired-reads input

enum PairedReadsInput
{
    PAIRED_READS_INPUT_NONE,
    PAIRED_READS_INPUT_ALL1ALL2,
    PAIRED_READS_INPUT_COUNT
};

static const string pairedReadsInputLabels[] =
{
    "none",
    "all1all2",
    "" // end marker
};


// options: generate cycle qualities off/pbe

enum GenerateCycleQual
{
    GENERATE_CYCLE_QUAL_OFF,
    GENERATE_CYCLE_QUAL_PBE,
    GENERATE_CYCLE_QUAL_COUNT,
};

static const string generateCycleQualLabels[] =
{
    "off",
    "PBE",
    "" // end marker
};



// Option container

enum BwtParameterIds
{
    PARAMETER_UNDEFINED = -1,
    PARAMETER_INPUT_FILENAME = 0,
    PARAMETER_OUTPUT_FILENAME,
    PARAMETER_MEMORY_LIMIT,
    PARAMETER_INPUT_FORMAT,
    PARAMETER_OUTPUT_FORMAT,
    PARAMETER_ALGORITHM,
    PARAMETER_INTERMEDIATE_FORMAT,
    PARAMETER_INTERMEDIATE_STORAGE_MEDIUM,
    PARAMETER_PARALLEL_PREFETCH,
    //    PARAMETER_PARALLEL_PROCESSING,
    PARAMETER_PROCESS_QUALITIES,
    PARAMETER_GENERATE_LCP,
    PARAMETER_ADD_REV_COMP,
    PARAMETER_REVERSE,
    PARAMETER_SUB_SEQUENCE_LENGTH,
    PARAMETER_PAIRED_READS_INPUT,
    PARAMETER_SINGLE_CYCLE,
    PARAMETER_CONCATENATE_OUTPUT,
    PARAMETER_SAP_ORDERING,
    PARAMETER_GENERATE_ENDPOSFILE,
    PARAMETER_GENERATE_CYCLE_BWT,
    PARAMETER_GENERATE_CYCLE_QUAL,
    PARAMETER_PAUSE_BETWEEN_CYCLES,
    PARAMETER_COUNT // end marker
};

} // namespace BeetlBwtParameters



class BwtParameters : public ToolParameters
{

public:
    BwtParameters()
    {
        using namespace BeetlBwtParameters;
        addEntry( PARAMETER_INPUT_FILENAME, "input filename", "--input", "-i", "Input file name or prefix", "", TYPE_STRING | REQUIRED );
        addEntry( PARAMETER_OUTPUT_FILENAME, "output filename", "--output", "-o", "Output file name or prefix", "outBWT", TYPE_STRING | REQUIRED );
        addEntry( PARAMETER_INPUT_FORMAT, "input format", "--input-format", "", "", "detect", TYPE_CHOICE | REQUIRED, inputFormatLabels );
        addEntry( PARAMETER_OUTPUT_FORMAT, "output format", "--output-format", "", "", "rle", TYPE_CHOICE | REQUIRED, outputFormatLabels );
        addEntry( PARAMETER_INTERMEDIATE_FORMAT, "intermediate format", "--intermediate-format", "", "", "", TYPE_CHOICE | REQUIRED | AUTOMATED, intermediateFormatLabels );
        //    addEntry( PARAMETER_INTERMEDIATE_STORAGE_MEDIUM, "intermediate storage medium", "--intermediate-medium", "", "[disk|ram] (multirle->ram, others->disk)", "", TYPE_CHOICE | REQUIRED | AUTOMATED, intermediateStorageMediumLabels );
        addEntry( PARAMETER_ALGORITHM, "algorithm", "--algorithm", "-a", "", "", TYPE_CHOICE | REQUIRED | AUTOMATED, algorithmLabels );
        addEntry( PARAMETER_MEMORY_LIMIT, "memory limit MB", "--memory-limit", "-M", "RAM constraint in MB", "smallest of ulimit -v and /proc/meminfo", TYPE_INT | REQUIRED );

        addEntry( PARAMETER_PROCESS_QUALITIES, "process qualities", "--qualities", "-q", "Ignore/Permute qualities", "ignore", TYPE_CHOICE, processQualitiesLabels );
        addEntry( PARAMETER_CONCATENATE_OUTPUT, "concatenate output", "--concatenate-output", "", "Concatenate BWT files at the end", "", TYPE_SWITCH );
        addEntry( PARAMETER_ADD_REV_COMP, "add reverse complement", "--add-rev-comp", "", "Add reverse complemented sequences", "", TYPE_SWITCH );
        addEntry( PARAMETER_REVERSE, "reverse", "--reverse", "", "Process cycles in reverse order", "", TYPE_SWITCH );
        addEntry( PARAMETER_SUB_SEQUENCE_LENGTH, "sub-sequence length", "--sub-sequence-length", "", "Split sequences into two sub-sequences. Useful for paired reads", "", TYPE_INT );
        addEntry( PARAMETER_PAIRED_READS_INPUT, "paired-reads input", "--paired-reads-input", "", "If your input file contains paired reads", "none", TYPE_CHOICE, pairedReadsInputLabels );
        addEntry( PARAMETER_SAP_ORDERING, "SAP ordering", "--sap-ordering", "", "Use SAP ordering (see SAP note below)", "", TYPE_SWITCH );
        addEntry( PARAMETER_GENERATE_ENDPOSFILE, "generate endPosFile", "--generate-end-pos-file", "", "Generate mapping between BWT '$' signs and sequence numbers", "", TYPE_SWITCH );
        addEntry( PARAMETER_GENERATE_LCP, "generate LCP", "--generate-lcp", "", "Generate Longest Common Prefix lengths (see LCP note below)", "", TYPE_SWITCH );
        addEntry( PARAMETER_GENERATE_CYCLE_BWT, "generate cycle BWT", "--cycle-bwt", "", "PBE=Generate cycle-by-cycle BWT with prediction-based encoding", "off", TYPE_CHOICE, generateCycleBwtLabels );
        addEntry( PARAMETER_GENERATE_CYCLE_QUAL, "generate cycle qualities", "--cycle-qual", "", "PBE=Generate cycle-by-cycle qualities zeroed at correctly-predicted bases", "off", TYPE_CHOICE, generateCycleQualLabels );
#ifdef _OPENMP
        addEntry( PARAMETER_PARALLEL_PREFETCH, "parallel prefetch", "--no-parallel-prefetch", "", "Disable parallel prefetch of cycle files", "", TYPE_SWITCH | AUTOMATED );
        //    addEntry( PARAMETER_PARALLEL_PROCESSING, "parallel processing", "--no-parallel-processing", "", "Disable parallel processing by letter", "", TYPE_SWITCH | AUTOMATED, parallelProcessingLabels );
#endif //ifdef _OPENMP
        //    addEntry( PARAMETER_, "", " --hw-constraints         File describing hardware constraints for speed estimates", "", TYPE_STRING );
        addEntry( PARAMETER_PAUSE_BETWEEN_CYCLES, "pause between cycles", "--pause-between-cycles", "", "Wait for a key press after each cycle", "", TYPE_SWITCH );

        addDefaultVerbosityAndHelpEntries();
    }

};


#endif //ifndef BEETL_BWT_PARAMETERS_HH
