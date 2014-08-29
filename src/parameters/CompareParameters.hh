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

#ifndef BEETL_COMPARE_PARAMETERS_HH
#define BEETL_COMPARE_PARAMETERS_HH

#include "libzoo/cli/ToolParameters.hh"

#include <string>


namespace BeetlCompareParameters
{

// options: mode splice/reference/metagenomics

enum Mode
{
    MODE_TUMOUR_NORMAL,
    MODE_SPLICE,
    MODE_REFERENCE,
    MODE_METAGENOMICS,
    MODE_COUNT,
};

static const string modeLabels[] =
{
    "tumour-normal",
    "splice",
    "reference",
    "metagenomics",
    "" // end marker
};


// options: report minlength off/on

enum ReportMinlength
{
    REPORT_MINLENGTH_OFF,
    REPORT_MINLENGTH_ON,
    REPORT_MINLENGTH_COUNT
};

static const string reportMinlengthLabels[] =
{
    "off",
    "on",
    "" // end marker
};


// options: input format

enum InputFormat
{
    INPUT_FORMAT_BWT_ASCII,
    INPUT_FORMAT_BWT_RLE,
    INPUT_FORMAT_DETECT,
    INPUT_FORMAT_COUNT
};

static const string inputFormatLabels[] =
{
    "bwt_ascii",
    "bwt_rle",
    "detect",
    "" // end marker
};



// Option container

enum CompareOptions
{
    COMPARE_OPTION_MODE,
    COMPARE_OPTION_REPORT_MINLENGTH,
    COMPARE_OPTION_COUNT
};


} // namespace BeetlCompareParameters


class CompareParameters : public ToolParameters
{
public:
    CompareParameters()
    {
        using namespace BeetlCompareParameters;
        addEntry( COMPARE_OPTION_MODE, "mode", "--mode", "-m", "See \"Mode\" note below", "", TYPE_CHOICE | REQUIRED, modeLabels );
        addEntry( -1, "input setA", "--inputA", "-a", "Input filename prefix for Set A (such as \"prefix-B0[0-6]\" are Set A's BWT files)", "", TYPE_STRING | REQUIRED );
        addEntry( -1, "input setB", "--inputB", "-b", "Input filename prefix for Set B (such as \"prefix-B0[0-6]\" are Set B's BWT files)", "", TYPE_STRING | REQUIRED );
        addEntry( -1, "output directory", "--output", "-o", "Output directory", "BeetlCompareOutput", TYPE_STRING | REQUIRED );
        addEntry( -1, "max length", "--max-length", "-k", "Maximal k-mer length (number of analysis cycles)", "100", TYPE_INT );
        addEntry( -1, "min occ", "--min-occ", "-n", "Minimum number of occurrences (coverage)", "2", TYPE_INT );
        addEntry( -1, "inputA format", "--inputA-format", "", "", "detect", TYPE_CHOICE | REQUIRED, inputFormatLabels );
        addEntry( -1, "inputB format", "--inputB-format", "", "", "detect", TYPE_CHOICE | REQUIRED, inputFormatLabels );
        addEntry( -1, "subset", "--subset", "", "Restrict computation to this suffix - Used for distributed computing", "", TYPE_STRING );
        addEntry( -1, "4-way distributed", "--4-way", "", "4-way distributed process number (0-3)", "", TYPE_INT );
        addEntry( -1, "generate seq num A", "--generate-seq-numA", "", "Propagate breakpoints to output read numbers. Requires {inputA}-end-pos file", "", TYPE_SWITCH );
        addEntry( -1, "generate seq num B", "--generate-seq-numB", "", "Propagate breakpoints to output read numbers. Requires {inputB}-end-pos file", "", TYPE_SWITCH );
        addEntry( -1, "memory limit MB", "--memory-limit", "-M", "RAM constraint in MB", "smallest of ulimit -v and /proc/meminfo", TYPE_INT | REQUIRED );
        addEntry( -1, "no comparison skip", "--no-comparison-skip", "", "Don't skip already processed comparisons (slower, but smoother output)", "", TYPE_SWITCH );
        addEntry( -1, "pause between cycles", "--pause-between-cycles", "", "Wait for a key press after each cycle", "", TYPE_SWITCH );
        addEntry( -1, "BWT in RAM", "--bwt-in-ram", "", "Keep BWT in RAM for faster processing", "", TYPE_SWITCH );
        addEntry( -1, "propagate sequence", "--propagate-sequence", "", "Propagate and output sequence with each BWT range (slower)", "", TYPE_SWITCH );

        //        addEntry( -1, "setB metadata", "--genome-metadata", "-c", "For Metagenomics mode only: Input filename \"extended\" prefix for Set B's metadata (for files \"prefix[0-6]\")", "${inputB}-C0", TYPE_STRING );
        addEntry( -1, "taxonomy", "--taxonomy", "-t", "For Metagenomics mode only: Input filename for Set B's taxonomy information", "", TYPE_STRING );
        addEntry( -1, "min kmer length", "--min-kmer-length", "-w", "For Metagenomics mode only: Minimum k-mer length", "50", TYPE_INT );
        addEntry( COMPARE_OPTION_REPORT_MINLENGTH, "report min length", "--report-min-length", "-d", "For Metagenomics mode only: Report the minimal needed word length for the different taxa in the database", "off", TYPE_CHOICE, reportMinlengthLabels );
        addEntry( -1, "mmap C files", "--mmap-c-files", "", "Memory-mapping -C0* files may lead to performance improvements", "", TYPE_SWITCH );

        addDefaultVerbosityAndHelpEntries();
    }

};


#endif //ifndef BEETL_COMPARE_PARAMETERS_HH
