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

#ifndef BEETL_UNBWT_PARAMETERS_HH
#define BEETL_UNBWT_PARAMETERS_HH

#include "libzoo/cli/ToolParameters.hh"

#include <string>


namespace BeetlUnbwtParameters
{

// options: input format

enum InputFormat
{
    INPUT_FORMAT_BWT_ASCII,
    INPUT_FORMAT_COUNT
};

static const string inputFormatLabels[] =
{
    "BWT_ASCII",
    "" // end marker
};


// options: output format

enum OutputFormat
{
    OUTPUT_FORMAT_FASTA,
    OUTPUT_FORMAT_FASTQ,
    OUTPUT_FORMAT_COUNT
};

static const string outputFormatLabels[] =
{
    "fasta",
    "fastq",
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


// options: decode direction off/on

enum DecodeDirection
{
    DECODE_DIRECTION_BACKWARD,
    DECODE_DIRECTION_FORWARD,
    DECODE_DIRECTION_COUNT,
};

static const string decodeDirectionLabels[] =
{
    "backward",
    "forward",
    "" // end marker
};


// options: use vector off/on

enum UseVector
{
    USE_VECTOR_OFF,
    USE_VECTOR_ON,
    USE_VECTOR_COUNT,
};

static const string useVectorLabels[] =
{
    "off",
    "on",
    "" // end marker
};



// Option container

enum UnbwtOptions
{
    //    PARAMETER_PROCESS_QUALITIES,
    PARAMETER_DECODE_DIRECTION,
    PARAMETER_USE_VECTOR,
    PARAMETER_COUNT // end marker
};

} // namespace BeetlUnbwtParameters


class UnbwtParameters : public ToolParameters
{
public:
    UnbwtParameters()
    {
        using namespace BeetlUnbwtParameters;
        addEntry( -1, "input filename prefix", "--input", "-i", "Input file name prefix (without -B0x)", "", TYPE_STRING | REQUIRED );
        addEntry( -1, "output filename", "--output", "-o", "Output file name", "outUnBWT.fasta", TYPE_STRING | REQUIRED );
        addEntry( -1, "input format", "--input-format", "", "Must be:", "detect", TYPE_CHOICE | REQUIRED, inputFormatLabels );
        addEntry( -1, "output format", "--output-format", "", "", "detect", TYPE_CHOICE | REQUIRED, outputFormatLabels );
        addEntry( PARAMETER_DECODE_DIRECTION, "decode direction", "--decode-direction", "-d", "", "backward", TYPE_CHOICE, decodeDirectionLabels );
        addEntry( PARAMETER_USE_VECTOR, "use vector", "--use-vector", "", "", "on", TYPE_CHOICE, useVectorLabels );

        addDefaultVerbosityAndHelpEntries();
    }

};


#endif //ifndef BEETL_UNBWT_PARAMETERS_HH
