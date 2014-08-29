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

#ifndef BEETL_CONVERT_PARAMETERS_HH
#define BEETL_CONVERT_PARAMETERS_HH

#include "libzoo/cli/ToolParameters.hh"

#include <string>

using std::string;


namespace BeetlConvertParameters
{

// options: input format

enum InputFormat
{
    INPUT_FORMAT_FASTA,
    INPUT_FORMAT_FASTQ,
    INPUT_FORMAT_CYC,
    INPUT_FORMAT_SEQ,
    INPUT_FORMAT_BCL,
    INPUT_FORMAT_RUNFOLDER,
    INPUT_FORMAT_BWT_ASCII,
    INPUT_FORMAT_BWT_RLE,
    INPUT_FORMAT_BWT_RLE53,
    INPUT_FORMAT_BWT_RLE_V3,
    INPUT_FORMAT_COUNT
};

static const string inputFormatLabels[] =
{
    "fasta",
    "fastq",
    "cyc",
    "seq",
    "bcl",
    "runFolder",
    "bwt_ascii",
    "bwt_rle",
    "bwt_rle53",
    "bwt_rle_v3",
    "" // end marker
};


// options: output format

enum OutputFormat
{
    OUTPUT_FORMAT_FASTA,
    OUTPUT_FORMAT_FASTQ,
    OUTPUT_FORMAT_CYC,
    OUTPUT_FORMAT_SEQ,
    OUTPUT_FORMAT_BCL,
    OUTPUT_FORMAT_BWT_ASCII,
    OUTPUT_FORMAT_BWT_RLE,
    OUTPUT_FORMAT_BWT_RLE53,
    OUTPUT_FORMAT_BWT_RLE_V2,
    OUTPUT_FORMAT_BWT_RLE_V3,
    OUTPUT_FORMAT_COUNT
};

static const string outputFormatLabels[] =
{
    "fasta",
    "fastq",
    "cyc",
    "seq",
    "bcl",
    "bwt_ascii",
    "bwt_rle",
    "bwt_rle53",
    "bwt_rle_v2",
    "bwt_rle_v3",
    "" // end marker
};

} // namespace BeetlConvertParameters



class ConvertParameters : public ToolParameters
{

public:
    ConvertParameters()
    {
        using namespace BeetlConvertParameters;
        addEntry( -1, "input filename", "--input", "-i", "Input file name or prefix", "", TYPE_STRING | REQUIRED );
        addEntry( -1, "output filename", "--output", "-o", "Output file name or prefix", "", TYPE_STRING | REQUIRED );
        addEntry( -1, "input format", "--input-format", "", "", "detect", TYPE_CHOICE | REQUIRED, inputFormatLabels );
        addEntry( -1, "output format", "--output-format", "", "", "detect", TYPE_CHOICE | REQUIRED, outputFormatLabels );
        addEntry( -1, "sequence length", "--sequence-length", "-l", "If specified, cut the end of longer sequences and pads the start of shorter ones", "", TYPE_INT );
        addEntry( -1, "remove padding", "--remove-padding", "", "For FastQ->FastQ only: Remove 'N' bases from the beginning and end of reads", "", TYPE_SWITCH );
        addEntry( -1, "use missing data from", "--use-missing-data-from", "", "e.g. for FASTA->FASTQ: use the qualities from this file", "", TYPE_STRING );
        addEntry( -1, "extract sequences", "--extract-sequences", "", "Input file containing the sequence numbers to extract (zero-based, one per line)", "", TYPE_STRING );

        addDefaultVerbosityAndHelpEntries();
    }

};


#endif //ifndef BEETL_CONVERT_PARAMETERS_HH
