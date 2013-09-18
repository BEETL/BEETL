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

#ifndef ALIGNMENT_PARAMETERS_HH
#define ALIGNMENT_PARAMETERS_HH

#include "libzoo/cli/ToolParameters.hh"

#include <string>

using std::string;


// options: input format

enum ReadsFormat
{
    READS_FORMAT_FASTA = 0,
    READS_FORMAT_FASTQ = 1
};

static const string formatLabels[] =
{
    "fasta",
    "fastq",
    "" // end marker
};

enum AlignmentType
{
    ALIGNMENT_TYPE_SW = 0,
    ALIGNMENT_TYPE_NO_INDELS = 1,
    ALIGNMENT_TYPE_STITCH = 2//where the correction strings are just stitched on...
};

static const string alignmentTypeLabels[] =
{
    "smith-waterman",
    "no-indels",
    "stitch",
    "" // end marker
};

// Option container

enum AlignmentParameterIds
{
    PARAMETER_UNDEFINED = -1,
    PARAMETER_INPUT_READS_FILE = 0,
    PARAMETER_INPUT_READS_FORMAT,
    PARAMETER_OUTPUT_READS_FORMAT,
    PARAMETER_INPUT_CORRECTIONS_FILE,
    PARAMETER_CORRECTED_READS_FILE,
    PARAMETER_ALIGNMENT_TYPE,
    PARAMETER_SW_INSERTION_PENALTY,
    PARAMETER_SW_DELETION_PENALTY,
    PARAMETER_SW_MISMATCH_PENALTY,
    PARAMETER_CORRECTION_QUALITY,
    PARAMETER_MIN_WITNESS_LENGTH,
    PARAMETER_TRIM_CORRECTED_READS
};

class AlignmentParameters : public ToolParameters
{
public:
    AlignmentParameters()
    {
        addEntry( PARAMETER_INPUT_READS_FILE , "input reads file", "--input-reads", "-i", "Original reads file", "", TYPE_STRING | REQUIRED );
        addEntry( PARAMETER_INPUT_READS_FORMAT , "input reads format", "--input-reads-format", "", "Original reads format", "", TYPE_CHOICE | REQUIRED, formatLabels );
        addEntry( PARAMETER_OUTPUT_READS_FORMAT , "output reads format", "--output-reads-format", "", "Output reads format", "", TYPE_CHOICE | REQUIRED, formatLabels );
        addEntry( PARAMETER_INPUT_CORRECTIONS_FILE , "input corrections file", "--input-corrections-file", "-c", "Input corrections file", "", TYPE_STRING | REQUIRED );
        addEntry( PARAMETER_ALIGNMENT_TYPE , "alignment type", "--alignment-type", "-a", "Alignment type", "", TYPE_CHOICE | REQUIRED, alignmentTypeLabels );
        addEntry( PARAMETER_CORRECTED_READS_FILE, "corrected reads output file", "--corrected-reads-file", "-o", "File to which corrected reads are written", "", TYPE_STRING | REQUIRED );
        addEntry( PARAMETER_SW_MISMATCH_PENALTY, "mismatch penalty", "--mismatch-penalty", "-m", "Mismatch penalty for Smith-Waterman alignment", "", TYPE_INT );
        addEntry( PARAMETER_SW_DELETION_PENALTY, "deletion penalty", "--deletion-penalty", "-d", "Deletion penalty for Smith-Waterman alignment", "", TYPE_INT );
        addEntry( PARAMETER_SW_INSERTION_PENALTY, "insertion penalty", "--insertion-penalty" , "-n", "Insertion penalty for Smith-Waterman alignment", "", TYPE_INT );
        addEntry( PARAMETER_CORRECTION_QUALITY, "correction quality", "--correction-quality" , "-q", "Correction string letter quality character", "", TYPE_STRING );
        addEntry( PARAMETER_MIN_WITNESS_LENGTH, "min witness length", "--min-witness-length" , "", "Minimum length of witness to use correction", "", TYPE_INT );
        addEntry( PARAMETER_TRIM_CORRECTED_READS, "trim corrected reads", "--trim" , "-t", "Trim corrected reads and quality strings to length of originals", "", TYPE_SWITCH );
        addDefaultVerbosityAndHelpEntries();
    }
};


#endif //ifndef ALIGNMENT_PARAMETERS_HH
