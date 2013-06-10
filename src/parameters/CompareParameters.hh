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


namespace BeetlCompareParameters
{

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
    COMPARE_OPTION_MODE,
    COMPARE_OPTION_REPORT_MINLENGTH,
    COMPARE_OPTION_COUNT // end marker
};

static const string compareOptionNames[] =
{
    "mode",
    "report minlength",
    "" // end marker
};

static const string *compareOptionPossibleValues[] =
{
    modeLabels,
    reportMinlengthLabels,
    NULL // end marker
};


} // namespace BeetlCompareParameters


class CompareParameters : public ToolParameters
{
    virtual const string **getOptionPossibleValues() const
    {
        return BeetlCompareParameters::compareOptionPossibleValues;
    }
public:
    virtual const string getOptionName( const unsigned i ) const
    {
        return BeetlCompareParameters::compareOptionNames[i];
    }
};


#endif //ifndef BEETL_COMPARE_PARAMETERS_HH
