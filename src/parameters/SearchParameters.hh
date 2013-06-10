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

#ifndef BEETL_SEARCH_PARAMETERS_HH
#define BEETL_SEARCH_PARAMETERS_HH

#include "ToolParameters.hh"

#include <string>


namespace BeetlSearchParameters
{

// Option container

enum SearchOptions
{
    SEARCH_OPTION_COUNT // end marker
};

static const string searchOptionNames[] =
{
    "" // end marker
};

static const string *searchOptionPossibleValues[] =
{
    NULL // end marker
};


} // namespace BeetlSearchParameters


class SearchParameters : public ToolParameters
{
    virtual const string **getOptionPossibleValues() const
    {
        return BeetlSearchParameters::searchOptionPossibleValues;
    }
public:
    virtual const string getOptionName( const unsigned i ) const
    {
        return BeetlSearchParameters::searchOptionNames[i];
    }
};


#endif //ifndef BEETL_SEARCH_PARAMETERS_HH
