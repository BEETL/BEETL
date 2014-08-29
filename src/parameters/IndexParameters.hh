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

#ifndef BEETL_INDEX_PARAMETERS_HH
#define BEETL_INDEX_PARAMETERS_HH

#include "libzoo/cli/ToolParameters.hh"

#include <string>


namespace BeetlIndexParameters
{

// Option container

enum SearchOptions
{
    SEARCH_OPTION_COUNT // end marker
};

} // namespace BeetlIndexParameters


class IndexParameters : public ToolParameters
{
public:
    IndexParameters()
    {
        using namespace BeetlIndexParameters;
        addEntry( -1, "input", "--input", "-i", "Input filename prefix (i.e. BWT files are \"prefix-B0[0-6]\")", "", TYPE_STRING | REQUIRED );
        addEntry( -1, "block size", "--block-size", "-b", "Interval between index points (smaller=faster but more RAM)", "", TYPE_INT );
        addEntry( -1, "force", "--force", "-f", "Overwrite any existing index files", "", TYPE_SWITCH );


        //        addEntry( -1, "output", "--output", "-o", "Output filename", "searchedKmers_positions", TYPE_STRING | REQUIRED );
        //        addEntry( -1, "pause between cycles", "--pause-between-cycles", "", "Wait for a key press after each cycle", "", TYPE_SWITCH );

        addDefaultVerbosityAndHelpEntries();
    }

};


#endif //ifndef BEETL_INDEX_PARAMETERS_HH
