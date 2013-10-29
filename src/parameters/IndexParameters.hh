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
        //        addEntry( -1, "BWT in RAM", "--bwt-in-ram", "", "Keep BWT in RAM for faster processing", "", TYPE_SWITCH );

        addDefaultVerbosityAndHelpEntries();
    }

};


#endif //ifndef BEETL_INDEX_PARAMETERS_HH
