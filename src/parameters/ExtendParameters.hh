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

#ifndef BEETL_EXTEND_PARAMETERS_HH
#define BEETL_EXTEND_PARAMETERS_HH

#include "libzoo/cli/ToolParameters.hh"


namespace BeetlExtendParameters
{

} // namespace BeetlExtendParameters


class ExtendParameters : public ToolParameters
{
public:
    ExtendParameters()
    {
        using namespace BeetlExtendParameters;
        addEntry( -1, "intervals filename", "--intervals", "-i", "Input file: intervals to extend", "", TYPE_STRING | REQUIRED );
        addEntry( -1, "bwt filename prefix", "--bwt-prefix", "-b", "Input BWT index files prefix", "", TYPE_STRING | REQUIRED );
        addEntry( -1, "sequence numbers output filename", "--output-seqnum", "-o", "Destination file to output sequence numbers", "", TYPE_STRING );
        addEntry( -1, "dollar positions output filename", "--output-dollar-pos", "-p", "Destination file to output BWT positions of dollar signs", "", TYPE_STRING );
        addEntry( -1, "propagate sequence", "--propagate-sequence", "", "Propagate and output sequence with each BWT range (slower)", "", TYPE_SWITCH );

        addDefaultVerbosityAndHelpEntries();
    }

};


#endif //ifndef BEETL_EXTEND_PARAMETERS_HH
