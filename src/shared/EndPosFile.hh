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

#ifndef INCLUDED_ENDPOSFILE
#define INCLUDED_ENDPOSFILE

#include "Types.hh"

#include <fstream>


//SequenceNumber EndPosFile_convertDollarNumToSequenceNum( const SequenceNumber dollarNum );
class EndPosFile
{
public:
    EndPosFile( const string &bwtFilenamePrefix );
    SequenceNumber convertDollarNumToSequenceNum( const SequenceNumber dollarNum );

private:
    std::ifstream file_;

    SequenceNumber sequenceGroupCount_;
    uint8_t sequenceCountInGroup_;
    uint8_t hasRevComp_;
    SequenceNumber dollarSignCount_;
};


#endif // INCLUDED_ENDPOSFILE
