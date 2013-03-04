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

#ifndef TOOL_HH
#define TOOL_HH


enum CompressionFormatType
{
    compressionASCII,
    compressionRunLength,
    compressionIncrementalRunLength,
    compressionHuffman
};


class Algorithm
{
public:
    virtual void run( void )  = 0; // run method, must be implemented by all tools

    virtual ~Algorithm() {}

};




#endif /* TOOL_HH */

