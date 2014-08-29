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

