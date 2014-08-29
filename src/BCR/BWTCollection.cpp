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


/******************************************************************************
* It assumes that:
* -  the length of each sequence is 100.
* -  the alphabet is $ACGNTZ
******************************************************************************/

#include "BWTCollection.hh"

#include "BCRexternalBWT.hh"


namespace SXSI
{
/**
 * Init bwt collection
 *
 * See BCRexternalBWT.h for more details.
 */
BWTCollection *BWTCollection::InitBWTCollection( const string &file1, const string &fileOut, const int mode, const CompressionFormatType outputCompression )
{
    BWTCollection *result =
        new BCRexternalBWT( file1, fileOut, mode, outputCompression );
    return result;
}
}
