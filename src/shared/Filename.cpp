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

#include "Filename.hh"

#include "libzoo/util/TemporaryFilesManager.hh"

using namespace std;


std::string Filename::str() const
{
    return str_;
}

std::string TmpFilename::str() const
{
    if ( fullPathStr_.empty() )
    {
        const string &tempPath = TemporaryFilesManager::get().tempPath_;
        if ( tempPath.empty() )
            fullPathStr_ = str_;
        else
            fullPathStr_ = tempPath + string( "/" ) + str_;
    }
    return fullPathStr_;
}
