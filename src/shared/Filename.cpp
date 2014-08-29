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
