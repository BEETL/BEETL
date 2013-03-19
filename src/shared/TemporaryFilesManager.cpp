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

#include "TemporaryFilesManager.hh"

#include "Logger.hh"

#include <cstdio>
#include <iostream>

using namespace std;


void TemporaryFilesManager::addFilename( const string &filename )
{
    filenames_.push_back( filename );
}

void TemporaryFilesManager::cleanupAllFiles()
{
    Logger::out( LOG_SHOW_IF_VERBOSE ) << "Removing " << filenames_.size() << " temporary files" << endl;

    for ( unsigned int i = 0; i < filenames_.size(); ++i )
    {
        const string &filename = filenames_[i];
        // cout << "Removing temporary file " << filename << endl;

        if ( remove( filename.c_str() ) != 0 )
            cerr << "TemporaryFilesManager: Error deleting file " << filename << endl;
    }
}
