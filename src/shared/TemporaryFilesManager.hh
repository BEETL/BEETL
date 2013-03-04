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

#ifndef TEMPORARY_FILES_MANAGER_HH
#define TEMPORARY_FILES_MANAGER_HH

#include <string>
#include <vector>

using std::string;
using std::vector;

class TemporaryFilesManager
{
private:
    TemporaryFilesManager() {}
    TemporaryFilesManager( TemporaryFilesManager const & ); // no impl to avoid copies of singleton
    void operator=( TemporaryFilesManager const & ); // no impl to avoid copies of singleton

public:
    static TemporaryFilesManager &get()
    {
        static TemporaryFilesManager singleton;
        return singleton;
    }

    void addFilename( const string &filename );
    void cleanupAllFiles();

private:
    vector<string> filenames_;
};

#endif // TEMPORARY_FILES_HH
