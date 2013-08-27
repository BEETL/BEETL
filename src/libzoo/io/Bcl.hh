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

#ifndef BCL_HH
#define BCL_HH

#include <fstream>
#include <string>
#include <vector>

using std::string;
using std::vector;


class BclRunFolder
{
public:
    BclRunFolder( const string &runFolder, const string &laneFormat, const string &tileFormat );
    unsigned int getCycleCount();
    bool getNextLaneAndTileNames( string &nextLane, string &nextTile );
    void initReader( const string &lane, const string &tile, const unsigned int firstCycle, const unsigned int lastCycle );
    bool getRead( vector< unsigned int > &bclValues, bool *passFilter = NULL );
    void reportProgress( float step = 0.01 );

private:
    void generateLanesAndTilesList( const string &laneFormat, const string &tileFormat );

    const string runFolder_;
    string lane_;
    string tile_;
    unsigned int cycleCount_;
    unsigned int firstCycle_;
    unsigned int lastCycle_;
    vector< std::istream * > bclFiles_;
    std::ifstream filterFile_;
    unsigned int readCount_;
    unsigned int currentReadNum_;
    unsigned int lastReportedReadNum_;
    vector< std::pair< string, string > > lanesAndTiles_;
    int laneAndTileIndex_;
};

#endif //ifndef BCL_HH
