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

#include "../shared/Tools.hh"

#include<cassert>
#include<cstdlib>
#include<cstring>
#include<fstream>
#include<iostream>
#include<map>
#include<sstream>
#include<stdint.h>
#include<string>
#include<vector>

using namespace std;


struct NCBINode
{
    string name_;
    int id_;
    int parentId_;
    int taxLevel_;
    int parentLevel_;
};

struct SequenceInformation
{
    unsigned short fileNum_;
    string tag_;
    string fileName_;
    uint64_t giNumber_;
    string giString_;

    SequenceInformation() : fileNum_( 0 ), giNumber_( 0 ) {}
};


vector<SequenceInformation> loadSequenceInfo( string headerFile, string fileCounter );

void findNCBITaxa( vector<SequenceInformation> &seqs, string namesDMP, string nodesDMP, string mergedDMP, string fileNumFile );

map<int, int> getGiTotaxId( string gitoTaxId, vector<int> giNumbers );

int getTaxonomicLevel( string s );

void printAncestors( map<int, NCBINode > &idToNode, int id );

//string levelNames[] = {"superkingdom","phylum" ,"class", "order", "family", "genus", "species"};

map<int, int> giToTaxId;


int main( int argc, char **argv )
{
    string namesFile;
    string nodesFile;
    string mergedFile;
    string headerFile;
    string output;
    string fileCountFile;
    string giToTaxIdFile;

    if ( argc < 5 )
        cout << "-nA namesfile, -nO nodesFile, [-nM mergedFile,] -nG giToTaxIdFile, -h headerFile, -f fileCountFile " << endl;
    for ( int i( 0 ) ; i < argc ; i++ )
    {
        if ( strcmp( argv[i], "-nA" ) == 0 )
            namesFile = argv[i + 1];
        else if ( strcmp( argv[i] , "-nO" ) == 0 )
            nodesFile = argv[i + 1];
        else if ( strcmp( argv[i] , "-nM" ) == 0 )
            mergedFile = argv[i + 1];
        else if ( strcmp( argv[i], "-h" ) == 0 )
            headerFile = argv[i + 1];
        else if ( strcmp( argv[i], "-f" ) == 0 )
            fileCountFile = argv[i + 1];
        else if ( strcmp( argv[i], "-nG" ) == 0 )
            giToTaxIdFile = argv[i + 1];
    }
    vector <SequenceInformation > seqs = loadSequenceInfo( headerFile, fileCountFile );
    vector<int> gis;
    gis.push_back( seqs[0].giNumber_ );
    for ( unsigned int i( 1 ); i < seqs.size(); i++ )
    {
        bool found = false;
        for ( unsigned int j( 0 ); j < gis.size(); j++ )
        {
            if ( gis[j] == ( int )seqs[i].giNumber_ )
            {
                found = true;
                break;
            }
        }
        if ( ! found )
        {
            gis.push_back( seqs[i].giNumber_ );
            //            clog << seqs[i].giNumber_ << "|" ;
        }
    }
    //    clog << endl;
    clog << "Gi Id size " << gis.size() << endl;
    giToTaxId = getGiTotaxId( giToTaxIdFile, gis );
    for ( unsigned int i ( 0 ); i < seqs.size() ; i++ )
    {
        if ( giToTaxId.find( ( int )seqs[i].giNumber_ ) == giToTaxId.end() )
            cout << "# No taxId for giNumber " << seqs[i].giNumber_ << " (seq " << i << ")" << endl;
    }

    clog << "got " << giToTaxId.size() << " GI IDs " << endl;
    findNCBITaxa( seqs, namesFile, nodesFile, mergedFile, output );
    return 0;
}

vector<SequenceInformation> loadSequenceInfo( string headerFile, string fileCounter )
{
    vector<SequenceInformation> seqsInfo;
    ifstream fileCount( fileCounter );
    string line;
    string fName;
    while  ( getline( fileCount, line ) )
    {
        SequenceInformation oneSeq;
        unsigned short fileCNum = ( unsigned short ) strtoul( line.substr( 0, line.find( "," ) ).c_str(), NULL, 0 );
        oneSeq.fileNum_ = fileCNum;
        fName = line.substr( line.find( "," ) + 1 );
        oneSeq.fileName_ = fName;
        seqsInfo.push_back( oneSeq );
    }
    ifstream head( headerFile );
    int tagCount( 0 );

    while ( getline( head, line ) )
    {
        string fileCount = line.substr( 0, line.find( "," ) );
        string tag = line.substr( line.find( "," ) + 1 );
        string f = "G_" + fileCount;
        string fr = "G_" + fileCount + "_rev";
        int revCount( 0 );
        int foundFile( 0 );

        for ( unsigned int i ( 0 ); i < seqsInfo.size(); i++ )
        {
            //      cout << seqsInfo[i].fileName_ << " f " << f << " fr "  <<fr <<endl;
            if ( seqsInfo[i].fileName_.compare( f ) == 0 || seqsInfo[i].fileName_.compare( fr ) == 0 )
            {
                revCount++;
                //	cout << seqsInfo[i].fileName_ << " f " << f << " fr "  <<fr <<endl;
                if ( i != seqsInfo[i].fileNum_ )
                    cerr << "wrong file number " << endl;
                //1,gi|15604717|ref|NC_000117.1| Chlamydia trachomatis D/UW-3/CX, complete genome
                seqsInfo[i].tag_ = tag;
                vector<string> tagSplit = splitString( line, "|" );
                seqsInfo[i].giNumber_ =  atol( tagSplit[1].c_str() );

                seqsInfo[i].giString_ = tagSplit[1];
                tagCount++;
                ++foundFile;
            }
        }
        if ( foundFile != 2 )
        {
            cerr << "Found " << foundFile << " File for " << line << endl;
            cerr << "file should be " << f << " or " << fr << endl;
        }
    }
    clog << "TagCount " << tagCount << endl;
    clog << "Got Sequenceinformation " << seqsInfo.size() << endl;
    return seqsInfo;
}

map<int, int> getGiTotaxId( string giToTaxIdFName, vector<int> giNumbers )
{
    map<int, bool> giNumbersPresent;
    for ( int giNum : giNumbers )
        giNumbersPresent[giNum] = true;

    ifstream giToTaxFile( giToTaxIdFName );
    string line;
    map<int, int> giToTaxId;
    clog << "GI numbers count: " << giNumbers.size() << endl;

    while ( getline( giToTaxFile, line ) )
    {
        istringstream iss( line );
        int giNumber = -1;// = atoi( lineVector[0].c_str() );
        int taxId = -1; //atoi( lineVector[1].c_str() );

        if ( iss >> giNumber >> taxId )
        {
            assert( taxId != -1 );
            if ( giNumbersPresent.find( giNumber ) != giNumbersPresent.end() )
            {
                giToTaxId[giNumber] = taxId;
                clog << "Found gi tax number->id (" << giToTaxId.size() << "/" << giNumbers.size() << "): " << giNumber << " " << taxId << endl;
            }
        }
    }

    clog << " got all " << giToTaxId.size() << " GIs" << endl;

    return giToTaxId;
}

void findNCBITaxa( vector<SequenceInformation> &seqs, string namesDMP, string nodesDMP, string mergedDMP, string outputInfo )
{
    ifstream ncbiNames( namesDMP );
    vector<string> lineVector;
    string line;
    //    map<string, int> nameToId;
    map<int, string> idToScientificName;
    map<int, NCBINode> idToNode;
    vector<  map<int, vector<unsigned short> > > taxLevelToFileCount;
    taxLevelToFileCount.resize( 7 );

    while ( getline( ncbiNames, line ) )
    {
        lineVector = splitString( line , "\t|\t" );
        string name = lineVector[1];
        //    cout << "ncbiName >"<< name <<"<" <<endl;
        int id = atoi( lineVector[0].c_str() );
        //        nameToId[name] = id;
        if ( line.find( "scientific" ) != string::npos )
            idToScientificName[id] = name;
    }
    ifstream ncbiNodes( nodesDMP );
    NCBINode node;
    while ( getline( ncbiNodes, line ) )
    {
        lineVector = splitString( line, "\t|\t" );

        int id = atoi( lineVector[0].c_str() );
        int parentId = atoi( lineVector[1].c_str() );
        node.name_ = idToScientificName[id];
        node.id_ = id;
        if ( lineVector.size() > 2 )
        {
            int taxLevel =  getTaxonomicLevel( lineVector[2] );
            node.taxLevel_ = taxLevel;
        }
        else
            cerr << line << endl;
        node.parentId_ = parentId;
        if ( idToNode.find( id ) != idToNode.end() )
            cerr << "already has id " << id << " " << idToScientificName[id] << endl;
        idToNode[id] = node;
    }
    clog << "got ncbi names " << endl;
    //get parent levels, and change the taxLevel under species to strain
    for ( map< int, NCBINode >::iterator it = idToNode.begin(); it != idToNode.end(); ++it )
    {
        NCBINode parent = idToNode[it->second.parentId_];
        if ( parent.taxLevel_ == getTaxonomicLevel( "species" ) )
            it->second.taxLevel_ = getTaxonomicLevel( "strain" );
        it->second.parentLevel_ = parent.taxLevel_;
    }
    clog << "node to id " << idToNode.size() << endl;

    /*  ifstream mergedStream( mergedDMP.c_str(), ios::in);
    while(mergedStream.good()) {
      getline(mergedStream, line);
      lineVector = split(line,"\t|\t");

      int oldId = atoi(lineVector[0].c_str());
      int newId = atoi(lineVector[1].c_str());
      cerr << oldId << " " << newId <<endl;

      NCBINode oldNode = idToNode[oldId];
      cerr << "name " << oldNode.name_ <<endl;
      idToNode[newId] = oldNode;
      }*/

    clog << "after merged out " << idToNode.size() << endl;
    clog << "seq size " << seqs.size() << endl;
    for ( unsigned int j ( 0 ) ; j < seqs.size(); j++ )
    {
        unsigned short fileNum = seqs[j].fileNum_;
        //1,gi|15604717|ref|NC_000117.1| Chlamydia trachomatis D/UW-3/CX, complete genome
        int gi = seqs[j].giNumber_;
        clog << "fileNum " << fileNum << " (" << j << "/" << seqs.size() << "): gi=" << gi << ", taxId=" << giToTaxId[gi] << endl;
        //        if ( giToTaxId[gi] != 0 )
        {
            cout << fileNum;
            printAncestors( idToNode, giToTaxId[gi] );
            cout << endl;
        }
    }
}


void printAncestors( map<int, NCBINode> &idToNode, int id )
{
    clog << "Ancestors for id " << id;
    vector<int> taxIds( 8 );

    if ( id != 0 )
    {
        const NCBINode *n;
        do
        {
            n = &idToNode[id];
            if ( n->taxLevel_ < 8 )
            {
                taxIds[ n->taxLevel_ ] = id;
            }
            clog << " => id " << id << "(lvl " << n->taxLevel_ << ")";
            id = n->parentId_;
        }
        while ( n->taxLevel_ != 0 );
        clog << endl;
    }

    for ( int k( 0 ) ; k < 8 ; ++k )
        cout << " " << taxIds[k];
}


int getTaxonomicLevel( string s )
{
    //string levelNames[] = {"superkingdom","phylum" ,"class", "order", "family", "genus", "species", "strain"};
    int level;
    if ( s.compare( "superkingdom" ) == 0 )
        level = 0;
    else if ( s.compare( "phylum" ) == 0 )
        level = 1;
    else if ( s.compare( "class" ) == 0 )
        level = 2;
    else if ( s.compare( "order" ) == 0 )
        level = 3;
    else if ( s.compare( "family" ) == 0 )
        level = 4;
    else if ( s.compare( "genus" ) == 0 )
        level = 5;
    else if ( s.compare( "species" ) == 0 )
        level = 6;
    else if ( s.compare( "strain" ) == 0 )
        level = 7;
    else
        level = 10;
    return level;
}

