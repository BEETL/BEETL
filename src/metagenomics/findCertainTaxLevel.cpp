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

#include "metaShared.hh"
#include "Krona.hh"
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

    bool inUse_;
    int magnitude_;

    NCBINode()
        : id_( 0 )
        , parentId_( 0 )
        , taxLevel_( 0 )
        , parentLevel_( 0 )
        , inUse_( 0 )
        , magnitude_( 0 )
    {}
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

vector< SequenceInformation > seqs;


vector<SequenceInformation> loadSequenceInfo( string headerFile, string fileCounter );

void findNCBITaxa( string nodesDMP, string mergedDMP, string fileNumFile );

void printAncestors( const int id, const SequenceInformation &seqItem );

map<int, int> giToTaxId;
map<int, string> idToScientificName;
map<int, NCBINode> idToNode;
map<int, vector<NCBINode> > parentToChildren;

int extraTaxId = -2;
map< string, int > extraTaxNameToIdMap;


void parseNamesFile( const string &namesDMP )
{
    ifstream ncbiNames( namesDMP );
    string line;
    vector<string> lineVector;

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
}

void parseGiTotaxId( const string &giToTaxIdFName, const vector<int> &giNumbers )
{
    map<int, bool> giNumbersPresent;
    for ( int giNum : giNumbers )
        giNumbersPresent[giNum] = true;

    ifstream giToTaxFile( giToTaxIdFName );
    string line;
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

    clog << " got " << giToTaxId.size() << " GenBank Ids" << endl;
}

string getPseudoScientificNameFromTag( string tag )
{
    vector<string> tagSplit = splitString( tag, "|" );
    string name = tagSplit[4];
    size_t pos = name.find( ',' );
    if ( pos != string::npos )
        name = name.substr( 0, pos );

    pos = name.find( " complete" );
    if ( pos != string::npos )
        name = name.substr( 0, pos );

    pos = name.find( " main" );
    if ( pos != string::npos )
        name = name.substr( 0, pos );

    pos = name.find( " chromosome" );
    if ( pos != string::npos )
        name = name.substr( 0, pos );

    pos = name.find( " genome" );
    if ( pos != string::npos )
        name = name.substr( 0, pos );

    pos = name.find( " strain " );
    if ( pos != string::npos )
        name = name.erase( pos, 7 );

    if ( name[0] == ' ' )
        name = name.substr( 1 );

    return name;
}

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
    seqs = loadSequenceInfo( headerFile, fileCountFile );
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
    parseGiTotaxId( giToTaxIdFile, gis );
    parseNamesFile( namesFile );

    vector<int> extraGis;
    for ( unsigned int i ( 0 ); i < seqs.size() ; i++ )
    {
        if ( giToTaxId.find( ( int )seqs[i].giNumber_ ) == giToTaxId.end() )
        {
            clog << "# No taxId for giNumber " << seqs[i].giNumber_ << " (seq " << i << ")" << endl;
            // This is due to an incorrect gi number in the NCBI database. Trying to find a correct one by name
            string name = getPseudoScientificNameFromTag( seqs[i].tag_ );

            clog << " Trying identification by name: \"" << seqs[i].tag_ << "\" / \"" << name << "\"" << endl;

            for ( auto & item : idToScientificName )
            {
                if ( !strcasecmp( name.c_str(), item.second.c_str() ) )
                {
                    int taxId = item.first;
                    clog << " Found taxId matching by name: " << taxId << endl;

                    giToTaxId[seqs[i].giNumber_] = taxId;
                    clog << "Faking gi tax number->id (" << giToTaxId.size() << "/" << gis.size() << "): " << seqs[i].giNumber_ << " " << taxId << endl;
                }
            }
        }
    }

    clog << "got " << giToTaxId.size() << " GI IDs " << endl;
    findNCBITaxa( nodesFile, mergedFile, output );

    // output newly created ids
    {
        ofstream ofs( "metaBeetlExtraNames.dmp" );
        for ( const auto & item : extraTaxNameToIdMap )
            ofs << item.second << "\t|\t" << item.first << "\t|\t" << "\t|\t" << "metaBeetl pseudo-scientific name" << "\t|\t" << endl;
    }

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

void markUsefulNodes()
{
    for ( auto & item : giToTaxId )
    {
        int taxId = item.second;
        while ( taxId != idToNode[taxId].parentId_ )
        {
            idToNode[taxId].inUse_ = true;
            ++idToNode[taxId].magnitude_;
            taxId = idToNode[taxId].parentId_;
        }
        ++idToNode[taxId].magnitude_;
        idToNode[taxId].inUse_ = true;
    }
}

void printFullTree_recursive( ofstream &output, const NCBINode &node, const double magnitude = 1.0, const int depth = 0 )
{
    string indent( depth, ' ' );
    //    output << indent << "{ " << id << endl;
    vector<NCBINode> &children = parentToChildren[node.id_];

    // krona output
    output << indent << "<node name=\"" << node.name_ << '(' << node.id_ << ")\">";
    output << indent << "<magnitude><val>" << node.magnitude_ << "</val></magnitude>";
    output << indent << "<rank><val>" << ( depth - 2 ) << '/' << node.taxLevel_ << "</val></rank>" << endl;

    // count useful children
    int usefulChildrenCount = 0;
    for ( auto & childNode : children )
    {
        if ( childNode.inUse_ ) // Only because id 1 has itself as its parent
            ++usefulChildrenCount;
    }

    for ( auto & childNode : children )
    {
        if ( childNode.inUse_ && childNode.id_ != node.id_ ) // Only because id 1 has itself as its parent
            printFullTree_recursive( output, childNode, magnitude / usefulChildrenCount, depth + 1 );
    }
    //    output << indent << "}" << endl;
    output << indent << "</node>" << endl;
}

void printFullTree( )
{
    ofstream treeStream( "tree_krona.html" );
    printKronaHeader( treeStream );

    printFullTree_recursive( treeStream, idToNode[1] );

    printKronaFooter( treeStream );
}

void findNCBITaxa( string nodesDMP, string mergedDMP, string outputInfo )
{
    vector<string> lineVector;
    string line;
    //    map<string, int> nameToId;
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
        node.inUse_ = false;
        idToNode[id] = node;
    }
    clog << "got ncbi names " << endl;

    markUsefulNodes();

    // Create reverse mapping
    for ( auto & node : idToNode )
    {
        parentToChildren[ node.second.parentId_ ].push_back( node.second );
    }

    /*
        //get parent levels, and change the taxLevel under species to strain
        for ( map< int, NCBINode >::iterator it = idToNode.begin(); it != idToNode.end(); ++it )
        {
            NCBINode parent = idToNode[it->second.parentId_];
            if ( parent.taxLevel_ == getTaxonomicLevel( "species" ) )
                it->second.taxLevel_ = getTaxonomicLevel( "strain" );
            else if ( parent.taxLevel_ == getTaxonomicLevel( "strain" ) )
                it->second.taxLevel_ = getTaxonomicLevel( "subsubspecies" );
            it->second.parentLevel_ = parent.taxLevel_;
        }
    */
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
            printAncestors( giToTaxId[gi], seqs[j] );
            cout << endl;
        }
    }

    printFullTree();
}


void printAncestors( const int lowestRankTaxId, const SequenceInformation &seqItem )
{
    int id( lowestRankTaxId );
    clog << "Ancestors for id " << id;
    vector<int> taxIds( taxLevelSize );
    vector<int> taxIds2; // whole taxonomy, including the "no rank" levels
    vector<int> taxIds3; // same as taxIds, plus all the "no rank" levels below "sub-species"
    bool knownRankReached = false;

    if ( id != 0 )
    {
        const NCBINode *n;
        do
        {
            n = &idToNode[id];
            taxIds2.push_back( id );
            if ( n->taxLevel_ < ( int )taxLevelSize )
            {
                taxIds[ n->taxLevel_ ] = id;
                //taxIds3[ n->taxLevel_ ] = id;
                knownRankReached = true;
            }
            else if ( !knownRankReached )
                taxIds3.push_back( id );
            clog << " => id " << id << "(lvl " << n->taxLevel_ << ")";
            id = n->parentId_;
        }
        while ( n->taxLevel_ != 0 );
        clog << endl;
    }

    unsigned int unusedLowestRanks = 0;
    while ( unusedLowestRanks < taxLevelSize && taxIds[ taxLevelSize - unusedLowestRanks - 1] == 0 )
        ++unusedLowestRanks;

    // Check if the lowest rank taxId is a node of the tax tree, in which case we create a sub-item in order to prevent also using it as a leaf
    if ( unusedLowestRanks != taxLevelSize ) // Equality only happens one for an annoying genome that we need to discard
    {
        if ( parentToChildren.find( lowestRankTaxId ) != parentToChildren.end() )
        {
            auto &childrenList = parentToChildren[ lowestRankTaxId ];
            bool atLeastOneInUse = false;
            for ( const auto & child : childrenList )
                atLeastOneInUse |= child.inUse_;
            if ( atLeastOneInUse )
            {
                clog << "TaxId " << lowestRankTaxId << " seems to be both a leaf and a node of the taxonomic tree. Creating a new sub-id to use as a leaf" << endl;
                string pseudoName = getPseudoScientificNameFromTag( seqItem.tag_ );
                if ( extraTaxNameToIdMap.find( pseudoName ) == extraTaxNameToIdMap.end() )
                {
                    extraTaxNameToIdMap[ pseudoName ] = extraTaxId;
                    taxIds3.insert( taxIds3.begin(), extraTaxId );
                    --extraTaxId;
                }
                else
                {
                    taxIds3.insert( taxIds3.begin(), extraTaxNameToIdMap[ pseudoName ] );
                }
                //childrenList.push_back
            }
        }
    }


#ifdef DEBUG_TAX_LEVELS
    cout << " " << taxIds2.size();
    cout << " " << ( taxLevelSize + taxIds3.size() - unusedLowestRanks );

    for ( int k( 0 ) ; k < ( int )taxLevelSize ; ++k )
        cout << " " << taxIds[k];

    cout << " |";
    for ( int k( taxIds2.size() - 1 ) ; k >= 0 ; --k )
        cout << " " << taxIds2[k] << "(" << idToNode[taxIds2[k]].name_ << "|" << idToNode[taxIds2[k]].taxLevel_ << ")";

    cout << " |XXX| ";
    for ( int k( 0 ) ; k < ( int )taxLevelSize - unusedLowestRanks ; ++k )
        cout << " " << taxIds[k] << "(" << idToNode[taxIds[k]].name_ << "|" << idToNode[taxIds[k]].taxLevel_ << ")";

    cout << " |X| ";
    for ( int k( taxIds3.size() - 1 ) ; k >= 0 ; --k )
        cout << " " << taxIds3[k] << "(" << idToNode[taxIds3[k]].name_ << "|" << idToNode[taxIds3[k]].taxLevel_ << ")";
#else

    // Final choice
    for ( uint k( 0 ) ; k < taxLevelSize - unusedLowestRanks ; ++k )
        cout << " " << taxIds[k];

    for ( int k( taxIds3.size() - 1 ) ; k >= 0 ; --k )
        cout << " " << taxIds3[k];

    for ( uint k( 0 ) ; k < unusedLowestRanks - taxIds3.size(); ++k )
        cout << " 0";

#endif // DEBUG_TAX_LEVELS
}
