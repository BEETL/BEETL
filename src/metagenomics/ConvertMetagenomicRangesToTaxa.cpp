#include "../shared/Alphabet.hh"
#include "../shared/Tools.hh"
#include "../shared/Types.hh"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace std;


#include "metaShared.hh"
#include "Krona.hh"


struct InputLine //metaBeetl.out.cycle23.subset2.5:BKPT+SGLT       5       T       22      0:0:1:0:0:0     0:0:1:0:0:0     36413   743151552       1       1
{
    string header;
    uint8_t pileNum;
    char pileChar;
    uint16_t cycle;
    string countA, countB;
    LetterNumber posA, posB;
    LetterNumber numA, numB;

    friend std::istream &operator>>( std::istream &is, InputLine &obj )
    {
        is
                >> obj.header
                >> obj.pileNum
                >> obj.pileChar
                >> obj.cycle
                >> obj.countA
                >> obj.countB
                >> obj.posA
                >> obj.posB
                >> obj.numA
                >> obj.numB;
        if ( obj.pileNum >= '0' && obj.pileNum <= '6' )
            obj.pileNum -= '0';
        return is;
    }

    friend std::ostream &operator<<( std::ostream &os, const InputLine &obj )
    {
        return os
               << "{ " << obj.header
               << ", " << ( obj.pileNum + '0' )
               << ", " << obj.pileChar
               << ", " << obj.cycle
               << ", " << obj.countA
               << ", " << obj.countB
               << ", " << obj.posA
               << ", " << obj.posB
               << ", " << obj.numA
               << ", " << obj.numB
               << " }";
    }
};

struct ShortInputLineWithId
{
    uint32_t id;
    //    uint8_t pileNum;
    LetterNumber posB;
    LetterNumber numB;
    ShortInputLineWithId()
        : id( 0 )
        , posB( 0 )
        , numB( 0 )
    {}

    ShortInputLineWithId( const uint32_t id, const InputLine &inputLine )
        : id( id )
        //        , pileNum( inputLine.pileNum )
        , posB( inputLine.posB )
        , numB( inputLine.numB )
    {}

    friend std::ostream &operator<<( std::ostream &os, const ShortInputLineWithId &obj )
    {
        return os
               << "{ " << obj.id
               //            << ", " << obj.pileNum
               << ", " << obj.posB
               << ", " << obj.numB
               << " }";
    }

    bool operator<( const ShortInputLineWithId &rhs ) const
    {
        //        return posB < rhs.posB;
        return posB < rhs.posB || ( posB == rhs.posB && numB > rhs.numB );
        //        return pileNum < rhs.pileNum || ( pileNum == rhs.pileNum && posB < rhs.posB );
    }
};

vector< vector< int> > fileNumToTaxIds;
void loadFileNumToTaxIds( const string &taxIdFilename )
{
    ifstream taxas( taxIdFilename );

    string line;
    //each line should contain the fileNumber followed by the taxIds split up with one space character
    while ( getline( taxas, line ) )
    {
        if ( line.empty() || line[0] == '#' )
            continue;

        istringstream iss( line );

        unsigned int fileNum;
        int taxId, lastTaxId = 0;
        vector<int> taxIds;

        iss >> fileNum;
        if ( fileNum != fileNumToTaxIds.size() )
            cout << "Wrong filenumber " << fileNum << " " << fileNumToTaxIds.size() << endl;

        while ( iss >> taxId )
        {
            if ( taxId )
                lastTaxId = taxId;
            else
                taxId = lastTaxId;

            taxIds.push_back( taxId );
        }
        //test if all TaxIds were found
        if ( taxIds.size() < taxLevelSize )
            cerr << "Tax Ids don't have enough taxonomic Information. Only " << taxIds.size() << " could be found (" << line << ")" << endl
                 << "Will add unknown taxa until size is right" << endl;
        else if ( taxIds.size() > taxLevelSize )
            cerr << "Tax Ids have too much taxonomic information (" << line << ")" << endl
                 << "Please note, that the taxonomic information about one file should be shown as: " << endl
                 << "FileNumber Superkingdom Phylum Order Family Genus Species Strain " << endl;
        taxIds.resize( taxLevelSize );

        fileNumToTaxIds.push_back( taxIds );
    }
}

void restrictCommonTaxaToIncludeFileNum( vector<int> &commonTaxa, const MetagFileNumRefType fileNum )
// 0=undefined, -1=not match
{
    vector< int > &taxIds = fileNumToTaxIds[fileNum];
    for ( unsigned int i = 0; i < taxLevelSize; ++i )
    {
        switch ( commonTaxa[i] )
        {
            case -1:
                break;
            case 0:
                commonTaxa[i] = taxIds[i];
                break;
            default:
                if ( commonTaxa[i] != taxIds[i] )
                    commonTaxa[i] = -1;
                break;
        }
    }
}


int main( int argc, char **argv )
{
    if ( argc != 3 )
    {
        cerr << "Usage: " << argv[0] << " fileNumToTaxTree databaseBwtPrefix" << endl;
        exit( -1 );
    }

    loadFileNumToTaxIds( argv[1] );

    InputLine inputLine;
    //vector<InputLine> inputLines;
    vector<ShortInputLineWithId> shortInputLines[alphabetSize];
    uint32_t id = 0;
    while ( cin >> inputLine )
    {
        //cout << inputLine << endl;
        //inputLines.push_back( inputLine );
        shortInputLines[inputLine.pileNum].push_back( ShortInputLineWithId( id++, inputLine ) );
    }

    map<int, int> lcaCounts;

    for ( int i = 0; i < alphabetSize; ++i )
    {
        ShortInputLineWithId lastItem;

        if ( !shortInputLines[i].empty() )
        {
            sort( shortInputLines[i].begin(), shortInputLines[i].end() );

            ostringstream cFilename;
            cFilename << argv[2] << "-C0" << i;
            ifstream cFile( cFilename.str() );
            int lowestCommonAncestor = 0;

            for ( auto & item : shortInputLines[i] )
            {
                if ( lastItem < item )
                {
                    ostringstream outStream;
                    outStream << "Pile " << i << ": " << item << " => ";

                    cFile.seekg( item.posB * sizeof( MetagFileNumRefType ) );
                    vector<int> commonTaxa( taxLevelSize, 0 );
                    for ( LetterNumber j = 0; j < item.numB; ++j )
                    {
                        MetagFileNumRefType fileNum;
                        cFile.read( reinterpret_cast<char *>( &fileNum ), sizeof( MetagFileNumRefType ) );
                        outStream << fileNum << " ";

                        restrictCommonTaxaToIncludeFileNum( commonTaxa, fileNum );
                    }

                    //                    cout << outStream.str() << " => ";
                    lowestCommonAncestor = 0;
                    for ( auto & item : commonTaxa )
                    {
                        //                        cout << item << " ";
                        if ( item > 0 )
                            lowestCommonAncestor = item;
                    }

                    lastItem = item;
                }
                else
                {
                    //                    cout << "same";
                }

                //                cout << " => " << lowestCommonAncestor << endl;
                ++( lcaCounts[lowestCommonAncestor] );
            }
        }
    }

    cout << "Counts: ";
    for ( auto & item : lcaCounts )
        cout << item.first << ":" << item.second << " ";
    cout << endl;

    return 0;
}
