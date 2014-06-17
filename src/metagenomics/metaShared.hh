#ifndef META_SHARED_HH
#define META_SHARED_HH

#include "../shared/Tools.hh"

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using namespace std;


struct BWTInformation
{
    unsigned short pileNum_;
    double readCount_;
    //  unsigned double readAverage_;
    unsigned short wordLength_;
    unsigned int taxId_;
    unsigned short taxLevel_;
    vector< unsigned int > charBCount_;
    vector< unsigned int > charACount_;
    vector<unsigned short> fileNumbers_;
    vector<int> readIds;
};

struct TaxInformation
{
    unsigned short taxLevel_;
    vector<int> files_;
    double *wordCountPerSize_;
    double *wordCountPerSizeOfChildren_;
    string name_;
    vector <uint64_t> bwtPositions_;
    vector <double>  wordCounts_;
    vector <unsigned short> wordLengths_;
    vector <unsigned short> pileNumbers_;
    int parentId_;
    double normalisedCount_;
    uint64_t seqLengthSum_;

    TaxInformation()
        : taxLevel_( 0 )
        , files_( 0 )
        , wordCountPerSize_( NULL )
        , wordCountPerSizeOfChildren_( NULL )
        , parentId_( 0 )
        , normalisedCount_( 0 )
        , seqLengthSum_( 0 )
    {}
};

struct FileInformation
{
    vector<unsigned> suffixPos_;
    vector< int > suffixCounts_;
    vector<unsigned short > suffixLengths_;
    vector<uint64_t > bwtPositions_;
    int sequenceLength_;
    vector<char> suffixChar_;
};



typedef map<int, TaxInformation> TAXMAP;
typedef map<int, FileInformation> FILEMAP;

typedef map <uint64_t, BWTInformation > BWTMAP;

typedef map <unsigned int, string> READMAP;



map<int, TaxInformation> loadTaxInformationForDatabase( string taxInfoFilename, int cycleSize, string ncbiNamesDMP )
{
    ifstream dbTax( taxInfoFilename );
    map<int, TaxInformation> taxInformation;
    string line;
    cerr << "loading NCBI taxonomy " << endl;
    ifstream ncbiNames( ncbiNamesDMP.c_str(), ios::in );
    map<int, string> idToScientificName;
    while ( ncbiNames.good() )
    {
        getline( ncbiNames, line );
        vector<string> lineVector = splitString( line , "\t|\t" );
        if ( lineVector.size() < 2 ) continue;
        string name = lineVector[1];
        //    cout << "ncbiName >"<< name <<"<" <<endl;
        int id = atoi( lineVector[0].c_str() );
        if ( line.find( "scientific" ) != string::npos )
            idToScientificName[id] = name;
    }

    cerr << "Got NCBI Taxonomy " << idToScientificName.size() << endl;

    while ( dbTax.good() )
    {
        getline( dbTax, line );
        vector<string> splitTax = splitString( line, " " );
        int fileNum = atoi( splitTax[0].c_str() );
        for ( unsigned int i( 1 ) ; i < splitTax.size(); ++i )
        {
            int taxId = atoi( splitTax[i].c_str() );
            //cout <<" taxId "  <<taxId <<endl;
            bool taxFound( false );
            for ( TAXMAP::iterator it( taxInformation.begin() ) ; it != taxInformation.end(); ++it )
            {
                if ( it->first == taxId )
                {
                    //	  cout << (i-1) << " "<< countForTaxLevel[i-1] <<endl;
                    taxFound = true;
                    break;
                }
            }
            if ( taxFound )
                taxInformation[taxId].files_.push_back( fileNum );
            else
            {
                TaxInformation tax;
                if ( i > 1 )
                {
                    tax.parentId_ = atoi( splitTax[i - 1].c_str() );
                    if ( tax.parentId_ == 0 )
                        tax.parentId_ =  atoi( splitTax[i - 2].c_str() );
                    //	  cerr << taxId << " parent " << tax.parentId_ <<endl;
                }
                else //set the root as a parent Id
                    tax.parentId_ = 1;

                tax.files_.push_back( fileNum );
                tax.taxLevel_ = i - 1;
                tax.name_ = idToScientificName[taxId];
                tax.seqLengthSum_ = 0;
                //countForTaxLevel[i - 1] += 1;
                //	cout <<countForTaxLevel[i-1] << endl;
                //initialise all word counts with zero
                //cycle 0 starts with two chars
                tax.normalisedCount_ = 0;
                tax.wordCountPerSize_ = new double[cycleSize + 2];
                tax.wordCountPerSizeOfChildren_ = new double[cycleSize + 2];

                for ( int i( 0 ) ; i < cycleSize + 1; ++i )
                    tax.wordCountPerSize_[i] = 0;
                taxInformation[taxId] = tax;
            }
        }
    }
    //set the root
    taxInformation[1].taxLevel_ = -1;
    taxInformation[1].wordCountPerSize_ = new double[cycleSize + 2];
    taxInformation[1].wordCountPerSizeOfChildren_ = new double[cycleSize + 2];
    taxInformation[1].name_ = "root";
    for ( int i( 0 ) ; i < cycleSize + 1; ++i )
        taxInformation[1].wordCountPerSize_[i] = 0;
    dbTax.close();
    return taxInformation;
}


int getTaxonomicLevel( string s )
{
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
    else if ( s.compare( "subspecies" ) == 0 )
        level = 7;
    else if ( s.compare( "strain" ) == 0 )
        level = 7;
    else if ( s.compare( "subsubspecies" ) == 0 )
        level = 8;
    else
        level = 100;
    return level;
}

#endif // META_SHARED_HH
