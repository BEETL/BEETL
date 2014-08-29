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

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;


void createSingeSeqFilesPlusReference( const vector<string> &fileNames, const string &savingDir );
string reverseComplementSequence( const string &s );

int main ( int argc, char *argv[] )
{
    if ( argc < 3 )
    {

        cerr << "This script takes an array of fasta files and a saving directory. This is the first step in generating a metagenome database for count words." << endl
             << "For each fasta file two single sequence files will be saved. One with the singe sequence from the file and one with the reverse complement. " << endl
             << "Plasmids and transposons will be omitted through a header search. " << endl
             << "Additional a file with the name headerFile.cvs will be created. There all new file numbers will be connected with their fasta header" << endl << endl ;


        cerr << "Ussage: " << endl
             << "./GenomeToSingleseq -f <input_files> -s output_dir " << endl
             << "-f list of fasta files " << endl
             << "-s output directory" << endl;
    }
    else
    {
        vector<string> inputFiles;
        string savingDir;
        for ( int i( 0 ) ; i < argc; i++ )
        {
            if ( argv[i][1] == 'f' )
            {
                for ( int j = i + 1; j < argc; ++j )
                {
                    if ( argv[j][0] != '-' )
                        inputFiles.push_back( argv[j] );
                    else
                        break;
                }
            }
            else if ( argv[i][1] == 's' )
            {
                savingDir = argv[i + 1];
            }
        }
        createSingeSeqFilesPlusReference( inputFiles, savingDir );
    }
}

void createSingeSeqFilesPlusReference( const vector<string> &fileNames, const string &savingDir )
{
    mkdir( savingDir.c_str(), 0750 );
    ofstream headerOut;
    string headerFile = savingDir + "/headerFile.csv";
    headerOut.open( headerFile.c_str(), ios::out );
    int genomeCount( 1 );
    for ( uint i( 0 ) ; i < fileNames.size(); i++ )
    {
        string fileName = fileNames[i];
        cout << "inputFile " << fileName << endl;

        ifstream inputFas;
        inputFas.open( fileName.c_str(), ios::in );
        string line;
        string sequence;
        //bool first(true);

        string name;

        //  cout << "inputFile "<< fileName << endl;
        stringstream gof;
        gof << savingDir << "/G_" << genomeCount;
        string  genomeOutFile = gof.str();
        stringstream rgof;
        rgof << savingDir  << "/G_" << genomeCount <<  "_rev";
        string revGenomeOutFile = rgof.str();
        cout << "output " << genomeOutFile << endl;
        while ( inputFas.good() )
        {
            getline( inputFas, line );

            if ( line[0] != '>' )
            {
                sequence += line;
            }
            else
            {
                name = line.substr( 1, line.length() - 1 );
                cout << "header " << name << endl;
            }
        }
        if ( name.find( "plasmid" ) == string::npos && name.find( "transposon" ) == string::npos )
        {
            ofstream genomeOut;
            genomeOut.open( genomeOutFile.c_str(), ios::out );
            genomeOut << ">" << name << endl;
            genomeOut << sequence << endl;

            ofstream revGenomeOut;
            revGenomeOut.open( revGenomeOutFile.c_str(), ios::out );
            revGenomeOut << ">Reverse " <<  name << endl;
            revGenomeOut << reverseComplementSequence( sequence ) << endl;
            genomeOut.close();
            revGenomeOut.close();
            headerOut << genomeCount << "," << name << endl;
            genomeCount++;
        }
        else
        {
            cerr << "did not change " << name << endl;
        }
        //    cout << "wrote files " << endl;
    }
    cout << "done\n";
}




string reverseComplementSequence( const string &s )
{
    string revChar = s;
    //  cout<< "s "  << s << endl;
    int j( 0 );
    for ( int i = s.length() - 1; i >= 0 ; i-- )
    {
        //      cout << "i "<<  i <<  endl;
        if ( s[i] == 'A' || s[i] == 'a' )
        {
            revChar[j] = 'T';
        }
        else if ( s[i] == 'G' || s[i] == 'g' )
        {
            revChar[j] = 'C';
        }
        else if ( s[i] == 'C' || s[i] == 'c' )
        {
            revChar[j] = 'G';
        }
        else if ( s[i] == 'T' || s[i] == 't' )
        {
            revChar[j] = 'A';
        }
        else if ( s[i] == 'N' || s[i] == 'n' )
        {
            revChar[j] = 'N';
        }
        else
        {
            revChar[j] = s[i] ;
            //  cout << "found " <<s[i] << "in the genome \n";
        }
        j++;
    }
    return revChar;
}
