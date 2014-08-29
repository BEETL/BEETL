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

#ifdef HAVE_SEQAN

#include "../shared/Alphabet.hh"

#include <iostream>
#include <seqan/index.h>

using namespace std;
using namespace seqan;


void getFileName( const string &stem, const char code, const int pile,
                  string &fileName )
{
    fileName = stem;
    fileName += '-';
    fileName += code;
    fileName += '0';
    assert( pile <= 9 );
    fileName += ( char )( 48 + pile );
    // cerr << "Made file name " << fileName << endl;
}


int main ( int numArgs, const char *args[] )
{
    if ( numArgs != 3 )
    {
        cerr << "Usage: " << args[0] <<  "outputprefix  fileToConvert" << endl;
        exit( EXIT_FAILURE );
    }
    cerr << "metabeetl-db-makeBWTSkew inputFile " << args[1] << endl;
    std::fstream fstrm;
    fstrm.open( args[1], ::std::ios_base::in | ::std::ios_base::binary );
    String<char> fasta_tag;
    String<char> fasta_seq;
    //Read the meta-information.
    readMeta( fstrm, fasta_tag, Fasta() );
    std::cout << "Tag: " << fasta_tag << "\n"; //prints "a test file"
    //Read the sequence.
    read( fstrm, fasta_seq, Fasta() );
    fstrm.close();
    for ( uint i( 0 ); i < length( fasta_seq ); i++ )
    {
        fasta_seq[i] = toupper( fasta_seq[i] );
        if ( strchr( "ACGNT", fasta_seq[i] ) == NULL )
        {
            cerr << "Found invalid character " << fasta_seq[i]
                 << " at position " << i << " in file "
                 << args[1] << endl;
            fasta_seq[i] = 'N';
        }
    }
    appendValue( fasta_seq, '$' );
    // std::cout <<"Seq: "<< fasta_seq << "\n"; //prints the sequence
    // ModifiedString< ModifiedString< String<char>, ModView< FunctorComplement<char> > >, ModReverse > myMod(fasta_seq);
    String<unsigned> sa;
    ///Build a suffix array using the Skew7 algorithm.
    resize( sa, length( fasta_seq ) );
    createSuffixArray( sa, fasta_seq, Skew7() );
    char lastChar( notInAlphabet );
    int fileNum;
    string fileName;
    FILE *pFile( NULL ), *pArrayFile( NULL );
    cout << length( fasta_seq ) << endl;
    for ( uint i( 0 ); i < length( fasta_seq ); i++ )
    {
        if ( sa[i] > length( fasta_seq ) - 1 )
        {
            cerr << "sa bigger " << endl;
        }
        char thisChar = fasta_seq[sa[i]];
        if ( thisChar != lastChar )
        {
            fileNum = whichPile[( int )thisChar];
            assert( fileNum != nv );
            string fileNameStem = "bwt_" + ( string )args[2];
            getFileName( fileNameStem, 'B', fileNum, fileName );
            cerr << "Opening new file " << fileName << endl;
            if ( pFile != NULL )
            {
                fclose( pFile );
                pFile = NULL;
            }
            pFile = fopen( fileName.c_str(), "w" );
            assert( pFile != NULL );
            getFileName( fileNameStem, 'A', fileNum, fileName );
            cerr << "Opening new file " << fileName << endl;
            if ( pArrayFile != NULL )
            {
                fclose( pArrayFile );
                pArrayFile = NULL;
            }
            pArrayFile = fopen( fileName.c_str(), "w" );
            assert( pArrayFile != NULL );
            lastChar = thisChar;
        }
        thisChar = ( ( sa[i] == 0 ) ? '$' : fasta_seq[sa[i] - 1] );
        fputc( thisChar, pFile );
        fwrite( &sa[i], sizeof( unsigned ), 1, pArrayFile );
        // cout << fasta_seq[sa[i]] << endl;
    } // ~for i
    cout << endl;
    if ( pFile != NULL ) fclose ( pFile );
    if ( pArrayFile != NULL ) fclose ( pArrayFile );
    // for (int i(0);i<10;i++)
    // cout << fasta_seq[sa[i]] << end
    return 0;
} // ~main


#else //ifdef HAVE_SEQAN

#warning Compiled without Seqan library. Missing metagenomics tool.
#include <iostream>

int main ( int numArgs, const char *args[] )
{
    std::cerr << "This tool is unavailable as BEETL was compiled without the Seqan library" << std::endl;
}

#endif //ifdef HAVE_SEQAN
