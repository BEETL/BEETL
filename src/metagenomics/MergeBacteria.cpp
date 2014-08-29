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
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <queue>
#include <vector>

#include "Types.hh"

using namespace std;


//#define DEBUG;

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


long getFileSize( FILE *file )
{
    long lCurPos, lEndPos;
    lCurPos = ftell( file );
    fseek( file, 0, 2 );
    lEndPos = ftell( file );
    fseek( file, lCurPos, 0 );
    return lEndPos;
}

struct SuffixInfo
{
    // FILE* pFile_;
    MetagFileNumRefType fileNum_;
    unsigned char bwtSymbol_;
    const char *pSeq_;
    unsigned loc_;
};

bool operator<( const SuffixInfo &lhs, const SuffixInfo &rhs )
{
    // assert(1==0);
    int compare( strcmp( lhs.pSeq_ + lhs.loc_, rhs.pSeq_ + rhs.loc_ ) );

#ifdef DEBUG

    cout << "compare l " << lhs.fileNum_ << " " << lhs.loc_ <<  "  " << strlen( lhs.pSeq_ ) << endl;
    cout << "compare r " << rhs.fileNum_ << " " << rhs.loc_ << "  " << strlen( rhs.pSeq_ ) << endl;
    cout << compare << endl;
    int i( 0 );
    while ( lhs.pSeq_[lhs.loc_ + i] == rhs.pSeq_[rhs.loc_ + i] ) i++;
    cout << i << " " << lhs.pSeq_[lhs.loc_ + i] << " " << rhs.pSeq_[rhs.loc_ + i] << endl;
#endif

    // if (compare==0) cout << "Hello" << endl;
    //       cout << lhs.fileNum_ << " " << lhs.loc_ << " " << rhs.fileNum_ << " " << rhs.loc_ << " " << lhs.pSeq_+lhs.loc_ << " " << (long)rhs.pSeq_ << endl;
    return ( ( compare == 0 ) ? ( lhs.fileNum_ > rhs.fileNum_ ) : ( compare > 0 ) );
    // return ((compare==0)?(lhs.fileNum_>rhs.fileNum_):(compare>0));

    // return (strcmp(lhs.pSeq_+lhs.loc_,rhs.pSeq_+rhs.loc_)>0);
}

long getFileSize( const string &name )
{
    ifstream file( name );

    // Find file size
    file.seekg( 0, ios::end );
    size_t fileSize = file.tellg();

    return fileSize;
}

void readFastaFile( const char *name, vector<char> &data )
{
    ifstream fastaFile( name, ios::in );
    size_t fileSize = getFileSize( name );
    data.reserve( data.size() + fileSize + 1 );

    string line;
    while ( getline( fastaFile, line ) )
    {
        if ( line[0] != '>' )
        {
            for ( unsigned int i( 0 ) ; i < line.length() ; i++ )
            {
                char c = toupper( line[i] );
                if ( strchr( "ACGNT", c ) == NULL )
                {
                    cout << "found wrong char " << c << " in " << name << " => replaced with 'N'." << endl;
                    c = 'N';
                }
                data.push_back( c );
            }
        }
    }
    fastaFile.close();
    data.push_back( '\0' );
} // ~readFastaFile

void readFileList( const char *fileName, std::vector<string> &files )
{
    ifstream fileList( fileName, ios::in );
    if ( ! fileList.is_open() )
    {
        cerr << "Cannot open " << fileName << " for reading." << endl;
        exit( EXIT_FAILURE );
    }
    string line;
    while ( getline( fileList, line ) )
    {
        files.push_back( line );
    }
    fileList.close();
}


int main ( int numArgs, const char *args[] )
{
    if ( numArgs != 4 )
    {
        cerr << "Usage: " << args[0] << " pileNum outputPrefix fileList" << endl;
        exit( EXIT_FAILURE );
    }

    int pileNum( atoi( args[1] ) );
    cerr << "Will merge BWTs from pile " << pileNum << endl;
    //assert((pileNum>=0)&&(pileNum<alphabetSize));

    cerr << "Output file will be named " << args[2] << endl;
    cerr << "Retrieving list of BWT files from " << args[3] << endl;

    vector<vector<char> > chromSeqs;
    vector<queue<unsigned> > suffixArrays;
    //  vector<vector<unsigned>::iterator> suffixArraysEnd;
    vector<queue<char> > bwts;
    //  vector<string::iterator> bwtsEnd;

    //  vector<vector<unsigned> > suffixstore;
    //vector<string> bwtstore;

    FILE *pFileChrom;

    FILE *pFileArray;
    FILE *pFileBWT;

    string fileNameStem( args[2] ), fileName;

    getFileName( fileNameStem, 'A', pileNum, fileName );
    cerr << "Will send suffix array entries to " << fileName << endl;
    pFileArray = fopen( fileName.c_str(), "w" );
    assert( pFileArray != NULL );

    getFileName( fileNameStem, 'B', pileNum, fileName );
    cerr << "Will send BWT entries to " << fileName << endl;
    pFileBWT = fopen( fileName.c_str(), "w" );
    assert( pFileBWT != NULL );

    getFileName( fileNameStem, 'C', pileNum, fileName );
    cerr << "Will send genome information to " << fileName << endl;
    pFileChrom = fopen( fileName.c_str(), "w" );
    assert( pFileChrom != NULL );


    SuffixInfo thisSuffix;
    priority_queue<SuffixInfo, vector<SuffixInfo> > pq;
    ofstream fileCounter( "filecounter.csv", ios::out );
    int fileN( 0 );
    vector<string> files;
    readFileList( args[3], files );
    cerr << "Collected " << files.size() << " files" << endl;

    for ( uint i( 0 ); i < files.size(); ++i )
    {
        //      cout << i <<endl;
        //read the fasta file in
        chromSeqs.push_back( vector<char>() );
        readFastaFile( files[i].c_str(), chromSeqs.back() );
        //readFastaFile( args[i], chromSeqs.back() );
        //readFastaFile( args[i], chromSeqs.back() );

        //cerr << "Sequence is of length " << chromSeqs.back().size() << endl;
        if ( i % 20 == 0 )
        {
            cerr << "Sequence " << i <<  " has length " << chromSeqs.back().size() << endl;
        }
        //fileNameStem = "bwt_" + ( string )args[i];
        fileNameStem = "bwt_" + files[i];
        getFileName( fileNameStem, 'A', pileNum, fileName );
        //      cerr << "Will open suffix array file " << fileName << endl;
        //read in the whole suffix array
        queue<unsigned> thisSuffArr;
        queue<char> bwtQue;
        suffixArrays.push_back( thisSuffArr );
        bwts.push_back( bwtQue );
        FILE *arrayInputStream = fopen( fileName.c_str(), "r" );
        if ( arrayInputStream != NULL )
        {
            long fileSize = getFileSize( arrayInputStream );
            //allocate the space
            unsigned *arrayBuf = new unsigned[fileSize / 4];
            assert( ( int )fread( arrayBuf, 4, fileSize / 4, arrayInputStream ) == ( fileSize / 4 ) );
            fileCounter << fileN << "," << files[i] << endl;
            //fileCounter << fileN << "," << args[i] << endl;
            fileN++;
            //      cout << "read suffix\n";
            for ( int k( 0 ) ; k < ( fileSize / 4 ); k++ )
            {
                thisSuffArr.push( arrayBuf[k] );
            }
            suffixArrays.back() = thisSuffArr;
            fclose( arrayInputStream );
        }
        if ( !suffixArrays.back().empty() )
        {
            //read the bwt
            getFileName( fileNameStem, 'B', pileNum, fileName );

            if ( i % 20 == 0 )
            {
                cerr << "Will open BWT file " << fileName << endl;
            }

            ifstream bwtIn( fileName );
            if ( bwtIn.good() )
            {
                long fileSize = getFileSize( fileName );
                //allocate the space
                vector<char> arrayBuf( fileSize );
                bwtIn.read( arrayBuf.data(), fileSize );

                for ( int k( 0 ) ; k < fileSize; ++k )
                {
                    bwtQue.push( arrayBuf[k] );
                }
                bwts.back() = bwtQue;
            }
        }
        //      cerr <"got bwt" <<endl;
        else
        {
            cerr << "Warning: could not open " << fileName << endl;
        }
    }
    //chromSeqs = whole sequence out of fasta file
    // for each bacterium initialise one suffix
    cerr << "Got " << chromSeqs.size() << " sequences " << endl;
    cerr << "Got " << suffixArrays.size() << " suffixe" << endl;
    fileCounter.close();
    for ( unsigned int i( 0 ); i < chromSeqs.size(); i++ )
    {
        if ( !suffixArrays[i].empty() )
        {
            //pointer to first char in vector<char> of the genome sequence
            thisSuffix.pSeq_ = &( chromSeqs[i][0] );
            thisSuffix.fileNum_ = i;
            thisSuffix.loc_ = suffixArrays[i].front();
            suffixArrays[i].pop();
            thisSuffix.bwtSymbol_ = bwts[i].front();
            bwts[i].pop();
            //      bwtIndices[i] += 1;
            pq.push( thisSuffix );
        }
    }
    // ~for

    int j( 0 );

#ifdef DEBUG
    const char *prevSuff( NULL );
    const char *thisSuff( NULL );
    unsigned char prevNum;
#endif

    while ( !pq.empty() )
    {
        thisSuffix = pq.top();

#ifdef DEBUG
        for ( int jj( 0 ); jj < 50 ; jj++ )
            cout << *( thisSuffix.pSeq_ + thisSuffix.loc_ + jj );
        cout << endl;

        thisSuff = thisSuffix.pSeq_ + thisSuffix.loc_;
        if ( prevSuff != NULL )
        {
            assert( ( strcmp( prevSuff, thisSuff ) < 1 ) ||
                    ( ( strcmp( prevSuff, thisSuff ) == 0 ) && ( prevNum < thisSuffix.fileNum_ ) ) );
        }
        prevSuff = thisSuff;
        prevNum = thisSuffix.fileNum_;
#endif

        pq.pop();
        assert( fwrite( &thisSuffix.loc_, sizeof( unsigned ), 1, pFileArray ) == 1 );
        assert( fwrite( &thisSuffix.bwtSymbol_, sizeof( unsigned char ), 1, pFileBWT ) == 1 );
        assert( fwrite( &thisSuffix.fileNum_, sizeof( MetagFileNumRefType ), 1, pFileChrom ) == 1 );

        if ( ( j % 1000000 ) == 0 ) cout << "." << endl;
        j++;
#ifdef DEBUG
        //   for (int j(0);j<50;j++)
        //    	cout << *(thisSuffix.pSeq_+thisSuffix.loc_+j);
        // cout << endl;
#endif
        //read 4 bytes

        if ( !suffixArrays[thisSuffix.fileNum_].empty() )
        {
            thisSuffix.loc_ = suffixArrays[thisSuffix.fileNum_].front();
            suffixArrays[thisSuffix.fileNum_].pop();
            //	cout << "suffix loc assigned " << suffixArrays[thisSuffix.fileNum_][suffixIndices[thisSuffix.fileNum_]];
            //	suffixIndices[thisSuffix.fileNum_] +=1;
            thisSuffix.bwtSymbol_ = bwts[thisSuffix.fileNum_].front();
            bwts[thisSuffix.fileNum_].pop();
            //	cout << "bwt assigned " << bwts[thisSuffix.fileNum_][bwtIndices[thisSuffix.fileNum_]];
            //	bwtIndices[thisSuffix.fileNum_] +=1;
            //open bwt

            //	  bwtInputStream = fopen(filesBWT[thisSuffix.fileNum_].c_str(),"r");
            //set right position
            //  fseek(bwtInputStream, bwtStreamPositions[thisSuffix.fileNum_], SEEK_SET);
            //read one byet
            //  assert(fread(&thisSuffix.bwtSymbol_,sizeof(unsigned char),1,
            //	       bwtInputStream)==1);
            //remember new position
            //  bwtStreamPositions[thisSuffix.fileNum_] += sizeof(unsigned char);
            //close bwt
            //	  fclose(bwtInputStream);
            pq.push( thisSuffix );
        }
    } // ~while


    /*  for (unsigned int i(0);i<files.size();i++)
      {
        if (files[i]!=NULL)
    {
      fclose(files[i]);
    }

        if (filesBWT[i]!=NULL)
    {
      fclose(filesBWT[i]);
    }
      }
    */
    fclose( pFileArray );
    fclose( pFileBWT );
    fclose( pFileChrom );

} // ~main
