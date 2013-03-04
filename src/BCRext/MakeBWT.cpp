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


// MakeBWT.cpp
// Function:
// 1. Take a set of chromosomes
// 2. Make one big $-separated string
// 3. Append its reverse complement
// 4. Build the BWT of the whole lot

#include "Alphabet.hh"
#include "Timer.hh"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <sys/time.h>

using namespace std;


//#define DEBUG 1

typedef char TranslationTable[256];
const char nr( '?' );


static const TranslationTable reverseCharASCII =
{
    nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,
    nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,
    nr,nr,nr,nr,'$',nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,
    nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,
    nr, /* next is 'A' */
    'T',nr,'G',nr,nr,nr,'C',nr,nr,nr,nr,nr,nr,'N',nr, /* next is 'P' */
    nr,nr,nr,nr,'A',nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,
    nr, /* next is 'a' */
    't',nr,'g',nr,nr,nr,'c',nr,nr,nr,nr,nr,nr,'n',nr, /* next is 'p' */
    nr,nr,nr,nr,'a',nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,
    nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,
    nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,
    nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,
    nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,
    nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,
    nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,
    nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,
    nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr,nr
};


#ifdef XXX
char reverseComp( const char c )
{
    return ( c=='A' )
           ? 'T'
           : ( ( c=='T' )
               ? 'A'
               : ( ( c=='C' )
                   ? 'G'
                   : ( ( c=='G' ) ? 'C' : '!' ) ) );
} // ~reverseComp;
#endif


vector<char> charArray;

typedef unsigned long long ulonglong;

vector<const char *> pointerArray;



const ulonglong charArraySize=6400000000; // needs to be > 2 x genome size


class LessThan
{
public:
    LessThan() : pChar( &charArray[0] ) {}

    bool operator()( const char *pa, const char *pb )
    {

        //  const char *pa(&charArray[a]),*pb(&charArray[b]);
        //    const char *pa(pChar+a),*pb(pChar+b)

        while ( 1 )
        {
            if ( *pa!=*pb )
            {
                return ( *pa<*pb );
            } // ~if
            else
            {
                if ( *pa=='$' )
                    return ( pa<pb );
                else
                {
                    pa++;
                    pb++;
                }
            } // ~else
        } // ~while
    } // ~operator()
private:
    const char *const pChar;
}; // ~LessThan



void getFileName( const string &stem, const char code, const int pile,
                  string &fileName )
{
    fileName=stem;
    fileName+='-';
    fileName+=code;
    fileName+='0';
    assert( pile<=9 );
    fileName+=( char )( 48+pile );
    //  cerr << "Made file name " << fileName << endl;
}



int main( int numArgs, char *args [] )
{
    if ( numArgs<3 )
    {
        cerr << "syntax: " << args[0]
             << " outputPrefix listOfFiles" << endl;
        exit( EXIT_FAILURE );
    } // ~if

    cerr << "Resizing char array to size " << charArraySize << endl;
    assert ( sizeof( ulonglong )==8 );
    assert ( reverseCharASCII[( int )'$']=='$' );

    charArray.resize( charArraySize );

    cerr << "Resizing char array to size " << charArraySize << endl;

    Timer timer;

    cerr << "Output BWT files will be named " << args[1] << endl;


    char *pStartChar( &charArray[0] );
    char *pThisChar( &charArray[0] );
    char *pCharToReverse;
    FILE *pFile;
    FILE *pArrayFile;

    cerr << "Reading files, time now: " << timer.timeNow();
    cerr << "Reading files, usage: " << timer << endl;

    for ( int thisFile( 2 ); thisFile<numArgs; thisFile++ )
    {
        cerr << "Opening file " << args[thisFile] << " at position " << ( pThisChar-pStartChar ) << endl;
        *pThisChar++='$';
        pFile=fopen( args[thisFile],"r" );
        assert( pFile!=NULL );

        assert( fgets( pThisChar, 1024, pFile )!=NULL );
        cerr << "Read first line " << pThisChar;

        while ( fgets( pThisChar, 1024, pFile )!=NULL )
        {
            pThisChar+=strlen( pThisChar )-1;
        }

        fclose( pFile );
    } // ~for


    cerr << "Reversing sequences, time now: " << timer.timeNow();
    cerr << "Reversing sequences, usage: " << timer << endl;

    cerr << "Reversing at position " << ( pThisChar-pStartChar ) << endl;

    pCharToReverse=pThisChar;
    *pThisChar++='$';
    int nonStandardChars( 0 );
    while ( ( --pCharToReverse )!=pStartChar )
    {
        *pThisChar=reverseCharASCII[( int )*pCharToReverse];

        if ( *pThisChar==nr )
        {
            *pCharToReverse='N';
            *pThisChar='N';
            ++nonStandardChars;
        } // ~if
        pThisChar++;
    } // ~while

    *pThisChar++='$';
    *pThisChar++='\0';

    cerr << "Read " << nonStandardChars
         << " non-ACGTN chars, converted those to N" << endl;
    //  cerr << pStartChar << endl;

    cerr << "Creating pointer array" << endl;

    charArray.resize( pThisChar-pStartChar );
    assert( ( pThisChar-1 )==&charArray.back() );

    string fileName;
    pFile=NULL;


    for ( int j( 0 ); j<alphabetSize; j++ )
    {
        cerr << "Starting pass " << j << ", time now: " << timer.timeNow();
        cerr << "Starting pass " << j << ", usage: " << timer << endl;

        getFileName( args[1],'B',j,fileName );
        cerr << "Opening new file " << fileName << endl;
        //    if (pFile!=NULL) fclose(pFile);
        pFile=fopen( fileName.c_str(),"w" );
        assert( pFile!=NULL );

        getFileName( args[1],'A',j,fileName );
        cerr << "Opening new file " << fileName << endl;
        //    if (pFile!=NULL) fclose(pFile);
        pArrayFile=fopen( fileName.c_str(),"w" );
        assert( pArrayFile!=NULL );

        pointerArray.clear();

        for ( ulonglong i( 1 ); i<charArray.size(); i++ )
        {
            //      if ((i%1000000)==0) cout << i << endl;
            //      cerr << i << charArray[i] << endl;
            if ( whichPile[( int )charArray[i]]==j )
                pointerArray.push_back( &charArray[i] );
        }

        //    pointerArray.resize(charArray.size()-2);
        //  for (int i(0);i<pointerArray.size();i++)
        //  pointerArray[i]=i+1;

        cerr << "Found " << pointerArray.size() << " suffixes beginning with "
             << alphabet[j] << endl;

        cerr << "Sorting suffixes" << endl;

        sort( pointerArray.begin(),pointerArray.end(),LessThan() );

        cerr << "Printing suffixes" << endl;

        //    char lastChar(nr);
        //  string fileName;
        //  int fileNum(0);

        //    pFile=NULL;
        ulonglong thisDiff;

        for ( ulonglong i( 0 ); i<pointerArray.size(); i++ )
        {
            fputc( *( pointerArray[i]-1 ),pFile );
            thisDiff=( ulonglong )( pointerArray[i]-pStartChar );
            fwrite( pointerArray[i],sizeof( ulonglong ),1,pArrayFile );
#ifdef DEBUG
            cerr << i << " " << pointerArray[i] << " " << charArray[pointerArray[i]-1]
                 << " " << &charArray[pointerArray[i]] << endl;
#endif
        } // ~for i

        fclose ( pFile );
        fclose ( pArrayFile );

    } // ~for j

    cerr << "Exiting, time now: " << timer.timeNow();
    cerr << "Exiting, usage: " << timer << endl;


} // ~main
