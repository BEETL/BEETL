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
 **
 ** Collection of basic tools and defines
 **/

#include "Tools.hh"

#include <cstdlib>

using namespace std;


double Tools::startTime;

void Tools::StartTimer()
{
    static const double in_seconds = 1.0 / static_cast<double>( CLOCKS_PER_SEC );
    startTime = clock() * in_seconds;
}

double Tools::GetTime()
{
    static const double in_seconds = 1.0 / static_cast<double>( CLOCKS_PER_SEC );
    double curTime = clock() * in_seconds;
    return curTime - startTime;
}



uchar *Tools::GetRandomString( unsigned min, unsigned max, unsigned &alphabetSize )
{
    unsigned len = rand() % ( max - min ) + min;
    alphabetSize = rand() % 26 + 1;
    uchar *temp = new uchar[len + 2];
    for ( unsigned i = 0; i < len; i++ )
        temp[i] = 97 + rand() % alphabetSize;
    temp[len] = 0u ;
    temp[len + 1] = '\0';
    return temp;
}







unsigned Tools::FloorLog2( ulong i )
{
    ulong b = 0;
    if ( i == 0 )
        return 0;
    while ( i )
    {
        b++;
        i >>= 1;
    }
    return b - 1;
}

//Creating table to find logn in small time
unsigned *Tools::MakeTable()
{
    unsigned *table = new unsigned[512];
    for ( unsigned i = 0; i < 9; i++ )
    {
        if ( i == 0 )
            table[i] = 0;
        if ( i >= 1 && i < ( 1 << 1 ) )
            table[i] = 1;
        if ( i >= ( 1 << 1 ) && i < ( 1 << 2 ) )
            table[i] = 2;
        if ( i >= ( 1 << 2 ) && i < ( 1 << 3 ) )
            table[i] = 3;
        if ( i >= ( 1 << 3 ) && i < ( 1 << 4 ) )
            table[i] = 4;
        if ( i >= ( 1 << 4 ) && i < ( 1 << 5 ) )
            table[i] = 5;
        if ( i >= ( 1 << 5 ) && i < ( 1 << 6 ) )
            table[i] = 6;
        if ( i >= ( 1 << 6 ) && i < ( 1 << 7 ) )
            table[i] = 6;
        if ( i >= ( 1 << 7 ) && i < ( 1 << 8 ) )
            table[i] = 7;
        if ( i >= ( 1 << 8 ) && i < ( 1 << 9 ) )
            table[i] = 8;
    }
    return table;
}

unsigned Tools::FastFloorLog2( unsigned i )
{
    unsigned *table = MakeTable();
    unsigned u;
    if ( i >> 24 )    u = 22 + table[ i >> 24] ;
    if ( i >> 16 )    u = 14 + table[ i >> 16] ;
    if ( i >> 8 )     u = 6 + table[ i >> 8] ;
    u =  table[i] - 1;
    delete [] table;
    return u;
}

unsigned Tools::CeilLog2( ulong i )
{
    unsigned j = FloorLog2( i );
    if ( ( ulong )( 1lu << j ) != i )
        return j + 1;

    return j;
}

uchar *Tools::GetFileContents( char *filename, ulong maxSize )
{
    std::ifstream::pos_type posSize;
    std::ifstream file ( ( char * )filename, std::ios::in | std::ios::binary | std::ios::ate );
    if ( file.is_open() )
    {
        posSize = file.tellg();
        ulong size = posSize;
        if ( maxSize != 0 && size > maxSize )
            size = maxSize;
        char *memblock = new char [size + 1];
        file.seekg ( 0, std::ios::beg );
        file.read ( memblock, size );
        memblock[size] = '\0';
        file.close();
        return ( uchar * )memblock;
    }
    else
        return 0;
}

void getFileName( const string &stem, const char code, const int pile,
                  string &fileName )
{
    fileName = stem;
    fileName += '-';
    fileName += code;
    fileName += '0';
    assert( pile <= 9 );
    fileName += ( char )( 48 + pile );
    //  cerr << "Made file name " << fileName << endl;
}

bool isValidFastaFile( const char *fileName )
{
    FILE *file;
    char probe_char;

    file = fopen( fileName, "r" );

    if ( file != NULL )
    {
        cerr << "Checking fasta file " << fileName << endl;
    }
    else
    {
        cerr << "Could not read fasta file " << fileName << endl;
        return 0;
    }

    // grep first char, should be a >, if not print error

    probe_char = getc( file );
    if ( probe_char != EOF && probe_char != '>' )
    {
        cerr << "Validation of fasta file " << fileName
             << " failed. Maybe wrong format?" << endl;
        return 0;
    }
    ungetc( probe_char, file );

    fclose( file );
    return 1;
}

bool isValidReadFile( const char *fileName )
{
    FILE *file;
    char probe_char;

    file = fopen( fileName, "r" );

    if ( file != NULL )
    {
        cerr << "Checking read file " << fileName << endl;
    }
    else
    {
        cerr << "Could not read file " << fileName << endl;
        return 0;
    }

    // grep first char, should be a ACGTN, if not print error

    probe_char = getc( file );
    // only basic check, no EOF and no fastq or fasta header
    if ( probe_char != EOF && ( probe_char == '@' || probe_char == '>' ) )
    {
        cerr << "Validation of read file " << fileName
             << " failed. Maybe wrong format?" << endl;
        return 0;
    }
    ungetc( probe_char, file );

    fclose( file );
    return 1;
}

void readWriteCheck( const char *fileName, const bool readWrite )
{
    FILE *file;
    file = fopen( fileName, readWrite ? "w" : "r" );
    string mode = readWrite ? "writing." : "reading.";
    if ( file == NULL )
    {
        cerr << "Could not open file " << fileName << " for "
             << mode << " Aborting." << endl;
        exit( EXIT_FAILURE );
    }
    fclose( file );
}

void checkIfEqual( const int arg1, const int arg2 )
{

    if ( arg1 != arg2 )
    {
        cerr << "Validation failed: " <<  arg1 << " != "
             << arg2 << ". Should be (==). Aborting." << endl;
        exit( EXIT_FAILURE );
    }
}

void checkIfNotEqual( const int arg1, const int arg2 )
{

    if ( arg1 == arg2 )
    {
        cerr << "Validation failed: " <<  arg1 << " == "
             << arg2 << ". Should be (!=). Aborting." << endl;
        exit( EXIT_FAILURE );
    }
}

bool hasPrefix( const string &fullString, const string &prefix )
{
    return fullString.compare( 0, prefix.size(), prefix ) == 0;
}

bool hasSuffix( const string &fullString, const string &suffix )
{
    return fullString.size() >= suffix.size()
           && fullString.compare( fullString.size() - suffix.size(), suffix.size(), suffix ) == 0;
}
