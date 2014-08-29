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
 **
 ** Collection of basic tools and defines
 **/

#include "Tools.hh"

#include "Alphabet.hh"
#include "libzoo/util/Logger.hh"

#include <cstdlib>
#include <dirent.h>
#include <sstream>
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>

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







uint8_t Tools::FloorLog2( uint64_t i )
{
    uint8_t b = 0;
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
uint8_t *Tools::MakeTable()
{
    uint8_t *table = new uint8_t[256];
    for ( unsigned i = 0; i < 256; i++ )
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

uint8_t Tools::FastFloorLog2( uint32_t i )
{
    static uint8_t *table = MakeTable();
    uint8_t u;
    if ( i >> 24 )       u = 22 + table[ i >> 24];
    else if ( i >> 16 )  u = 14 + table[ i >> 16];
    else if ( i >> 8 )   u =  6 + table[ i >>  8];
    else                 u = table[i] - 1;
    return u;
}

uint8_t Tools::CeilLog2( uint64_t i )
{
    uint8_t j = FloorLog2( i );
    if ( ( uint64_t )( 1lu << j ) != i )
        return j + 1;

    return j;
}

uchar *Tools::GetFileContents( char *filename, size_t maxSize )
{
    std::ifstream::pos_type posSize;
    std::ifstream file ( ( char * )filename, std::ios::in | std::ios::binary | std::ios::ate );
    if ( file.is_open() )
    {
        posSize = file.tellg();
        size_t size = posSize;
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

/*
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
*/

bool isValidFastaFile( const char *fileName )
{
    FILE *file;
    int probe_char;

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
        fclose( file );
        return 0;
    }
    ungetc( probe_char, file );

    fclose( file );
    return 1;
}

bool isValidReadFile( const char *fileName )
{
    FILE *file;
    int probe_char;

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
        fclose( file );
        return 0;
    }
    ungetc( probe_char, file );

    fclose( file );
    return 1;
}

bool readWriteCheck( const char *fileName, const bool checkWrite, const bool failIfError )
{
    FILE *file;
    file = fopen( fileName, checkWrite ? "w" : "r" );
    if ( file == NULL && failIfError )
    {
        string mode = checkWrite ? "writing." : "reading.";
        cerr << "Could not open file " << fileName << " for "
             << mode << " Aborting." << endl;
        exit( EXIT_FAILURE );
    }
    if ( file )
        fclose( file );
    return ( file != NULL );
}

int isDirectoryEmpty( const string &dirname ) // return: -1=directory doesn't exist, 0=not empty, 1=empty
{
    int n = 0;
    struct dirent *d;
    DIR *dir = opendir( dirname.c_str() );
    if ( dir == NULL ) //Not a directory or doesn't exist
        return -1;
    while ( ( d = readdir( dir ) ) != NULL )
    {
        if ( ++n > 2 )
            break;
    }
    closedir( dir );
    if ( n <= 2 ) //Directory Empty
        return 1;
    else
        return 0;
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

vector<string> splitString ( string s, const string &token )
{
    vector<string> vs;
    while ( s.find( token ) != string::npos )
    {
        vs.push_back( s.substr( 0, s.find( token ) ) );
        s = s.substr( s.find( token ) + ( token.length() ) );
    }
    vs.push_back( s );
    return vs;
}

bool isBwtFileCompressed( const string &filename )
{
    ifstream bwtFile( filename.c_str() );
    for ( int i = 0; i < 10; ++i )
    {
        char c = 'A';
        bwtFile.get( c );
        switch ( toupper( c ) )
        {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'N':
            case '$':
                break;
            default:
                return true;
        }
    }
    return false;
}

void detectInputBwtProperties( const string &prefix, vector<string> &filenames, bool &isBwtCompressed, string &availableFileLetters )
{
    // Detect {prefix}-B0* files
    for ( unsigned i = 0; i < alphabetSize; ++i )
    {
        stringstream filename;
        filename << prefix << "-B0" << i;
        Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Probing " << filename.str() << "..." << endl;
        if ( access( filename.str().c_str(), R_OK ) == -1 )
            break;
        filenames.push_back( filename.str() );
        Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Discovered " << filename.str() << endl;
    }

    // Detect {prefix}-*01 files
    for ( char c = 'A'; c <= 'Z' ; ++c )
    {
        stringstream filename;
        filename << prefix << "-" << c << "01";
        Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "Probing " << filename.str() << "..." << endl;
        if ( access( filename.str().c_str(), R_OK ) == 0 )
        {
            Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "Discovered " << filename.str() << endl;
            availableFileLetters += c;
        }
    }

    // Detect if files are run-length-encoded based on first few bytes
    isBwtCompressed = false;
    if ( !filenames.empty() )
    {
        for ( unsigned int j = 1; j < filenames.size() && !isBwtCompressed; ++j )
        {
            if ( isBwtFileCompressed( filenames[j] ) )
                isBwtCompressed = true;
        }

        // When BWT0x are RLE-encoded: Check that BWT0 file has the same encoding as others, to detect the old situation where BWT0 was always ASCII-encoded
        if ( isBwtCompressed )
        {
            bool isBwt0Compressed = isBwtFileCompressed( filenames[0] );
            if ( !isBwt0Compressed )
            {
                cerr << "ERROR: " << filenames[0] << " is ASCII-encoded, whereas the other BWT files are using a different encoding." << endl;
                cerr << " In this version of BEETL, all the files must have the same encoding." << endl;
                cerr << " You can convert this file with beetl-convert --input-format=bwt_ascii --output-format=bwt_rle." << endl;
                exit( -1 );
            }
        }
    }
    if ( isBwtCompressed )
    {
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "BWT files detected as RLE compressed" << endl;
    }
    else
    {
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE ) Logger::out() << "BWT files detected as ASCII" << endl;
    }
}

int safeRename( const string &from, const string &to )
{
    // renames/moves files even across partitions
    if ( rename( from.c_str(), to.c_str() ) )
    {
        Logger_if( LOG_SHOW_IF_VERY_VERBOSE )
        {
            perror( ( "Info: BCRexternalBWT: Error \"renaming\" file " + from + " to " + to ).c_str() );
            cerr << "Using mv command instead." << endl;
        }
        string cmd = "mv -f \""  + from + "\" \"" + to + "\"";
        system( cmd.c_str() );
    }
    return 0;
}

void pauseBetweenCycles()
{
    static int skip = 0;
    if ( skip )
    {
        clog << "Iteration complete. Still continuing for " << skip << " iteration." << endl;
        --skip;
    }
    else
    {
        fflush( 0 );
        clog << "Iteration complete" << endl;
        clog << " Press Return to continue, or enter a number of cycles to continue for..." << endl;
        string input;
        getline( cin, input );
        stringstream ss( input );
        ss >> skip;
        if ( skip )
            --skip;
    }
}

void readProcSelfStat( int &out_pid, int &out_num_threads, int &out_processor )
{
    using std::ios_base;
    using std::ifstream;
    using std::string;

    int tid = syscall( SYS_gettid ); //gettid();

    // 'file' stat seems to give the most reliable results
    ostringstream oss;
    oss << "/proc/" << tid << "/stat";
    ifstream stat_stream( oss.str().c_str(), ios_base::in );

    // dummy vars for leading entries in stat that we don't care about
    //
    string /*pid,*/ comm, state, ppid, pgrp, session, tty_nr, tpgid;
    string flags, minflt, cminflt, majflt, cmajflt, utime, stime, cutime;
    string cstime, priority, nice, /*num_threads,*/ itrealvalue, starttime, vsize, rss;
    string rsslim, startcode, endcode, startstack, kstkesp, kstkeip, signal, blocked;
    string sigignore, sigcatch, wchan, nswap, cnswap, exit_signal, /*processor,*/ rt_priority;
    string policy, delayacct_blkio_ticks;


    stat_stream >> out_pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >> tpgid
                >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime >> stime >> cutime
                >> cstime >> priority >> nice >> out_num_threads >> itrealvalue >> starttime >> vsize >> rss
                >> rsslim >> startcode >> endcode >> startstack >> kstkesp >> kstkeip >> signal >> blocked
                >> sigignore >> sigcatch >> wchan >> nswap >> cnswap >> exit_signal >> out_processor >> rt_priority
                >> policy >> delayacct_blkio_ticks;


    stat_stream.close();

    out_pid = tid;
    //   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    //   vm_usage     = vsize / 1024.0;
    //   resident_set = rss * page_size_kb;
}



shared_ptr<istream> openInputFileOrDashAsCin( const string &filename )
{
    shared_ptr<istream> input;
    if ( filename == "-" )
        input.reset( &cin, emptyDeleter() );
    else
        input.reset( new ifstream( filename.c_str() ) );
    return input;
}
