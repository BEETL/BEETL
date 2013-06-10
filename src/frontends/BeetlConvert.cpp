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

#include "BeetlConvert.hh"

#include "BwtReader.hh"
#include "BwtWriter.hh"
#include "Common.hh"
#include "LetterCount.hh"
#include "SeqReader.hh"
#include "TransposeFasta.hh"
#include "config.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/stat.h>

using namespace std;


void printUsage()
{
    cout << "Parameters:" << endl;
    cout << "    --input (-i)             Input file name or prefix" << endl;
    cout << "    --output (-o)            Output file name or prefix" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "    --input-format           = autodetect [fastq|fasta|seq|cyc|bcl|bwt_ascii|bwt_rle]" << endl;
    cout << "    --output-format          = autodetect [fastq|fasta|seq|cyc|bwt_ascii|bwt_rle]" << endl;
    cout << "    --help (-h)              Help" << endl;
    cout << endl;
    cout << "Formats:" << endl;
    cout << "    FASTQ     : 4 lines per read ('@'read name, nucleotides, '+', bases)" << endl;
    cout << "    FASTA     : 2 lines per read ('>'read name, nucleotides)" << endl;
    cout << "    seq       : 1 line  per read (nucleotides)" << endl;
    cout << "    cyc       : 1 file per cycle (byte K in file cycN  = Nth nucleotide of Kth read)" << endl;
    cout << "    BCL       : Conversion of a single cycle corresponding to 1 BCL file - only works with output-format=cyc" << endl;
    cout << "    BWT_ASCII : ASCII sequence of BWT-reordered nucleotides" << endl;
    cout << "    BWT_RLE   : Run-length-encoded version of BWT_ASCII, where bits 0-3 = binary-encoded nucleotide and bits 4-7 = count-1" << endl;
    cout << endl;
}

struct BeetlConvertArguments
{
    string argInput;
    string argOutput;
    string argInputFormat;
    string argOutputFormat;

    BeetlConvertArguments() {}
};


void convert( const BeetlConvertArguments &args )
{
    cout << "Conversion from " << args.argInputFormat << " to " << args.argOutputFormat << " (" << args.argInput << " -> " << args.argOutput << ")" << endl;

    if ( args.argInputFormat == args.argOutputFormat )
    {
        cerr << "Error: same input and output formats" << endl;
        exit ( 1 );
    }
    else if ( args.argInputFormat == "fasta" )
    {
        if ( args.argOutputFormat == "fastq" )
        {
            // FASTA -> FASTQ
            cerr << "Error: fasta->fastq needs extra qualities" << endl;
            exit ( 2 );
        }
        else if ( args.argOutputFormat == "seq" )
        {
            // FASTA -> SEQ
            ifstream inputStream( args.argInput.c_str() );
            ofstream outputStream( args.argOutput.c_str() );
            string str1, str2;
            while ( getline( inputStream, str1 ) &&
                    getline( inputStream, str2 ) )
            {
                assert( !str1.empty() && str1[0] == '>' );
                outputStream << str2 << "\n";
            }
            return;
        }
        else if ( args.argOutputFormat == "cyc" )
        {
            // FASTA -> CYC
            SeqReaderFile *pReader( SeqReaderFile::getReader( fopen( args.argInput.c_str(), "rb" ) ) );
            TransposeFasta trasp;
            trasp.init( pReader );
            trasp.convert( args.argInput, args.argOutput );
            delete pReader;
            return;
        }
        else if ( args.argOutputFormat == "bwt_ascii" || args.argOutputFormat == "bwt_rle" )
        {
            // FASTA -> BWT_*
            cerr << "Error: This is not a simple file conversion. Try \"beetl bwt\"" << endl;
            exit ( 1 );
        }
    }
    else if ( args.argInputFormat == "fastq" )
    {
        if ( args.argOutputFormat == "fasta" )
        {
            // FASTQ -> FASTA
            ifstream inputStream( args.argInput.c_str() );
            ofstream outputStream( args.argOutput.c_str() );
            string str1, str2, str3, str4;
            while ( getline( inputStream, str1 ) &&
                    getline( inputStream, str2 ) &&
                    getline( inputStream, str3 ) &&
                    getline( inputStream, str4 ) )
            {
                assert( !str1.empty() && str1[0] == '@' );
                outputStream << ">" << str1.substr( 1 ) << "\n";
                outputStream << str2 << "\n";
            }
            return;
        }
        else if ( args.argOutputFormat == "seq" )
        {
            // FASTQ -> SEQ
            ifstream inputStream( args.argInput.c_str() );
            ofstream outputStream( args.argOutput.c_str() );
            string str1, str2, str3, str4;
            while ( getline( inputStream, str1 ) &&
                    getline( inputStream, str2 ) &&
                    getline( inputStream, str3 ) &&
                    getline( inputStream, str4 ) )
            {
                assert( !str1.empty() && str1[0] == '@' );
                outputStream << str2 << "\n";
            }
            return;
        }
        else if ( args.argOutputFormat == "cyc" )
        {
            // FASTQ -> CYC
            SeqReaderFile *pReader( SeqReaderFile::getReader( fopen( args.argInput.c_str(), "rb" ) ) );
            TransposeFasta trasp;
            trasp.init( pReader );
            trasp.convert( args.argInput, args.argOutput );
            delete pReader;
            return;
        }
        else if ( args.argOutputFormat == "bwt_ascii" || args.argOutputFormat == "bwt_rle" )
        {
            // FASTQ -> BWT_*
            cerr << "Error: This is not a simple file conversion. Try \"beetl bwt\"" << endl;
            exit ( 1 );
        }
    }
    else if ( args.argInputFormat == "seq" )
    {
        if ( args.argOutputFormat == "fasta" )
        {
            // SEQ -> FASTA
            ifstream inputStream( args.argInput.c_str() );
            ofstream outputStream( args.argOutput.c_str() );
            int seqNum = 0;
            string str;
            while ( getline( inputStream, str ) )
            {
                outputStream << ">seq" << seqNum++ << "\n";
                outputStream << str << "\n";
            }
            return;
        }
        else if ( args.argOutputFormat == "fastq" )
        {
            // SEQ -> FASTQ
            cerr << "Error: seq->fastq needs extra qualities" << endl;
            exit ( 2 );
        }
        else if ( args.argOutputFormat == "cyc" )
        {
            // SEQ -> CYC
            SeqReaderFile *pReader( SeqReaderFile::getReader( fopen( args.argInput.c_str(), "rb" ) ) );
            TransposeFasta trasp;
            trasp.init( pReader );
            trasp.convert( args.argInput, args.argOutput );
            delete pReader;
            return;
        }
        else if ( args.argOutputFormat == "bwt_ascii" || args.argOutputFormat == "bwt_rle" )
        {
            // SEQ -> BWT_*
            cerr << "Error: This is not a simple file conversion. Try \"beetl bwt\"" << endl;
            exit ( 1 );
        }
    }
    else if ( args.argInputFormat == "cyc" )
    {
        if ( args.argOutputFormat == "fasta" )
        {
            // CYC -> FASTA
            TransposeFasta trasp;
            trasp.convertFromCycFileToFastaOrFastq( args.argInput, args.argOutput );
            return;
        }
        else if ( args.argOutputFormat == "fastq" )
        {
            // CYC -> FASTQ
            TransposeFasta trasp;
            trasp.convertFromCycFileToFastaOrFastq( args.argInput, args.argOutput );
            return;
        }
        else if ( args.argOutputFormat == "seq" )
        {
            // CYC -> SEQ
        }
        if ( args.argOutputFormat == "bcl" )
        {
            // CYC -> BCL
            // Single file conversion if the file exists, or multiple if %d appears in the middle
            int cycleNum = -1; // -1 = single file conversion

            for ( ;; )
            {
                // Open input file
                char inputFilename[1024];
                assert( args.argInput.size() < 1020 && "Input filename too long" );
                if ( cycleNum == -1 )
                    strcpy( inputFilename, args.argInput.c_str() );
                else
                    sprintf( inputFilename, args.argInput.c_str(), cycleNum );
                ifstream is( inputFilename, ios_base::binary );
                if ( !is.good() )
                {
                    if ( cycleNum == -1 )
                    {
                        cycleNum = 0;
                        continue;
                    }
                    else
                        break;
                }

                // Open input qual file
                string qualFilename = string( inputFilename ) + ".qual";
                ifstream isQual( qualFilename.c_str(), ios_base::binary );

                // Open output file
                char outputFilename[1024];
                assert( args.argOutput.size() < 1020 && "Output filename too long" );
                if ( cycleNum == -1 )
                    strcpy( outputFilename, args.argOutput.c_str() );
                else
                    sprintf( outputFilename, args.argOutput.c_str(), cycleNum + 1 );
                ofstream os( outputFilename, ios_base::binary );

                // Write header bytes
                struct stat fileStats;
                assert ( stat( inputFilename, &fileStats ) == 0 );
                unsigned int byteCount = fileStats.st_size;
                os.write( reinterpret_cast<char *>( &byteCount ), 4 );

                // Convert
                clog << "Converting: " << inputFilename << " + " << qualFilename << " -> " << outputFilename << endl;
                char c, q;
                vector<unsigned char> baseToBin( 256, 0xFF );
                baseToBin[( int )'A'] = 0;
                baseToBin[( int )'C'] = 1;
                baseToBin[( int )'G'] = 2;
                baseToBin[( int )'T'] = 3;
                while ( is.get( c ) )
                {
                    assert ( isQual.get( q ) );
                    unsigned char b = baseToBin[( int )c];
                    assert( b <= 3 );
                    assert( q >= 33 && q < ( 33 + 64 ) );
                    os.put( b | ( ( q - 33 ) << 2 ) );
                }

                if ( cycleNum == -1 )
                    break;
                else
                    ++cycleNum;
            }
            return;
        }
        else if ( args.argOutputFormat == "bwt_ascii" || args.argOutputFormat == "bwt_rle" )
        {
            // CYC -> BWT_*
            cerr << "Error: This is not a simple file conversion. Try \"beetl bwt\"" << endl;
            exit ( 1 );
        }
    }
    else if ( args.argInputFormat == "bcl" )
    {
        if ( args.argOutputFormat == "cyc" )
        {
            // BCL -> CYC
            string tileName = args.argInput.substr( 1 + args.argInput.find_last_of( "/" ) );
            tileName.erase( tileName.find_last_of( ".bcl" ) - 3 );

            string filterFilename = args.argInput;
            try
            {
                filterFilename.erase( filterFilename.find_last_of( "/" ) );
                filterFilename.erase( filterFilename.find_last_of( "/" ) + 1 );
            }
            catch ( const std::exception &e )
            {
                // Ignore... the file will get ignored anyway
            }
            filterFilename += tileName + ".filter";

            clog << "Converting tile " << tileName << " + " << filterFilename << endl;

            // Single file conversion if the file exists, or multiple if %d appears in the middle
            int cycleNum = -1; // -1 = single file conversion

            for ( ;; )
            {
                // Open input file
                char inputFilename[1024];
                assert( args.argInput.size() < 1020 && "Input filename too long" );
                if ( cycleNum == -1 )
                    strcpy( inputFilename, args.argInput.c_str() );
                else
                    sprintf( inputFilename, args.argInput.c_str(), cycleNum + 1 );
                ifstream is( inputFilename, ios_base::binary );
                if ( !is.good() )
                {
                    if ( cycleNum == -1 )
                    {
                        cycleNum = 0;
                        continue;
                    }
                    else
                        break;
                }

                // Open output file
                ostringstream oss;
                oss << args.argOutput;
                if ( cycleNum >= 0 )
                    oss << "." << cycleNum;
                string outputFilename = oss.str();
                ofstream os( outputFilename.c_str(), ios_base::binary );

                // Open output quality file
                ostringstream ossQual;
                ossQual << args.argOutput << ".qual";
                if ( cycleNum >= 0 )
                    ossQual << "." << cycleNum;
                string outputFilenameQual = ossQual.str();
                ofstream osQual( outputFilenameQual.c_str(), ios_base::binary );

                // Open input filter file
                ifstream filterFile( filterFilename.c_str(), ios_base::binary );

                // Skip header bytes of each file
                unsigned int buf4a, buf4b;
                filterFile.read( reinterpret_cast<char *>( &buf4a ), 4 );
                is.read( reinterpret_cast<char *>( &buf4b ), 4 );
                if ( buf4a != buf4b )
                {
                    // Try to read an extra 8 bytes, corresponding to a different version of the .filter format
                    filterFile.read( reinterpret_cast<char *>( &buf4a ), 4 );
                    filterFile.read( reinterpret_cast<char *>( &buf4a ), 4 );
                    if ( buf4a != buf4b )
                    {
                        cerr << "Error parsing filter file: File lengths are not matching" << endl;
                        exit( -1 );
                    }
                }

                // Convert
                clog << "Converting: " << inputFilename << " -> " << outputFilename << " + " << outputFilenameQual << endl;
                char c, filterVal;
                const char binToBase[4] = { 'A', 'C', 'G', 'T' };
                while ( is.get( c ) )
                {
                    if ( !filterFile.get( filterVal ) || filterVal != '\0' )
                    {
                        os.put( binToBase[c & 3] );
                        osQual.put( 33 + ( ( ( unsigned char )c ) >> 2 ) );
                    }
                }

                if ( cycleNum == -1 )
                    break;
                else
                    ++cycleNum;
            }
            return;
        }
    }
    else if ( args.argInputFormat == "bwt_ascii" )
    {
        if ( args.argOutputFormat == "bwt_rle" )
        {
            // BWT_ASCII -> BWT_RLE
            BwtReaderBase *pReader = new BwtReaderASCII( args.argInput );
            BwtWriterBase *pWriter = new BwtWriterRunLength( args.argOutput );

            while ( pReader->readAndSend( *pWriter, 1000000000 ) > 0 ) {}

            delete pReader;
            delete pWriter;
            return;
        }
        else if ( args.argOutputFormat == "fasta" || args.argOutputFormat == "fastq" || args.argOutputFormat == "seq" || args.argOutputFormat == "cyc" )
        {
            // BWT_ASCII -> FASTA|FASTQ|SEQ|CYC
            cerr << "Error: This is not a simple file conversion. Try \"beetl unbwt\"" << endl;
            exit ( 1 );
        }
    }
    else if ( args.argInputFormat == "bwt_rle" )
    {
        if ( args.argOutputFormat == "bwt_ascii" )
        {
            // BWT_RLE -> BWT_ASCII
            BwtReaderBase *pReader = new BwtReaderRunLength( args.argInput );
            BwtWriterBase *pWriter = new BwtWriterASCII( args.argOutput );

            while ( pReader->readAndSend( *pWriter, 1000000000 ) > 0 ) {}

            delete pReader;
            delete pWriter;
            return;
        }
        else if ( args.argOutputFormat == "fasta" || args.argOutputFormat == "fastq" || args.argOutputFormat == "seq" || args.argOutputFormat == "cyc" )
        {
            // BWT_RLE -> FASTA|FASTQ|SEQ|CYC
            cerr << "Error: This is not a simple file conversion. Try \"beetl unbwt\"" << endl;
            exit ( 1 );
        }
    }

    cerr << "Error: unknown file conversion" << endl;
    exit ( 1 );
}


int main( const int argc, const char **argv )
{
    // Generated using: http://patorjk.com/software/taag/#p=display&f=Soft&t=BEETL%20convert
    cout << ",-----.  ,------.,------.,--------.,--.                                                    ,--.   " << endl;
    cout << "|  |) /_ |  .---'|  .---''--.  .--'|  |        ,---. ,---. ,--,--,,--.  ,--.,---. ,--.--.,-'  '-. " << endl;
    cout << "|  .-.  \\|  `--, |  `--,    |  |   |  |       | .--'| .-. ||      \\\\  `'  /| .-. :|  .--''-.  .-' " << endl;
    cout << "|  '--' /|  `---.|  `---.   |  |   |  '--.    \\ `--.' '-' '|  ||  | \\    / \\   --.|  |     |  |   " << endl;
    cout << "`------' `------'`------'   `--'   `-----'     `---' `---' `--''--'  `--'   `----'`--'     `--'   " << endl;
    cout << "Version " << PACKAGE_VERSION << endl;
    cout << endl;

    cout << "Command called:" << endl << "   ";
    for ( int i = 0; i < argc; ++i )
    {
        cout << " " << argv[i];
    }
    cout << "\n" << endl;

    BeetlConvertArguments args;

    for ( int i = 1; i < argc; ++i )
    {
        if ( isNextArgument( "-h", "--help", argc, argv, i ) )
        {
            printUsage();
            exit( 0 );
        }
        else if ( isNextArgument( "-i", "--input"                  , argc, argv, i, &args.argInput               ) ) {}
        else if ( isNextArgument( ""  , "--input-format"           , argc, argv, i, &args.argInputFormat         ) ) {}
        else if ( isNextArgument( "-o", "--output"                 , argc, argv, i, &args.argOutput              ) ) {}
        else if ( isNextArgument( ""  , "--output-format"          , argc, argv, i, &args.argOutputFormat        ) ) {}
        else
        {
            cerr << "Error: Invalid parameter: " << argv[i] << endl;
            printUsage();
            exit( 1 );
        }
    }

    // Checking for required parameters
    if ( args.argInput.empty() || args.argOutput.empty() )
    {
        cerr << "Error: Missing arguments: --input and --output are required.\n" << endl;
        printUsage();
        exit( 1 );
    }

    // Auto-detection of missing arguments
    if ( args.argInputFormat.empty() )
    {
        args.argInputFormat = detectFileFormat( args.argInput );
    }
    checkFileFormat( args.argInput, args.argInputFormat );

    if ( args.argOutputFormat.empty() )
    {
        args.argOutputFormat = detectFileFormat( args.argOutput );
    }
    checkFileFormat( args.argOutput, args.argOutputFormat );

    //    checkIfAlreadyExistingFile( args.argOutput );

    // Convert args to lower case
    std::transform( args.argInputFormat.begin(), args.argInputFormat.end(), args.argInputFormat.begin(), ::tolower );
    std::transform( args.argOutputFormat.begin(), args.argOutputFormat.end(), args.argOutputFormat.begin(), ::tolower );

    convert( args );

    return 0;
}
