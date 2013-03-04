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

#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;


void printUsage()
{
    cout << "Parameters:" << endl;
    cout << "    --input (-i)             Input file name or prefix" << endl;
    cout << "    --output (-o)            Output file name or prefix" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "    --input-format           = autodetect [fastq|fasta|seq|cyc|bwt_ascii|bwt_rle]" << endl;
    cout << "    --output-format          = autodetect [fastq|fasta|seq|cyc|bwt_ascii|bwt_rle]" << endl;
    cout << "    --help (-h)              Help" << endl;
    cout << endl;
    cout << "Formats:" << endl;
    cout << "    FASTQ     : 4 lines per read ('@'read name, nucleotides, '+', bases)" << endl;
    cout << "    FASTA     : 2 lines per read ('>'read name, nucleotides)" << endl;
    cout << "    seq       : 1 line  per read (nucleotides)" << endl;
    cout << "    cyc       : 1 file per cycle (byte K in file cycN  = Nth nucleotide of Kth read)" << endl;
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
            SeqReaderFile *pReader( SeqReaderFile::getReader( fopen( args.argInput.c_str(),"rb" ) ) );
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
            SeqReaderFile *pReader( SeqReaderFile::getReader( fopen( args.argInput.c_str(),"rb" ) ) );
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
            SeqReaderFile *pReader( SeqReaderFile::getReader( fopen( args.argInput.c_str(),"rb" ) ) );
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
        else if ( args.argOutputFormat == "bwt_ascii" || args.argOutputFormat == "bwt_rle" )
        {
            // CYC -> BWT_*
            cerr << "Error: This is not a simple file conversion. Try \"beetl bwt\"" << endl;
            exit ( 1 );
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
    for ( int i=0; i < argc; ++i )
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

    convert( args );

    return 0;
}
