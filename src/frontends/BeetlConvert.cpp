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
#include "libzoo/cli/Common.hh"
#include "libzoo/io/Bcl.hh"
#include "LetterCount.hh"
#include "SeqReader.hh"
#include "TransposeFasta.hh"
#include "config.h"
#include "parameters/ConvertParameters.hh"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/stat.h>

#include "../redist/gzstream.hh"

using namespace std;
using namespace BeetlConvertParameters;


ConvertParameters params;
static const char binToBase[4] = { 'A', 'C', 'G', 'T' };

void printUsage()
{
    params.printUsage();

    cout << "Formats:" << endl;
    cout << "    FASTQ     : 4 lines per read ('@'read name, nucleotides, '+', bases)" << endl;
    cout << "    FASTA     : 2 lines per read ('>'read name, nucleotides)" << endl;
    cout << "    seq       : 1 line  per read (nucleotides)" << endl;
    cout << "    cyc       : 1 file per cycle (byte K in file cycN  = Nth nucleotide of Kth read)" << endl;
    cout << "    BCL       : Conversion of a single cycle corresponding to 1 BCL file - only works with output-format=cyc" << endl;
    cout << "    RunFolder : Illumina Run Folder containing BCL files" << endl;
    cout << "    BWT_ASCII : ASCII sequence of BWT-reordered nucleotides" << endl;
    cout << "    BWT_RLE   : Run-length-encoded version of BWT_ASCII, where bits 0-3 = binary-encoded nucleotide and bits 4-7 = count-1" << endl;
    cout << endl;
}

void outputSequenceConstrainedWithSequenceLength( ofstream &outputStream, string &str2, const ParameterEntry &sequenceLength )
{
    // Warning: may modify str2
    if ( sequenceLength.isSet() && str2.size() != sequenceLength )
    {
        if ( str2.size() < sequenceLength )
            outputStream << string( sequenceLength - str2.size(), 'N' ); // pre-pad with 'N'
        else
            str2.erase( sequenceLength );
    }
    outputStream << str2 << '\n';
}

void outputLineUsingMissingDataFrom( ofstream &outputStream, ifstream &missingDataFile, ostringstream &ss )
{
    string str;
    if ( missingDataFile.good() && getline( missingDataFile, str ) )
        outputStream << str << '\n';
    else
        outputStream << ss.str() << '\n';
}

void skipLineInMissingDataFile( ifstream &missingDataFile )
{
    string str;
    if ( missingDataFile.good() )
        getline( missingDataFile, str );
}

void launchBeetlConvert()
{
    cout << "Conversion from " << params.getStringValue( "input format" ) << " to " << params.getStringValue( "output format" ) << " (" << params.getStringValue( "input filename" ) << " -> " << params.getStringValue( "output filename" ) << ")" << endl;

    ifstream missingDataFile;
    if ( params["use missing data from"].isSet() )
    {
        missingDataFile.open( params.getStringValue( "use missing data from" ).c_str() );
    }

    if ( params.getStringValue( "input format" ) == params.getStringValue( "output format" ) && params[ "input format" ] != INPUT_FORMAT_FASTQ )
    {
        cerr << "Error: same input and output formats" << endl;
        exit ( 1 );
    }
    else if ( params["input format"] == INPUT_FORMAT_FASTA )
    {
        if ( params["output format"] == OUTPUT_FORMAT_FASTQ )
        {
            // FASTA -> FASTQ
            cerr << "Error: fasta->fastq needs extra qualities" << endl;
            exit ( 2 );
        }
        else if ( params["output format"] == OUTPUT_FORMAT_SEQ )
        {
            // FASTA -> SEQ
            ifstream inputStream( params.getStringValue( "input filename" ).c_str() );
            ofstream outputStream( params.getStringValue( "output filename" ).c_str() );
            string str1, str2;
            while ( getline( inputStream, str1 ) &&
                    getline( inputStream, str2 ) )
            {
                assert( !str1.empty() && str1[0] == '>' );
                outputStream << str2 << '\n';
            }
            return;
        }
        else if ( params["output format"] == OUTPUT_FORMAT_CYC )
        {
            // FASTA -> CYC
            SeqReaderFile *pReader( SeqReaderFile::getReader( fopen( params.getStringValue( "input filename" ).c_str(), "rb" ) ) );
            TransposeFasta trasp;
            trasp.init( pReader );
            trasp.convert( params.getStringValue( "input filename" ), params.getStringValue( "output filename" ), false );
            delete pReader;
            return;
        }
        else if ( params["output format"] == OUTPUT_FORMAT_BWT_ASCII || params["output format"] == OUTPUT_FORMAT_BWT_RLE )
        {
            // FASTA -> BWT_*
            cerr << "Error: This is not a simple file conversion. Try \"beetl bwt\"" << endl;
            exit ( 1 );
        }
    }
    else if ( params["input format"] == INPUT_FORMAT_FASTQ )
    {
        if ( params["output format"] == OUTPUT_FORMAT_FASTA )
        {
            // FASTQ -> FASTA
            ifstream inputStream( params.getStringValue( "input filename" ).c_str() );
            ofstream outputStream( params.getStringValue( "output filename" ).c_str() );
            string str1, str2, str3, str4;
            while ( getline( inputStream, str1 ) &&
                    getline( inputStream, str2 ) &&
                    getline( inputStream, str3 ) &&
                    getline( inputStream, str4 ) )
            {
                assert( !str1.empty() && str1[0] == '@' );
                outputStream << ">" << str1.substr( 1 ) << '\n';
                outputSequenceConstrainedWithSequenceLength( outputStream, str2, params["sequence length"] );
            }
            return;
        }
        else if ( params["output format"] == OUTPUT_FORMAT_SEQ )
        {
            // FASTQ -> SEQ
            ifstream inputStream( params.getStringValue( "input filename" ).c_str() );
            ofstream outputStream( params.getStringValue( "output filename" ).c_str() );
            string str1, str2, str3, str4;
            while ( getline( inputStream, str1 ) &&
                    getline( inputStream, str2 ) &&
                    getline( inputStream, str3 ) &&
                    getline( inputStream, str4 ) )
            {
                assert( !str1.empty() && str1[0] == '@' );
                outputSequenceConstrainedWithSequenceLength( outputStream, str2, params["sequence length"] );
            }
            return;
        }
        else if ( params["output format"] == OUTPUT_FORMAT_CYC )
        {
            // FASTQ -> CYC
            SeqReaderFile *pReader( SeqReaderFile::getReader( fopen( params.getStringValue( "input filename" ).c_str(), "rb" ) ) );
            TransposeFasta trasp;
            trasp.init( pReader );
            trasp.convert( params.getStringValue( "input filename" ), params.getStringValue( "output filename" ), false );
            delete pReader;
            return;
        }
        else if ( params["output format"] == OUTPUT_FORMAT_BWT_ASCII || params["output format"] == OUTPUT_FORMAT_BWT_RLE )
        {
            // FASTQ -> BWT_*
            cerr << "Error: This is not a simple file conversion. Try \"beetl bwt\"" << endl;
            exit ( 1 );
        }
        else if ( params["output format"] == OUTPUT_FORMAT_FASTQ )
        {
            // FASTQ -> FASTQ, only available with remove-padding
            if ( params["remove padding"] == false )
            {
                cerr << "Error: FASTQ->FASTQ is only available with --remove-padding" << endl;
                exit( 1 );
            }
            ifstream inputStream( params.getStringValue( "input filename" ).c_str() );
            ofstream outputStream( params.getStringValue( "output filename" ).c_str() );
            string str1, str2, str3, str4;
            while ( getline( inputStream, str1 ) &&
                    getline( inputStream, str2 ) &&
                    getline( inputStream, str3 ) &&
                    getline( inputStream, str4 ) )
            {
                assert( !str1.empty() && str1[0] == '@' );
                outputStream << str1 << '\n';
                assert( str2.size() == str4.size() && "Bases and Qualities must have the same length" );
                if ( !str2.empty() )
                {
                    int firstValidIndex = 0;
                    int lastValidIndex = str2.size() - 1;
                    while ( firstValidIndex < str2.size() && str2[firstValidIndex] == 'N' )
                        ++firstValidIndex;
                    while ( lastValidIndex >= 0 && str2[lastValidIndex] == 'N' )
                        --lastValidIndex;
                    if ( firstValidIndex > lastValidIndex )
                    {
                        clog << "Warning: Read full of 'N': " << str1 << endl;
                        str2 = str2[0];
                        str4 = str4[0];
                    }
                    else
                    {
                        str2 = str2.substr( firstValidIndex, lastValidIndex - firstValidIndex + 1 );
                        str4 = str4.substr( firstValidIndex, lastValidIndex - firstValidIndex + 1 );
                    }
                }
                outputStream << str2 << '\n';
                outputStream << str3 << '\n';
                outputStream << str4 << '\n';
            }
            return;
        }
    }
    else if ( params["input format"] == INPUT_FORMAT_SEQ )
    {
        if ( params["output format"] == OUTPUT_FORMAT_FASTA )
        {
            // SEQ -> FASTA
            ifstream inputStream( params.getStringValue( "input filename" ).c_str() );
            ofstream outputStream( params.getStringValue( "output filename" ).c_str() );
            int seqNum = 0;
            string str;
            while ( getline( inputStream, str ) )
            {
                // line 1
                ostringstream ss;
                ss << ">seq" << seqNum++;
                outputLineUsingMissingDataFrom( outputStream, missingDataFile, ss );

                // line 2
                outputStream << str << '\n';
                skipLineInMissingDataFile( missingDataFile );
            }
            return;
        }
        else if ( params["output format"] == OUTPUT_FORMAT_FASTQ )
        {
            // SEQ -> FASTQ
            if ( !params["use missing data from"].isSet() )
            {
                cerr << "Error: seq->fastq needs extra qualities (see --use-missing-data-from)" << endl;
                exit( 2 );
            }
            ifstream inputStream( params.getStringValue( "input filename" ).c_str() );
            ofstream outputStream( params.getStringValue( "output filename" ).c_str() );
            int seqNum = 0;
            string str;
            while ( getline( inputStream, str ) )
            {
                // line 1
                ostringstream ss1;
                ss1 << ">seq" << seqNum++;
                outputLineUsingMissingDataFrom( outputStream, missingDataFile, ss1 );

                // line 2
                outputStream << str << '\n';
                skipLineInMissingDataFile( missingDataFile );

                // line 3
                outputStream << "+\n";
                skipLineInMissingDataFile( missingDataFile );

                // line 4
                ostringstream ss4; // empty as we don't know the qualities
                outputLineUsingMissingDataFrom( outputStream, missingDataFile, ss4 );
            }
            return;
        }
        else if ( params["output format"] == OUTPUT_FORMAT_CYC )
        {
            // SEQ -> CYC
            SeqReaderFile *pReader( SeqReaderFile::getReader( fopen( params.getStringValue( "input filename" ).c_str(), "rb" ) ) );
            TransposeFasta trasp;
            trasp.init( pReader );
            trasp.convert( params.getStringValue( "input filename" ), params.getStringValue( "output filename" ), false );
            delete pReader;
            return;
        }
        else if ( params["output format"] == OUTPUT_FORMAT_BWT_ASCII || params["output format"] == OUTPUT_FORMAT_BWT_RLE )
        {
            // SEQ -> BWT_*
            cerr << "Error: This is not a simple file conversion. Try \"beetl bwt\"" << endl;
            exit ( 1 );
        }
    }
    else if ( params["input format"] == INPUT_FORMAT_CYC )
    {
        if ( params["output format"] == OUTPUT_FORMAT_FASTA )
        {
            // CYC -> FASTA
            TransposeFasta trasp;
            trasp.convertFromCycFileToFastaOrFastq( params.getStringValue( "input filename" ), params.getStringValue( "output filename" ), false );
            return;
        }
        else if ( params["output format"] == OUTPUT_FORMAT_FASTQ )
        {
            // CYC -> FASTQ
            TransposeFasta trasp;
            trasp.convertFromCycFileToFastaOrFastq( params.getStringValue( "input filename" ), params.getStringValue( "output filename" ), false );
            return;
        }
        else if ( params["output format"] == OUTPUT_FORMAT_SEQ )
        {
            // CYC -> SEQ
        }
        if ( params["output format"] == OUTPUT_FORMAT_BCL )
        {
            // CYC -> BCL
            // Single file conversion if the file exists, or multiple if %d appears in the middle
            int cycleNum = -1; // -1 = single file conversion

            for ( ;; )
            {
                // Open input file
                char inputFilename[1024];
                assert( params.getStringValue( "input filename" ).size() < 1020 && "Input filename too long" );
                if ( cycleNum == -1 )
                    strcpy( inputFilename, params.getStringValue( "input filename" ).c_str() );
                else
                    sprintf( inputFilename, params.getStringValue( "input filename" ).c_str(), cycleNum );
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
                assert( params.getStringValue( "output filename" ).size() < 1020 && "Output filename too long" );
                if ( cycleNum == -1 )
                    strcpy( outputFilename, params.getStringValue( "output filename" ).c_str() );
                else
                    sprintf( outputFilename, params.getStringValue( "output filename" ).c_str(), cycleNum + 1 );
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
        else if ( params["output format"] == OUTPUT_FORMAT_BWT_ASCII || params["output format"] == OUTPUT_FORMAT_BWT_RLE )
        {
            // CYC -> BWT_*
            cerr << "Error: This is not a simple file conversion. Try \"beetl bwt\"" << endl;
            exit ( 1 );
        }
    }
    else if ( params["input format"] == INPUT_FORMAT_BCL || params["input format"] == INPUT_FORMAT_BCL_GZ )
    {
        if ( params["output format"] == OUTPUT_FORMAT_CYC )
        {
            // BCL -> CYC
            string tileName = params.getStringValue( "input filename" ).substr( 1 + params.getStringValue( "input filename" ).rfind( "/" ) );

            // Convert names like s_1_xx.bcl to s_1_00xx.bcl
            int parsedLaneNum = 0;
            int parsedTileNum = 0;
            if ( sscanf( tileName.c_str(), "s_%d_%d.bcl", &parsedLaneNum, &parsedTileNum ) == 2 )
            {
                char buf[100];
                sprintf( buf, "s_%d_%04d", parsedLaneNum, parsedTileNum );
                tileName = string( buf );
            }
            else
            {
                // otherwise just remove the extension
                tileName.erase( tileName.rfind( ".bcl" ) );
            }

            string filterPath = params.getStringValue( "input filename" );
            string filterFilename;
            for ( int subdirdepth = 0; subdirdepth <= 2; ++subdirdepth )
            {
                try
                {
                    size_t lastSlashPos = filterPath.rfind( "/" );
                    if ( filterPath == "" || filterPath == "." )
                        filterPath = "..";
                    else if ( filterPath == "/" )
                        filterPath = "/";
                    else if ( lastSlashPos == 0 )
                        filterPath = "/";
                    else if ( lastSlashPos == string::npos )
                        filterPath = "";
                    else
                        filterPath.erase( lastSlashPos );
                }
                catch ( const std::exception &e )
                {
                    // Ignore... the file will get ignored anyway
                }
                filterFilename = ( filterPath.empty() ? "" : ( filterPath + "/" ) ) + tileName + ".filter";
                if ( readWriteCheck( filterFilename.c_str(), false, false ) )
                    break;
            }

            clog << "Converting tile " << tileName << " + " << filterFilename << endl;

            // Single file conversion if the file exists, or multiple if %d appears in the middle
            int cycleNum = -1; // -1 = single file conversion

            for ( ;; )
            {
                // Open input file
                char inputFilename[1024];
                assert( params.getStringValue( "input filename" ).size() < 1020 && "Input filename too long" );
                if ( cycleNum == -1 )
                    strcpy( inputFilename, params.getStringValue( "input filename" ).c_str() );
                else
                    sprintf( inputFilename, params.getStringValue( "input filename" ).c_str(), cycleNum + 1 );

                istream *is;
                if ( params["input format"] == INPUT_FORMAT_BCL_GZ )
                {
                    is = new igzstream( inputFilename );
                }
                else
                {
                    is = new ifstream( inputFilename, ios_base::binary );
                }
                if ( !is->good() )
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
                oss << params.getStringValue( "output filename" );
                if ( cycleNum >= 0 )
                    oss << "." << cycleNum;
                string outputFilename = oss.str();
                ofstream os( outputFilename.c_str(), ios_base::binary );

                // Open output quality file
                ostringstream ossQual;
                ossQual << params.getStringValue( "output filename" ) << ".qual";
                if ( cycleNum >= 0 )
                    ossQual << "." << cycleNum;
                string outputFilenameQual = ossQual.str();
                ofstream osQual( outputFilenameQual.c_str(), ios_base::binary );

                // Open input filter file
                ifstream filterFile( filterFilename.c_str(), ios_base::binary );

                // Skip header bytes of each file
                unsigned int buf4a, buf4b;
                filterFile.read( reinterpret_cast<char *>( &buf4a ), 4 );
                is->read( reinterpret_cast<char *>( &buf4b ), 4 );
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
                while ( is->get( c ) )
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
    else if ( params["input format"] == INPUT_FORMAT_RUNFOLDER )
    {
        if ( params["output format"] == OUTPUT_FORMAT_CYC )
        {
            // RunFolder -> CYC
            BclRunFolder bclRunFolder( params.getStringValue( "input filename" ), "", "" );
            uint cycleCount = bclRunFolder.getCycleCount();

            // Open output files
            vector<ofstream *> outCyc, outCycQual;
            for ( uint i = 0; i < cycleCount; ++i )
            {
                int cycleNum = i;

                // Open output bases file
                ostringstream oss;
                oss << params.getStringValue( "output filename" ) << "." << cycleNum;
                string outputFilename = oss.str();
                outCyc.push_back( new ofstream( outputFilename.c_str(), ios_base::binary ) );

                // Open output quality file
                ostringstream ossQual;
                ossQual << params.getStringValue( "output filename" ) << ".qual" << "." << cycleNum;
                string outputFilenameQual = ossQual.str();
                outCycQual.push_back( new ofstream( outputFilenameQual.c_str(), ios_base::binary ) );
            }

            string nextLane;
            string nextTile;
            while ( bclRunFolder.getNextLaneAndTileNames( nextLane, nextTile ) )
            {
                cout << "Processing lane " << nextLane << " tile " << nextTile << endl;
                bclRunFolder.initReader( nextLane, nextTile, 1, cycleCount );
                vector< unsigned int > bclValues;
                int count = 0;
                while ( bclRunFolder.getRead( bclValues ) )
                {
                    //                    cout << "value: " << bclValues[0] << endl;
                    ++count;
                    for ( uint i = 0; i < cycleCount; ++i )
                    {
                        char c = bclValues[i];
                        outCyc[i]->put( binToBase[c & 3] );
                        outCycQual[i]->put( 33 + ( ( ( unsigned char )c ) >> 2 ) );
                    }
                }
                cout << "processed " << count << " reads" << endl;
            }

            // Closing all files
            for ( uint i = 0; i < cycleCount; ++i )
            {
                delete outCyc[i];
                delete outCycQual[i];
            }
            outCyc.clear();
            outCycQual.clear();

            cout << "done" << endl;
            return;
        }
    }
    else if ( params["input format"] == INPUT_FORMAT_BWT_ASCII )
    {
        if ( params["output format"] == OUTPUT_FORMAT_BWT_RLE )
        {
            // BWT_ASCII -> BWT_RLE
            BwtReaderBase *pReader = new BwtReaderASCII( params.getStringValue( "input filename" ) );
            BwtWriterBase *pWriter = new BwtWriterRunLength( params.getStringValue( "output filename" ) );

            while ( pReader->readAndSend( *pWriter, 1000000000 ) > 0 ) {}

            delete pReader;
            delete pWriter;
            return;
        }
        else if ( params["output format"] == OUTPUT_FORMAT_FASTA || params["output format"] == OUTPUT_FORMAT_FASTQ || params["output format"] == OUTPUT_FORMAT_SEQ || params["output format"] == OUTPUT_FORMAT_CYC )
        {
            // BWT_ASCII -> FASTA|FASTQ|SEQ|CYC
            cerr << "Error: This is not a simple file conversion. Try \"beetl unbwt\"" << endl;
            exit ( 1 );
        }
    }
    else if ( params["input format"] == INPUT_FORMAT_BWT_RLE )
    {
        if ( params["output format"] == OUTPUT_FORMAT_BWT_ASCII )
        {
            // BWT_RLE -> BWT_ASCII
            BwtReaderBase *pReader = new BwtReaderRunLength( params.getStringValue( "input filename" ) );
            BwtWriterBase *pWriter = new BwtWriterASCII( params.getStringValue( "output filename" ) );

            while ( pReader->readAndSend( *pWriter, 1000000000 ) > 0 ) {}

            delete pReader;
            delete pWriter;
            return;
        }
        else if ( params["output format"] == OUTPUT_FORMAT_FASTA || params["output format"] == OUTPUT_FORMAT_FASTQ || params["output format"] == OUTPUT_FORMAT_SEQ || params["output format"] == OUTPUT_FORMAT_CYC )
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

    if ( !params.parseArgv( argc, argv ) || params["help"] == 1 || !params.chechRequiredParameters() )
    {
        printUsage();
        exit( 1 );
    }

    // Auto-detection of missing arguments
    if ( !params["input format"].isSet() )
    {
        const string &filename = params.getStringValue( "input filename" );
        string fileFormat = detectFileFormat( filename );
        if ( fileFormat.empty() )
        {
            cerr << "Error: file format not recognised for " << filename << endl;
            exit( -1 );
        }
        params["input format"] = fileFormat;
    }
    checkFileFormat( params.getStringValue( "input filename" ), params["input format"] );

    if ( !params["output format"].isSet() )
    {
        const string &filename = params.getStringValue( "output filename" );
        string fileFormat = detectFileFormat( filename );
        if ( fileFormat.empty() )
        {
            cerr << "Error: file format not recognised for " << filename << endl;
            exit( -1 );
        }
        params["output format"] = fileFormat;
    }
    checkFileFormat( params.getStringValue( "output filename" ), params["output format"] );

    // Use default parameter values where needed
    params.commitDefaultValues();

    //    checkIfAlreadyExistingFile( params.getStringValue("output filename") );

    // Check for unsupported cases
    if ( params["sequence length"].isSet() )
    {
        if ( ! ( params["input format"] == INPUT_FORMAT_FASTQ && ( params["output format"] == OUTPUT_FORMAT_SEQ || params["output format"] == OUTPUT_FORMAT_FASTA ) ) )
        {
            cerr << "Error: --sequence-length is currently only implemented for FASTQ->SEQ and FASTQ->FASTA conversions." << endl;
            exit( 1 );
        }
    }

    // Check for unsupported cases
    if ( params["remove padding"] == true )
    {
        if ( ! ( params["input format"] == INPUT_FORMAT_FASTQ && params["output format"] == OUTPUT_FORMAT_FASTQ ) )
        {
            cerr << "Error: --remove-padding is currently only implemented for FASTQ->FASTQ conversions." << endl;
            exit( 1 );
        }
    }

    launchBeetlConvert();

    return 0;
}
