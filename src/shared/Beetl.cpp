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

#include "Beetl.hh" // own declarations

#include "Algorithm.hh"   // framework class declaration
#include "BCRext.hh"    // interface to BCRext
#include "BWTCollection.hh" // interface to BCR
#include "CountWords.hh" // interface to countwords

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

using namespace std;


//#define BEETL_ID "$Id$"
const string BEETL_ID( "1.0" );


int main( int numArgs, char **args )
{

    if ( numArgs < 2 )
    {
        print_usage( args[0] );
        exit( EXIT_SUCCESS );
    }

    if ( strcmp( args[1], COMMAND_BCR ) == 0 )
    {

        bcrMode = 0; // explicitly set default mode to zero (ie build BWT from sequence)
        CompressionFormatType bcrCompression( compressionASCII );
        // start at 2, 0 is the executable, 1 the command name parsed above
        for ( int i = 2; i < numArgs; i++ )

            if ( args[i][0] == '-' ) // only flags here "-X etc."
            {

                switch ( args[i][1] )
                {
                    case 'i':
                        // next param should be the filename, checking now
                        isArgumentOrExit( i + 1, numArgs );
                        // don't check file here! In -m 1 mode this is a prefix name, not an actual file
                        // so the file doesn't need to exist for input to be valid - TC 4.9.12
                        // fileIsReadableOrExit(args[i + 1]);
                        bcrFileIn = args[i + 1]; // this should be the name
                        cout << "-> input file is " << bcrFileIn << endl;
                        break;
                    case 'o':
                        isArgumentOrExit( i + 1, numArgs );
                        bcrFileOut = args[i + 1];
                        cout << "-> output prefix is " << bcrFileOut << endl;
                        break;
                    case 'm':
                        isArgumentOrExit( i + 1, numArgs );
                        bcrMode = atoi( args[i + 1] );
                        if ( bcrMode > 2 || bcrMode < 0 )
                        {
                            cerr << bcrMode << " is no valid bcr mode " << endl;
                            exit( EXIT_FAILURE );
                        }
                        cout << "-> working mode set to \""
                             << bcrModes[bcrMode]
                             << "\""
                             << endl;
                        break;
                    case 'a':
                        bcrCompression = compressionASCII;;
                        cout << "-> writing ASCII encoded output"
                             << endl;
                        break;
                    case 'h':
                        bcrCompression = compressionHuffman;
                        cout << "Huffman encoding not yet supported, sorry."
                             << endl;
                        exit( EXIT_FAILURE );
                        //cout << "-> writing huffman encoded output"
                        //        << endl;
                        break;
                    case 'r':
                        bcrCompression = compressionRunLength;
                        cout << "-> writing runlength encoded output"
                             << endl;
                        break;
                    case 't':
                        bcrCompression = compressionIncrementalRunLength;
                        cout << "-> writing incremental runlength encoded output"
                             << endl;
                        break;
                    default:
                        cout << "!! unknown flag \""
                             << args[i][1] << "\"" << endl;
                        print_usage( args[0], COMMAND_BCR );
                        exit( EXIT_FAILURE );
                }
            }

        // check if all arguments are given
        if ( bcrFileIn.length() > 0 && bcrMode >= 0 )
        {

            if ( bcrFileOut.length() == 0 )
            {
                bcrFileOut = bcrFileOutPrefixDefault;
            }

            if ( bcrMode == 0 )
            {
                fileIsReadableOrExit( bcrFileIn.c_str() );
            }

            // Redundant and produces confusing warning - SeqReader now checks for and handles fasta, fastq and raw
            // sequence files, the latter two trigger a warning here. TC 4.9.12
            //   if (bcrMode == 0)
            //  {
            //     isValidFastaFile(bcrFileIn.c_str());
            // }

            // created new tool object
            Algorithm *pBCR
                = new BCR( bcrMode, bcrFileIn, bcrFileOut, bcrCompression );

            // run previous main method
            pBCR->run();

            // clean up
            delete pBCR;

            // die
            exit( EXIT_SUCCESS );
        }
        else
        {
            // something wrong happened
            print_usage( args[0], COMMAND_BCR );
            exit( EXIT_FAILURE );
        }

    }
    else if ( strcmp( args[1], COMMAND_BCR_EXT ) == 0 )
    {


        // set defaults for BCRext mode
        bcrExtAsciiOutput = false; // use normal ASCII alphabet as output
        bcrExtHuffmanOutput = false; // use huffman encoding as compression
        bcrExtRunlengthOutput = true; // use RunLength encoding [default]
        bcrExtImplicitSort = false; // do implicit sort of input sequences


        for ( int i = 2; i < numArgs; i++ )
        {
            if ( args[i][0] == '-' ) // only flags here "-X etc."
            {

                string thisArg( ( string )args[i] );
                if ( thisArg == "-i" )
                {
                    // next param should be the filename, checking now
                    isArgumentOrExit( i + 1, numArgs );
                    fileIsReadableOrExit( args[i + 1] );
                    bcrExtFileIn = args[i + 1]; // this should be the name
                    cout << "-> input file is " << bcrExtFileIn << endl;
                }
                else if ( thisArg == "-p" )
                {
                    isArgumentOrExit( i + 1, numArgs );
                    bcrExtFileOutPrefix = args[i + 1];
                    cout << "-> output prefix set to "
                         << bcrExtFileOutPrefix << endl;
                }
                else if ( thisArg == "-s" )
                {
                    cout << "-> using SeqFile input"
                         << endl;
                    bcrExtUseSeq = true;
                }
                else if ( thisArg == "-a" )
                {
                    bcrExtAsciiOutput = true;
                    bcrExtRunlengthOutput = false;
                    bcrExtHuffmanOutput = false;
                    cout << "-> writing ASCII encoded output"
                         << endl;
                }
                else if ( thisArg == "-h" )
                {
                    bcrExtAsciiOutput = false;
                    bcrExtRunlengthOutput = false;
                    bcrExtHuffmanOutput = true;
                    cout << "-> writing huffman encoded output"
                         << endl;
                }
                else if ( thisArg == "-r" )
                {
                    bcrExtAsciiOutput = false;
                    bcrExtRunlengthOutput = true;
                    bcrExtHuffmanOutput = false;
                    cout << "-> writing runlength encoded output"
                         << endl;
                }
                else if ( thisArg == "-sap" )
                {
                    bcrExtImplicitSort = true;
                    cout << "-> perform implicit sort of input sequences"
                         << endl;
                }
                else
                {
                    cout << "!! unknown flag \""
                         << args[i][1] << "\"" << endl;
                    print_usage( args[0], COMMAND_BCR_EXT );
                    exit( EXIT_FAILURE );
                }
            } // ~if begins with -
        } // ~for

        //   cout << bcrExtImplicitSort << bcrExtHuffmanOutput << bcrExtRunlengthOutput << endl;

        if ( bcrExtFileIn.empty() )
        {
            print_usage( args[0], COMMAND_BCR_EXT );
            exit( EXIT_FAILURE );
        }

        /*
                if (!bcrExtUseSeq && ! isValidFastaFile(bcrExtFileIn.c_str())) {
                    print_usage(args[0],COMMAND_BCR_EXT);
        //            exit(EXIT_FAILURE);
                }
                if (bcrExtUseSeq && !isValidReadFile(bcrExtFileIn.c_str())) {
                    print_usage(args[0],COMMAND_BCR_EXT);
                    exit(EXIT_FAILURE);
                }
        */

        // check if all arguments are given
        if ( bcrExtRunlengthOutput || bcrExtAsciiOutput || bcrExtHuffmanOutput )

            // check if its a fasta file when fasta flag set
            //&& ( isValidFastaFile(bcrExtFileIn.c_str()) ||

            // is this a valid .seq file? only when no fasta flag
            //(bcrExtUseSeq && isValidReadFile(bcrExtFileIn.c_str()))) )
        {

            if ( bcrExtFileOutPrefix.length() == 0 )
            {
                bcrExtFileOutPrefix = bcrExtFileOutPrefixDefault;
            }

            // created new tool object
            Algorithm *pBCRext = new BCRext( bcrExtHuffmanOutput,
                                             bcrExtRunlengthOutput,
                                             bcrExtAsciiOutput,
                                             bcrExtImplicitSort,
                                             bcrExtUseSeq,
                                             bcrExtFileIn,
                                             bcrExtFileOutPrefix );

            // run previous main method
            pBCRext->run();

            // clean up
            delete pBCRext;

            // die
            exit( EXIT_SUCCESS );
        }
        else
        {
            // something wrong happened
            print_usage( args[0], COMMAND_BCR_EXT );
            exit( EXIT_FAILURE );
        }

    }
    else if ( strcmp( args[1], COMMAND_COUNTWORDS ) == 0 )
    {
        vector<string> filesA, filesB, filesC;
        string ncbiTax;
        int minWord( 0 );
        string prefix;// prefix for output files
        for ( int i = 2; i < numArgs; i++ )

            if ( args[i][0] == '-' ) // only flags here "-X etc."
            {
                //cout << "blah" << endl;
                switch ( args[i][1] )
                {
                    case 'r':
                        cout << "-> reference genome mode for set B " << endl;
                        whichHandler = 'r';
                        break;
                    case 's' :
                        cout << "-> splice mode set " << endl;
                        whichHandler = 's';
                        break;
                    case 'm' :
                        cout << "-> metagenomic mode set " << endl;
                        whichHandler = 'm';
                        break;
                    case 'a':
                        while ( args[++i][0] != '-' )
                        {
                            //cout << args[i] << " fred " << endl;
                            fileIsReadableOrExit( args[i] );
                            filesA.push_back( args[i] );
                            cout << "-> input file A is "
                                 << filesA.back()
                                 << endl;
                        }
                        i--;
                        break;
#ifdef OLD
                        // next param should be the filename, checking
                        isArgumentOrExit( i + 1, numArgs );
                        fileIsReadableOrExit( args[i + 1] );
                        countWordsInputA = args[i + 1]; // should be the name
                        //   cout << "-> input file A is "
                        //   << countWordsInputA
                        //     << endl;
                        break;
#endif//OLD
                    case 'b':
                        while ( ( ++i ) != numArgs )
                        {
                            if ( strcmp( args[i], "-c" ) == 0 )
                            {
                                i--;
                                break;
                            }
                            fileIsReadableOrExit( args[i] );
                            filesB.push_back( args[i] );
                            //       cout << "-> input file B is "
                            //   << filesB.back()
                            //   << endl;
                        }
                        //        i--;
                        break;
                    case 'c' :
                        if ( whichHandler != 'm' )
                        {
                            cerr << "Merged set C not needed Here!" << endl;
                            break;
                        }
                        while ( ( ++i ) != numArgs )
                        {
                            if ( args[i][0] == '-' )
                            {
                                i--;
                                break;
                            }
                            fileIsReadableOrExit( args[i] );
                            filesC.push_back( args[i] );
                            cout << "-> input file C is "
                                 << filesC.back()
                                 << endl;
                        }
                        break;
#ifdef OLD
                        // next param should be the filename, checking now
                        isArgumentOrExit( i + 1, numArgs );
                        fileIsReadableOrExit( args[i + 1] );
                        countWordsInputB = args[i + 1]; // should be the name
                        cout << "-> input file B is "
                             << countWordsInputB
                             << endl;
                        break;
#endif
                    case 't':
                        if ( whichHandler != 'm' )
                        {
                            cerr << "The taxonomic file is not needed without the -m option." << endl;
                            break;
                        }
                        cout << "-> Got taxonomic information " << args[i + 1] << endl ;
                        ncbiTax = args[i + 1];
                        i++;
                        break;
                    case 'd':
                        if ( whichHandler != 'm' )
                        {
                            cerr << "Option -d is not available without merged input. " << endl;
                            break;
                        }
                        cout << "Will test database " << endl;
                        testDatabase = true;
                        break;
                    case 'w':
                        minWord = atoi( args[i + 1] );
                        cout << "-> minimum word length w in metagenomic mode set to \""
                             << minWord
                             << "\""
                             << endl;
                        break;
                    case 'p':
                        prefix  = args[i + 1];
                        break;
                    case 'k':
                        isArgumentOrExit( i + 1, numArgs );
                        maxLengthK = atoi( args[i + 1] );
                        if ( maxLengthK < 0 )
                        {
                            cerr << "!! "
                                 << maxLengthK
                                 << " is no valid length "
                                 << endl;
                            exit( EXIT_FAILURE );
                        }
                        cout << "-> max kmer length set to \""
                             << maxLengthK
                             << "\""
                             << endl;
                        break;
                    case 'n':
                        isArgumentOrExit( i + 1, numArgs );
                        minimalOccurencesN = atoi( args[i + 1] );
                        if ( minimalOccurencesN < 0 )
                        {
                            cerr << "!! "
                                 << minimalOccurencesN
                                 << " is no valid value "
                                 << endl;
                            exit( EXIT_FAILURE );
                        }
                        cout << "-> min occurences n set to \""
                             << minimalOccurencesN
                             << "\""
                             << endl;
                        break;
                    case 'A':
                        cout << "-> assuming set A is compressed" << endl;
                        compressedInputA = true;
                        break;
                    case 'B':
                        cout << "-> assuming set B is compressed" << endl;
                        compressedInputB = true;
                        break;
                    case 'C':
                        cout << "-> assuming set A&B are compressed" << endl;
                        compressedInputA = true;
                        compressedInputB = true;
                        break;
                    default:
                        cout << "!! unknown flag \"" << args[i][1]
                             << "\"" << endl;
                        print_usage( args[0], COMMAND_COUNTWORDS );
                        exit( EXIT_FAILURE );
                }
            }

        // check for required arguments
        if ( ( maxLengthK > 0 ) && ( minimalOccurencesN > 0 ) &&
             ( !filesA.empty() ) && ( filesA.size() == filesB.size() ) )
        {

            // create new tool object
            Algorithm *pcountWords = new CountWords( compressedInputA,
                    compressedInputB, whichHandler, minimalOccurencesN,
                    maxLengthK, filesA, filesB, filesC, ncbiTax, testDatabase, minWord, "" );

            // run the "main" method
            pcountWords->run();

            // clean up
            delete pcountWords;

            // closing time
            exit( EXIT_SUCCESS );
        }
        else
        {
            // oops
            print_usage( args[0], COMMAND_COUNTWORDS );
            exit( EXIT_FAILURE );
        }

    }
    else
    {
        cerr << "!! \"" << args[1] << "\": unknown command" << endl;
        print_usage( args[0] );
        exit( EXIT_FAILURE );
    }
    return 0;
}

void fileIsReadableOrExit( string filename )
{

    FILE *pFile;
    // test file for read access
    pFile = fopen( filename.c_str(), "r" );

    if ( pFile != NULL )
    {
        fclose( pFile );
        return;
    }
    else
    {
        cerr << "!! \"" << filename << "\" is NOT readable!" << endl;
        exit( EXIT_FAILURE );
    }
}

void isArgumentOrExit( int num, int numArgs )
{
    if ( ( num ) > ( numArgs - 1 ) )
    {
        cerr << "!! CLI parsing error. Wrong number of arguments?" << endl;
        exit( EXIT_FAILURE );
    }
}


void print_usage( char *args, const char *command )
{
    cerr << endl << "- This is the BEETL software library -" << endl
         << endl
         // Tony 13.6.12 - BEETL_ID is not informative now we have moved to git
         //            << "Framework version" << endl
         //            << BEETL_ID << endl
         << endl;

    if ( command == 0 )
        cerr << "Included in this framework are the following algorithms" << endl
             << endl
             << endl;

    if ( command == 0 || string( command ) == COMMAND_BCR_EXT )
        cerr << "-> BCRext - command \"" << COMMAND_BCR_EXT << "\"" << endl
             << "========================================================" << endl
             << "improved version of the original algorithm" << endl
             << "uses significantly less RAM (a.k.a. none) but depends heavily on I/O" << endl
             << endl
             << "Usage: " << args << " "
             << COMMAND_BCR_EXT << " -i <read file> -p <output file prefix> [-h -r -a] [-s] [-sap]" << endl
             << endl
             << "-i <file>:\tinput file in fasta format" << endl
             << "-s:\t\tuse .seq input files instead of fasta (each line one sequence)" << endl
             << "-p <string>:\toutput file names will start with \"prefix\"" << endl
             << "-a:\t\toutput ASCII encoded files" << endl
             << "-r:\t\toutput runlength encoded files [recommended]" << endl
             << "-h:\t\toutput Huffman encoded files" << endl
             << "-sap:\t\tperform implicit permutation of collection to obtain more compressible BWT"
             << endl
             << endl;

    if ( command == 0 || string( command ) == COMMAND_BCR )
        cerr << "-> BCR - command \"" << COMMAND_BCR << "\"" << endl
             << "========================================================" << endl
             << "original algorithm to construct the BWT of a set of reads" << endl
             << "needs approximately 14GB of RAM for 1 billion reads" << endl
             << endl
             << "Usage: " << args << " "
             << COMMAND_BCR << " -i <fasta or seq read file> -o <output file> [-r -t -a] -m <[0,1,2]>" << endl
             << endl
             << "-i <file>:\tinput set of reads [if mode = 1 set the prefix of the BWT files, normally BCR-B0]" << endl
             << "-o <file>:\toutput file" << endl
             << "-m <n>:\t\tmode = 0 --> BCR " << endl
             << "\t\tmode = 1 --> unBCR " << endl
             << "\t\tmode = 2 --> Backward search + Locate SeqID (Uses extra file called \"searchedKmers\" as input)" << endl
             << "-a:\t\toutput ASCII encoded files in BCR mode [mode = 0]" << endl
             << "-r:\t\toutput run length encoded files in BCR mode [mode = 0]" << endl
             << "-t:\t\toutput incremental runlength encoded files [experimental]" << endl
             << endl
             << endl;

    if ( command == 0 || string( command ) == COMMAND_COUNTWORDS )
        cerr << "-> countWords - command \"" << COMMAND_COUNTWORDS << "\"" << endl
             << "========================================================" << endl
             << "find all words of length at least k that occur" << endl
             << "at least n times in string set A and never in string set B" << endl
             << endl
             << "Usage: " << args << " "
             << COMMAND_COUNTWORDS << " [-A -B -C] [-s -r -m] -k <n> -n <n> -p <file_prefix> -a <set A> -b <set B> "
             << "[-c <C part of set B> -t file -d -o databaseTestOut ]" << endl
             << endl
             << "-A:\t\tassume BWT files for set A are in compressed format" << endl
             << "-B:\t\tassume BWT files for set B are in compressed format" << endl
             << "-C:\t\tassume BWT files for sets A and B are in compressed format" << endl
             << "-s:\t\tassume set B are reads of a reference" << endl
             << "-r:\t\tassume set B is a reference genome" << endl
             << "-m:\t\tassume set B are merged reference genomes" << endl
             << "-k <n>:\t\tmaximal length" << endl
             << "-n <n>:\t\tminimum number of occurences (coverage)" << endl
             << "-p <prefix>:\tprefix out output files" << endl
             << "-a <file>:\tinput set A" << endl
             << "-b <file>:\tinput set B" << endl
             << "-c <file>:\tC output of merging, only needed in -m Mode" << endl
             << "-t file:\ttaxonomy information for the merged input, only needed in -m Mode" << endl
             << "-w n\t\tminimal length k-mer" << endl
             << "-d:\t\tflag to test the minimal needed word length for the different taxa in the databse. Only available with -m" << endl
             << endl
             << endl;

    if ( !command )
        cerr << "If you had fun using these algorithms you may cite:" << endl
             << "---------------------------------------------------" << endl
             << "Markus J. Bauer, Anthony J. Cox and Giovanna Rosone" << endl
             << "Lightweight BWT Construction for Very Large String Collections. " << endl
             << "Proceedings of CPM 2011, pp.219-231, doi: 10.1007/978-3-642-21458-5_20" << endl
             << "[Description of BWT construction algorithms in 'bcr' and 'ext' modes]" << endl << endl
             << "Markus J. Bauer, Anthony J. Cox and Giovanna Rosone" << endl
             << "Lightweight algorithms for constructing and inverting the BWT of string collections "
             << endl << "Theoretical Computer Science, doi: 10.1016/j.tcs.2012.02.002" << endl
             << "[As above plus description of BWT inversion in 'bcr' mode]" << endl << endl
             << "Anthony J. Cox, Markus J. Bauer, Tobias Jakobi and Giovanna Rosone" << endl
             << "Large-scale compression of genomic sequence databases with the Burrows-Wheeler transform"
             << endl
             << "Bioinformatics, doi: 10.1093/bioinformatics/bts173" << endl
             << "[Description of '-sap' compression boosting strategy in 'ext' mode]" << endl << endl
             << "BEETL web page:" << endl
             << "---------------" << endl
             << "http://beetl.github.com/BEETL/" << endl << endl
             << "BEETL user group:" << endl
             << "-----------------" << endl
             << "http://tech.groups.yahoo.com/group/BEETL/" << endl << endl;
}
