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

#include "BCRexternalBWT.hh"

#include "BWTCollection.hh"
#include "Filename.hh"
#include "Tools.hh"
#include "TransposeFasta.hh"
#include "parameters/BwtParameters.hh"
#include "parameters/SearchParameters.hh"
#include "parameters/UnbwtParameters.hh"
#include "libzoo/util/Logger.hh"
#include "libzoo/util/TemporaryFilesManager.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <sys/stat.h>

using namespace std;
using SXSI::BWTCollection;


#define SIZEBUFFER 1024
#define DIMBLOCK 2048

////////////////////////////////////////////////////////////////////////////
// Class BCRexternalBWT

/**
 * Constructor inits
 */
BCRexternalBWT::BCRexternalBWT ( const string &file1, const string &fileOutput, const int mode, const CompressionFormatType outputCompression, ToolParameters *toolParams )
    : toolParams_( toolParams, emptyDeleter() )
    , bwtParams_( 0 )
    , unbwtParams_( 0 )
    , searchParams_( 0 )
{
    const char *intermediateCycFiles = "cyc.";
    if ( mode == 0 )
    {
        using namespace BeetlBwtParameters;
        bwtParams_ = dynamic_pointer_cast<BwtParameters>( toolParams_ );
        if ( !bwtParams_ )
        {
            // Legacy mode: old code provides only outputCompression. We create a BwtParameters structure to let it work.
            bwtParams_.reset( new BwtParameters );
            switch ( outputCompression )
            {
                case compressionASCII:
                    ( *bwtParams_ )[PARAMETER_INTERMEDIATE_FORMAT] = INTERMEDIATE_FORMAT_ASCII;
                    ( *bwtParams_ )[PARAMETER_OUTPUT_FORMAT] = OUTPUT_FORMAT_ASCII;
                    break;
                case compressionRunLength:
                    ( *bwtParams_ )[PARAMETER_INTERMEDIATE_FORMAT] = INTERMEDIATE_FORMAT_RLE;
                    ( *bwtParams_ )[PARAMETER_OUTPUT_FORMAT] = OUTPUT_FORMAT_RLE;
                    break;
                case compressionIncrementalRunLength:
                    ( *bwtParams_ )[PARAMETER_INTERMEDIATE_FORMAT] = INTERMEDIATE_FORMAT_MULTIRLE;
                    ( *bwtParams_ )[PARAMETER_OUTPUT_FORMAT] = OUTPUT_FORMAT_RLE;
                    break;
                case compressionHuffman:
                    ( *bwtParams_ )[PARAMETER_INTERMEDIATE_FORMAT] = INTERMEDIATE_FORMAT_HUFFMAN;
                    ( *bwtParams_ )[PARAMETER_OUTPUT_FORMAT] = OUTPUT_FORMAT_HUFFMAN;
                    break;
                default:
                    assert( false && "shouldn't reach here" );
            }
        }

        std::cerr << "Start BCR encode\n";

        outputCompression_ = ( CompressionFormatType )( -1 ); // todo: remove this variable altogether
        cerr << "Compression format for intermediate BWT files: " << bwtParams_->getStringValue( PARAMETER_INTERMEDIATE_FORMAT ) << endl;
        cerr << "Compression format for BWT output files: "       << bwtParams_->getStringValue( PARAMETER_OUTPUT_FORMAT ) << endl;

        //added by GIOVANNA
        if ( BUILD_SA == 1 )
            std::cerr << "Compute also the SA by using the BCR (BCR_SA)\n";
        else
            std::cerr << "Compute only the BWT by using the BCR  (BCR_SA)\n";

        if ( bwtParams_->getValue( PARAMETER_GENERATE_ENDPOSFILE ) || BUILD_SA )
        {
            // we make sure that this file does not exist, to avoid reading an old version
            Filename fileEndPos( bwtParams_->getStringValue( "output filename" ).c_str(), "-end-pos" );
            remove( fileEndPos );
        }

        int result = buildBCR( file1, intermediateCycFiles, bwtParams_.get() );
        checkIfNotEqual( result, 0 );
        bool hasProcessedQualities = ( result == 2 );

        if ( bwtParams_->getValue( PARAMETER_CONCATENATE_OUTPUT ) == true )
        {
            //Store the entire BWT from alphabetSize-files
            storeEntireBWT( fileOutput );

            if ( bwtParams_->getValue( PARAMETER_GENERATE_LCP ) == true )
            {
                storeEntireLCP( fileOutput );
            }
        }

        //Do we want compute the extended suffix array (position and number of sequence)?
        if ( BUILD_SA == 1 )  //To store the SA
        {
            storeEntirePairSA( fileOutput.c_str() );
            storeEntireSAfromPairSA( fileOutput.c_str() );
        }

        if ( verboseEncode == 1 )
        {
            if ( ( bwtParams_->getValue( PARAMETER_GENERATE_LCP ) == true ) && ( BUILD_SA == 1 ) )
            {
                std::cerr << "Store the files containing the BWT, LCP and SA in a single file\n";
                TmpFilename fnSA(       fileOutput, ".sa" );
                TmpFilename fnPairSA(   fileOutput, ".pairSA" );
                TmpFilename fnLCP(      fileOutput, ".lcp" );
                Filename fileOutRes( fileOutput, ".txt" );

                FILE *InFileBWT = fopen( fileOutput.c_str(), "rb" );
                if ( InFileBWT == NULL )
                {
                    std::cerr << "Entire BWT file: Error opening "  << fileOutput << std::endl;
                    exit ( EXIT_FAILURE );
                }


                FILE *InFilePairSA = fopen( fnPairSA, "rb" );
                if ( InFilePairSA == NULL )
                {
                    std::cerr << "Entire Pairs SA file: Error opening " << fnPairSA.str() << std::endl;
                    exit ( EXIT_FAILURE );
                }

                FILE *InFileSA = fopen( fnSA, "rb" );
                if ( InFileSA == NULL )
                {
                    std::cerr << "Entire SA file: Error opening " << fnSA << std::endl;
                    exit ( EXIT_FAILURE );

                }

                FILE *InFileLCP = fopen( fnLCP, "rb" );
                if ( InFileLCP == NULL )
                {
                    std::cerr << "Entire LCP file: Error opening " << fnLCP << std::endl;
                    exit ( EXIT_FAILURE );

                }

                ofstream outFile( fileOutRes.str() );
                if ( outFile.bad() )
                {
                    std::cerr << "Error opening output \"" << fileOutRes << "\" file" << std::endl;
                    exit ( EXIT_FAILURE );
                }

                uchar *bufferBWT = new uchar[SIZEBUFFER];
                uchar *bufferLCP = new uchar[SIZEBUFFER];
                ElementType *buffer = new ElementType[SIZEBUFFER];
                LetterNumber *bufferNChar = new LetterNumber[SIZEBUFFER];

                while ( ( !feof( InFileBWT ) ) && ( !feof( InFileSA ) ) && ( !feof( InFilePairSA ) ) && ( !feof( InFileLCP ) ) )
                {
                    LetterNumber numcharBWT = fread( bufferBWT, sizeof( uchar ), SIZEBUFFER, InFileBWT );
                    LetterNumber numcharPairSA = fread( buffer, sizeof( ElementType ), SIZEBUFFER, InFilePairSA );
                    LetterNumber numcharSA = fread( bufferNChar, sizeof( LetterNumber ), SIZEBUFFER, InFileSA );
                    LetterNumber numcharLCP = fread( bufferNChar, sizeof( LetterNumber ), SIZEBUFFER, InFileLCP );
                    //std::cerr << "Char read: " << numcharBWT  << "\t" << numcharLCP  << "\t" << numcharSA << "\t" << numcharPairSA  << "\n";
                    outFile << "bwt\tlcp\tpos\tnumSeq\tSA\n";
                    if ( ( numcharPairSA != numcharSA ) || ( numcharLCP != numcharSA ) || ( numcharLCP != numcharBWT ) )
                        std::cerr << "Error: number in BWT in Pair SA in SA and in LCP\n";
                    else
                    {
                        for ( LetterNumber i = 0; i < numcharSA; i++ )
                        {
                            outFile << bufferBWT[i]
                                    << '\t' << bufferLCP[i]
                                    << '\t' << buffer[i].sa
                                    << '\t' << buffer[i].numSeq
                                    << '\t' << bufferNChar[i]
                                    << endl;
                        }
                    }
                }
                delete[] buffer;
                delete[] bufferNChar;

                fclose( InFilePairSA );
                fclose( InFileSA );
                fclose( InFileLCP );
            }
        }

        TemporaryFilesManager::get().cleanupAllFiles();

        Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Removing/Renaming the BWT segments\n";
        for ( int g = 0 ; g < alphabetSize; g++ )
        {
            TmpFilename filename( g );
            if ( deletePartialBWT == 1 )
            {
                if ( remove( filename ) != 0 )
                    perror( ( "BCRexternalBWT: Error deleting file " + filename.str() ).c_str() );
            }
            else //rename the aux bwt file
            {
                Filename newFilename( fileOutput, "-B0", g, "" );
                safeRename( filename, newFilename );

                if ( hasProcessedQualities )
                {
                    TmpFilename qualFilename( "", g, ".qual" );
                    Filename newQualFilename( fileOutput, "-Q0", g, "" );
                    safeRename( qualFilename, newQualFilename );
                }
            }
        }

        if ( bwtParams_->getValue( PARAMETER_GENERATE_LCP ) == true )
        {
            std::cerr << "Removing/Renaming the LCP segment files\n";
            for ( AlphabetSymbol g = 0 ; g < alphabetSize; g++ )
            {
                TmpFilename filenameIn( "", g, ".lcp" );
                if ( deletePartialLCP == 1 )
                {
                    if ( remove( filenameIn ) != 0 )
                        perror( ( "BCRexternalBWT: Error deleting file " + filenameIn.str() ).c_str() );
                }
                else   //rename the aux lcp file
                {
                    Filename newFilename( fileOutput, "-L0", g, "" );
                    safeRename( filenameIn, newFilename );
                }
            }
        }

        /*  std::cerr << "Removing/Renaming the SA segments\n";
                for (AlphabetSymbol g = 0 ; g < alphabetSize; g++) {
                    Filename filenameIn( "sa_", g );
                    if (deletePartialSA == 1)  {
                        if (remove(filenameIn)!=0)
                            std::cerr << "BCRexternalBWT: Error deleting file" << std::endl;
                    }
                    else //renome the aux bwt file
                    {
                        Filename newfilename( fileOutput, filenameIn.str() );
                        if(rename(filenameIn, newfilename))
                            std::cerr  <<"BCRexternalBWT: Error renaming file" << std::endl;
                    }
                }
        */

        Logger_if( LOG_FOR_DEBUGGING )
        {
            if ( /*bwtParams_->getValue( PARAMETER_INTERMEDIATE_STORAGE_MEDIUM ) == INTERMEDIATE_STORAGE_MEDIUM_RAM*/
                bwtParams_->getValue( PARAMETER_INTERMEDIATE_FORMAT ) == INTERMEDIATE_FORMAT_MULTIRLE
            )
            {
                Logger::out() << "RAM file lengths:";
                extern vector< vector<unsigned char> > ramFiles; // declared in BwtWriter; todo: move those to another a new header file
                size_t totalFileLengths = 0;
                for ( unsigned int i = 0; i < ramFiles.size(); ++i )
                {
                    Logger::out() << " " << i << ":" << ramFiles[i].size();
                    totalFileLengths += ramFiles[i].size();
                }

                Logger::out() << std::endl;
                Logger::out() << "  total RAM used for files: " << totalFileLengths << std::endl;
            }
        }

    }
    else if ( mode == 1 )
    {
        using namespace BeetlUnbwtParameters;
        unbwtParams_ = dynamic_pointer_cast<UnbwtParameters>( toolParams_ );
        if ( unbwtParams_ == NULL )
        {
            // Legacy mode: old code was set using #defines. We create a BwtParameters structure to reflect the default values.
            unbwtParams_.reset( new UnbwtParameters );
            ( *unbwtParams_ )[PARAMETER_DECODE_DIRECTION] = DECODE_DIRECTION_BACKWARD;
            ( *unbwtParams_ )[PARAMETER_USE_VECTOR] = USE_VECTOR_ON;
        }

        std::cerr << "Start BCR decode\n";
        const char *fileOutBwt = "";
        int result = unbuildBCR( file1.c_str(), fileOutBwt, intermediateCycFiles, fileOutput.c_str() );
        checkIfEqual( result, 1 );
    }
    else if ( mode == 2 )
    {
        using namespace BeetlSearchParameters;
        searchParams_ = dynamic_pointer_cast<SearchParameters>( toolParams_ );
        if ( searchParams_ == NULL )
        {
            // Legacy mode: old code was set using #defines. We create a BwtParameters structure to reflect the default values.
            searchParams_.reset( new SearchParameters );
        }

        std::cerr << "Start Locate Function:\n";
        std::cerr << "Backward Search and Recovery of the number of sequences\n";
        const char *fileOutBwt = "";

        vector<string> kmers;
        char kmer[1000];
        SequenceLength lenKmer = 0;

        FILE *InFileKmer = fopen( "searchedKmers", "rb" );
        if ( InFileKmer == NULL )
        {
            std::cerr << "Error opening \"searchedKmers\" file" << std::endl;
            exit( EXIT_FAILURE );
        }
        while ( fgets( kmer, sizeof( kmer ), InFileKmer ) )
        {
            char *tmp = strchr( kmer, '\n' );
            if ( tmp )
                *tmp = '\0';
            tmp = strchr( kmer, '\r' );
            if ( tmp )
                *tmp = '\0';

            if ( ( strcmp( kmer, "\r" ) != 0 ) && ( strcmp( kmer, "\n" ) != 0 )
                 && ( strcmp( kmer, "\0" ) != 0 ) )
            {
                kmers.push_back( kmer );
                lenKmer = strlen( kmer );
            }
        }
        fclose( InFileKmer );



        vector <int> seqID;
        int result = SearchAndLocateKmer( file1.c_str(), fileOutBwt, intermediateCycFiles, kmers, lenKmer, seqID );
        checkIfEqual( result, 1 );
        std::cerr << "\nBCRexternalBWT: We have located all kmers, Now we store the positions of the found kmers" << endl;
        if ( seqID.empty() )
            std::cerr << "BCRexternalBWT: None of the k-mers occur in the collection" << endl;
        else
        {
            Filename newfilename( fileOutput );
            FILE *FilePosKmers = fopen( newfilename, "wb" );
            if ( FilePosKmers == NULL )
            {
                std::cerr << "BCRexternalBWT: could not open file " << newfilename << "!" << std::endl;
                exit ( EXIT_FAILURE );
            }
            fprintf( FilePosKmers, "kmer_ID \t N_kmer \t pos_in_the_SA\n" );
            SequenceNumber numTotKmers = 0;
            for ( SequenceNumber g = 0 ; g < kmers.size(); g++ )
            {
                for ( LetterNumber  j = FirstVector[g].posN ; j <= LastVector[g].posN; j++ ) //For each position between first and last
                {
                    fprintf( FilePosKmers, "%u \t %u \t %d\n", LastVector[g].seqN, numTotKmers, seqID[numTotKmers] );
                    numTotKmers++;
                }
            }
            fclose( FilePosKmers );

            /* if (verboseDecode == 1) {
                    SequenceNumber numTotKmers=0;
                    for (SequenceNumber g = 0 ; g < kmers.size(); g++) {
                        std::cerr << "\nk-mer of index " << LastVector[g].seqN << ": "<< kmers[LastVector[g].seqN] << ". Number of occurrences " << LastVector[g].posN-FirstVector[g].posN+1 << std::endl;
                        for (LetterNumber  j = FirstVector[g].posN ; j <= LastVector[g].posN; j++) { //For each position between first and last
                            std::cerr << "Number " << numTotKmers << " pos in the SA=\t"<< j << "\t SeqId=\t" << seqID[numTotKmers] <<  std::endl;
                            numTotKmers++;
                        }
                    }
                }
                */

        }
    }
    else
        std::cerr << "Mode Error" << endl;
}

int BCRexternalBWT::SearchAndLocateKmer ( char const *file1, char const *fileOutBwt, char const *fileOut, vector<string> kmers, SequenceLength lenKmer, vector <int> &seqID )
{
    LetterNumber freq[256];  //contains the distribution of the symbols.
    int resultInit = initializeUnbuildBCR( file1, fileOutBwt, freq );
    checkIfEqual( resultInit, 1 );

    std::cerr << "Frequency"  << "\n";
    for ( AlphabetSymbol i = 0; i < 255; i++ )
        if ( freq[i] > 0 )
        {
            std::cerr << i << "\t" << freq[i] << "\t" << ( int )alpha[i] << "\t" << ( int )alphaInverse[( int )alpha[i]] << "\n";
        }

    assert( unbwtParams_ || searchParams_ );
    if ( ( unbwtParams_ && unbwtParams_->getValue( BeetlUnbwtParameters::PARAMETER_USE_VECTOR ) == BeetlUnbwtParameters::USE_VECTOR_ON )
         || searchParams_ )
    {
        resultInit = computeVectorUnbuildBCR( file1, fileOutBwt, freq );
        checkIfEqual( resultInit, 1 );
    }

    std::cerr << "backwardSearchManyBCR\n";

    backwardSearchManyBCR( file1, fileOutBwt, fileOut, kmers, lenKmer );

    for ( SequenceNumber g = 0 ; g < kmers.size(); g++ )
    {
        std::cerr << "The number of occurrences of " << kmers[LastVector[g].seqN] << " is \t" << LastVector[g].posN - FirstVector[g].posN + 1 << "\n";
    }

    if ( verboseDecode == 1 )
    {
        std::cerr << "First and Last: "  <<  "\n";
        std::cerr << "Q  ";
        for ( SequenceNumber g = 0 ; g < kmers.size(); g++ )
        {
            std::cerr << ( int )FirstVector[g].pileN << " " << ( int )LastVector[g].pileN << "\t";
        }
        std::cerr << std::endl;
        std::cerr << "P  ";
        for ( SequenceNumber g = 0 ; g < kmers.size(); g++ )
        {
            std::cerr << FirstVector[g].posN  << " " << LastVector[g].posN  << "\t";
        }
        std::cerr << std::endl;
        std::cerr << "N  ";
        for ( SequenceNumber g = 0 ; g < kmers.size(); g++ )
        {
            std::cerr << FirstVector[g].seqN  << " " << LastVector[g].seqN  << "\t" ;
        }
        std::cerr << std::endl;
    }

    vector<int> tmpSeqId = recoverNSequenceForward( file1, fileOutBwt, kmers.size() );
    seqID.swap( tmpSeqId );
    if ( seqID.empty() )
        std::cerr << "\nSearchAndLocateKmer: No k-mer occurs in the collection";

    //result = recoverNSequenceForwardSequentially(file1, fileOutBwt, kmers.size());
    //assert (result ==1);

    //Free the memory
    for ( AlphabetSymbol j = 0 ; j < sizeAlpha; j++ )
    {
        delete [] tableOcc[j];
        tableOcc[j] = NULL;
    }
    delete [] tableOcc;

    delete[] alphaInverse;

    return 1;
}

//Computes the rank-inverse function for many sequences by using the vector and update posN with the number of symbols that it read.
//Computes the position of the i-th occurrences of the symbol toFindSymbolp[h] in the BWT.
//posN[h] is the number of occurrences of the symbol toFindSymbol[h] that I have to find in BWT corresponding to the i-th occurrence of the symbol in F.
int BCRexternalBWT::rankInverseManyByVector ( char const *file1, char const *fileOutBwt, SequenceNumber numKmersInput, uchar *toFindSymbols )
{
    uchar *buf = new uchar[SIZEBUFFER];

    //Timer timer;

    SequenceNumber j = 0;
    while ( j < numKmersInput )
    {
        //std::cerr << "===j= " << j << " vectTriple[j].pileN " << (int)vectTriple[j].pileN << " vectTriple[j].seqN " << vectTriple[j].seqN <<"\n";
        //We work into one BWT-partial at the time.
        AlphabetSymbol currentPile = vectTriple[j].pileN;
        //#ifdef DEBUG
        // std::cerr << "===Current BWT-partial= " << (int)currentPile << "\n";
        //#endif
        Filename newfilename( file1, "-B0", currentPile, "" );
        FILE *InFileBWT = fopen( newfilename, "rb" );
        if ( InFileBWT == NULL )
        {
            std::cerr << "rankInverseManyByVector: BWT file " << ( int )j << ": Error opening " << newfilename << std::endl;
            exit ( EXIT_FAILURE );
        }
        //contaAperturaFile++;

        SequenceNumber k = j;
        //For each tripla in the same current pile, we have to find the position of toRead-th occurrences of the symbol toFindSymbol
        while ( ( k < numKmersInput ) && ( vectTriple[k].pileN == currentPile ) )
        {
            uchar toFindSymbol = toFindSymbols[k];
            //std::cerr << "===k= " << k << " vectTriple[k].pileN " << (int)vectTriple[k].pileN << " vectTriple[k].seqN " << vectTriple[k].seqN << " toFindSymbol " << toFindSymbol <<"\n";
            if ( toFindSymbol != terminatorChar )
            {

                //Update each tripla, so posN is the number that we should read in a particular block into BWT-partial
                LetterNumber readChar = 0;
                LetterNumber toRead = vectTriple[k].posN;
                LetterNumber numBlock = 0;
                //Find the block
                //std::cerr << "toRead "<< toRead << " numBlock " << numBlock << " currentPile " << (int)currentPile << " toFindSymbol " << toFindSymbol<< "\n";
                //cerr << "UPDATE_POS...";
                //timer.timeNow();
                int result = update_Pos_Pile_Blocks( &toRead, &numBlock, currentPile, toFindSymbol );
                //cerr << "done." << timer << endl;
                checkIfEqual( result, 1 );
                readChar = numBlock * DIMBLOCK;
                //We have read readChar by using vectorOcc, now we have to read toRead symbols from the block numBlock
                //We move the file pointer in the position where the numBlock block starts.
                fseek ( InFileBWT, numBlock * DIMBLOCK, SEEK_SET );
                //Find the occurrences in the found block
                LetterNumber num = 0, num_read = 0;
                while ( ( !feof( InFileBWT ) ) && ( toRead > 0 ) )
                {
                    num_read = fread( buf, sizeof( uchar ), SIZEBUFFER, InFileBWT );
                    num = 0;
                    while ( ( num < num_read ) && ( toRead > 0 ) )
                    {
                        if ( buf[num] == toFindSymbol )
                            toRead--;
                        readChar++;       //it is the number of read symbols
                        num++;
                    }
                    if ( toRead < 0 )
                        std::cerr << "rankInverseManyByVector: position of the symbol not found" << "\n";
                }
                if ( toRead > 0 )
                {
                    std::cerr << "*Error rankInverseManyByVector: we should read " << toRead << " characters yet in " << newfilename << " file!\n";
                    exit ( EXIT_FAILURE );
                }
                //Update the value of posN
                vectTriple[k].posN = readChar;
            }
            k++;
        }
        j = k;
        fclose( InFileBWT );
    }

    delete [] buf;

    return 1;
}

//Computes the rank-inverse function and returns the number of symbols that it read.
//Computes the position of the i-th occurrences of the symbol toFindSymbol in the BWT.
//toRead is the number of occurrences of the symbol toFindSymbol that I have to find in BWT corresponding to the i-th occurrence of the symbol in F.
LetterNumber BCRexternalBWT::findRankInBWT ( char const *file1, char const *fileOutBwt, AlphabetSymbol currentPile, LetterNumber toRead, uchar toFindSymbol )
{
    assert( string( fileOutBwt ) == "" && "todo: remove this parameter if it is always null" );

    Filename newfilename( file1, "-B0", currentPile, "" );

    uchar *buf = new uchar[SIZEBUFFER];

    FILE *InFileBWT = fopen( newfilename, "rb" );
    if ( InFileBWT == NULL )
    {
        std::cerr << "findRankInBWT: could not open file " << newfilename << " !" << std::endl;
        exit ( EXIT_FAILURE );
    }

    LetterNumber num = 0, num_read = 0, readChar = 0;
    //#ifdef DEBUG
    // std::cerr << "***FindRankInBWT: we must to read " << toRead << " occurrences of the symbol " << toFindSymbol << "!\n";
    //#endif

    while ( ( !feof( InFileBWT ) ) && ( toRead > 0 ) )
    {
        num_read = fread( buf, sizeof( uchar ), SIZEBUFFER, InFileBWT );
        num = 0;
        while ( ( num < num_read ) && ( toRead > 0 ) )
        {
            if ( buf[num] == toFindSymbol )
                toRead--;
            readChar++;       //it is the number of read symbols
            num++;
        }
        if ( toRead < 0 )
            std::cerr << "findRankInBWT: position of the symbol not found" << "\n";
    }

    if ( toRead > 0 )
    {
        std::cerr << "Error findRankInBWT: we should read " << toRead << " characters yet in " << newfilename << " file!\n";
        exit ( EXIT_FAILURE );
    }
    fclose( InFileBWT );
    delete [] buf;
    return readChar;
}

//Computes the rank-inverse function and returns the number of symbols that it read.
//Computes the position of the i-th occurrences of the symbol toFindSymbol in the BWT.
//toRead is the number of occurrences of the symbol toFindSymbol that I have to find in BWT corresponding to the i-th occurrence of the symbol in F.
LetterNumber BCRexternalBWT::findRankInBWTbyVector ( char const *file1, char const *fileOutBwt, AlphabetSymbol currentPile, LetterNumber toRead, uchar toFindSymbol )
{
    assert( string( fileOutBwt ) == "" && "todo: remove this parameter if it is always null" );

    Filename newfilename( file1, "-B0", currentPile, "" );

    uchar *buf = new uchar[SIZEBUFFER];

    FILE *InFileBWT = fopen( newfilename, "rb" );
    if ( InFileBWT == NULL )
    {
        std::cerr << "findRankInBWTbyVector: could not open file " << newfilename << " !" << std::endl;
        exit ( EXIT_FAILURE );
    }

    LetterNumber readChar = 0;
    LetterNumber numBlock = 0;

    int result = update_Pos_Pile_Blocks( &toRead, &numBlock, currentPile, toFindSymbol );
    checkIfEqual( result, 1 );
    readChar = numBlock * DIMBLOCK;

    fseek ( InFileBWT, numBlock * DIMBLOCK, SEEK_SET );

    LetterNumber num = 0, num_read = 0;
    while ( ( !feof( InFileBWT ) ) && ( toRead > 0 ) )
    {
        num_read = fread( buf, sizeof( uchar ), SIZEBUFFER, InFileBWT );
        num = 0;
        while ( ( num < num_read ) && ( toRead > 0 ) )
        {
            if ( buf[num] == toFindSymbol )
                toRead--;
            readChar++;       //it is the number of read symbols
            num++;
        }
        if ( toRead < 0 )
            std::cerr << "findRankInBWTbyVector: position of the symbol not found" << "\n";
    }

    if ( toRead > 0 )
    {
        std::cerr << "*Error findRankInBWTbyVector: we should read " << toRead << " characters yet in " << newfilename << " file!\n";
        exit ( EXIT_FAILURE );
    }
    fclose( InFileBWT );
    delete [] buf;

    return readChar;
}

//Computes the rank function and returns the number of symbols that it read.
//The rank function computes the number char less than the symbol c from the starting position (startPos) in the BWT to the position pos (endPos) in the BWT.
//Here, we compute the number of occurrences of each symbol from from the starting position (startPos) in the BWT to the position pos (endPos) in the BWT.
//The startPos is the position of the File pointer InFileBWT, the endPos depends on toRead
//In the original definition of the rank, startPos corresponds to the position 1 and endPos corresponds to the previous symbol.
//Here, we work by using \sigma partial BWTs.
//toRead is the number of symbols that I have to read before to find the symbol in B corresponding to the symbol in F.
LetterNumber BCRexternalBWT::rankManySymbols( FILE &InFileBWT, LetterNumber *counters, LetterNumber toRead, uchar *foundSymbol )
{
    LetterNumber numchar, cont = 0; //cont is the number of symbols already read!
    uchar *buffer = new uchar[SIZEBUFFER];

    //it reads toRead symbols from the fp file (Partial BWT)
    while ( toRead > 0 )            //((numchar!=0) && (toRead > 0)) {
    {
        if ( toRead <= SIZEBUFFER )    //Read toRead characters
        {
            numchar = fread( buffer, sizeof( uchar ), toRead, &InFileBWT );
            // we should always read/write the same number of characters
            checkIfEqual( numchar, toRead );
            *foundSymbol = buffer[numchar - 1];   //The symbol of the sequence k.  It is the symbol in the last position in the partial BWT that we have read.
        }
        else     //Read sizebuffer characters
        {
            numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, &InFileBWT );
            // we should always read/write the same number of characters
            checkIfEqual( numchar, SIZEBUFFER );
        }

        //For each symbol in the buffer, it updates the number of occurrences into counters
        for ( LetterNumber r = 0; r < numchar; r++ )
            counters[alpha[( int )buffer[r]]]++;  //increment the number of letter symbol into counters


        cont   += numchar;  //number of read symbols
        toRead -= numchar;  //number of remaining symbols to read
        if ( ( numchar == 0 ) && ( toRead > 0 ) ) //it means that we have read 0 character, but there are still toRead characters to read
        {
            std::cerr << "rankManySymbols: read 0 character, but there are still " << toRead << " characters to read  " << std::endl;
            exit ( EXIT_FAILURE );
        }
    }
    delete [] buffer;

    return cont;
}

LetterNumber BCRexternalBWT::rankManySymbolsByVector( FILE &InFileBWT, LetterNumber *counters, LetterNumber toRead, uchar *foundSymbol, uchar *foundQual, FILE *InFileBWTQual )
{
    const LetterNumber offset = toRead;
    LetterNumber numchar, count = 0; //count is the number of symbols already read!
    static uchar *bufferBlock = new uchar[DIMBLOCK];

    //it reads toRead symbols from the fp file (Partial BWT)
    while ( toRead > 0 )            //((numchar!=0) && (toRead > 0)) {
    {
        if ( toRead <= DIMBLOCK )    //Read toRead characters
        {
            numchar = fread( bufferBlock, sizeof( uchar ), toRead, &InFileBWT );
            checkIfEqual( numchar, toRead ); // we should always read/write the same number of characters

            *foundSymbol = bufferBlock[numchar - 1];   //The symbol of the sequence k.  It is the symbol in the last position in the partial BWT that we have read.
        }
        else     //Read sizebuffer characters
        {
            std::cerr << "rankManySymbolsByVector: Error to read is" << toRead << std::endl;
            exit ( EXIT_FAILURE );
            //numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,&InFileBWT);
            //assert(numchar == SIZEBUFFER); // we should always read/write the same number of characters
            //aggiorna counters dalla tabella del vettori
        }

        //For each symbol in the buffer, it updates the number of occurrences into counters
        for ( LetterNumber r = 0; r < numchar; r++ )
            counters[alpha[( int )bufferBlock[r]]]++;  //increment the number of letter symbol into counters


        count  += numchar;  //number of read symbols
        toRead -= numchar;  //number of remaining symbols to read
        if ( ( numchar == 0 ) && ( toRead > 0 ) ) //it means that we have read 0 character, but there are still toRead characters to read
        {
            std::cerr << "rankManySymbolsByVector: read 0 character, but there are still " << toRead << " characters to read  " << std::endl;
            exit ( EXIT_FAILURE );
        }
    }
    //    delete [] bufferBlock;

    if ( foundQual )
    {
        if ( offset > 1 )
            fseek( InFileBWTQual, offset - 1, SEEK_CUR );
        numchar = fread( foundQual, sizeof( uchar ), 1, InFileBWTQual );
        checkIfEqual( numchar, 1 );

        size_t pos1 = ftell( &InFileBWT );
        size_t pos2 = ftell( InFileBWTQual );
        checkIfEqual( pos1, pos2 );
    }

    return count;
}


int BCRexternalBWT::computeNewPositionForBackSearch( char const *file1, char const *fileOutBwt, uchar symbol )
{
    //Last = C[c] + rank (c, Last)    --> (vectTriple[1].pileN, vectTriple[1].posN)

    //First = C[c] + rank (c, First - 1) + 1    --> (vectTriple[0].pileN, vectTriple[0].posN)
    //So we write:
    vectTriple[0].posN --;   //So we compute rank until position First - 1

    uchar foundSymbol = '\0';
    LetterNumber toRead = 0;
    LetterNumber *counters = new LetterNumber[sizeAlpha];  //it counts the number of each symbol into the i-Pile-BWT
    SequenceNumber j = 0;
    while ( j < 2 )
    {
        for ( AlphabetSymbol i = 0 ; i < sizeAlpha; i++ )
            counters[i] = 0;

        AlphabetSymbol currentPile = vectTriple[j].pileN;
        //if (verboseDecode == 1)
        // std::cerr << "\n===Current BWT-partial= " << (int)currentPile << "\n";
        Filename newfilename( file1, "-B0", currentPile, "" );

        FILE *InFileBWT = fopen( newfilename, "rb" );
        if ( InFileBWT == NULL )
        {
            std::cerr << "computeNewPositionForBackSearch: BWT file " << ( int )j << ": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        SequenceNumber k = j;
        LetterNumber cont = 0;   //number of the read symbols
        LetterNumber numberRead = 0;
        //uchar symbol;
        //SequenceLength lenCheck=0;
        while ( ( k < 2 ) && ( vectTriple[k].pileN == currentPile ) )
        {
            //The symbol for the sequences seqN in F[posN]  is the symbol
            //symbol = alphaInverse[vectTriple[k].pileN];
            //Now, I have to find the new symbol, it is in B[pileN] in position posN and so I can update pileN and posN

            //For any character (of differents sequences) in the same pile
            //symbol = '\0';
            //cont is the number of symbols already read!
            toRead = vectTriple[k].posN - cont;
            numberRead = rankManySymbols( *InFileBWT, counters, toRead, &foundSymbol );
            checkIfEqual( toRead, numberRead );
            cont += numberRead;

            //I have to update the value in vectTriple[k].posN, it must contain the position of the symbol in F
            //Symbol is
            //newSymb[vectTriple[k].seqN] = symbol;   //it is not useful here
            //PosN is
            vectTriple[k].posN = counters[alpha[( int )symbol]];
            for ( AlphabetSymbol g = 0 ; g < currentPile; g++ )  //I have to count in each pile g= 0... (currentPile-1)-pile
            {
                vectTriple[k].posN = vectTriple[k].posN + tableOcc[g][alpha[( int )symbol]];
            }
            //pileN is
            vectTriple[k].pileN = alpha[( int )symbol];
            k++;
        }
        fclose( InFileBWT );

        j = k;
    }


    delete [] counters;

    //First = c[c] + rank (c, First - 1) + 1
    vectTriple[0].posN ++;  //We must to sum 1 to first

    return 1;
}

int BCRexternalBWT::computeNewPositionForBackSearchByVector( char const *file1, char const *fileOutBwt, uchar symbol )
{
    //Last = C[c] + rank (c, Last)    --> (vectTriple[1].pileN, vectTriple[1].posN)

    //First = C[c] + rank (c, First - 1) + 1    --> (vectTriple[0].pileN, vectTriple[0].posN)
    //So we write:
    vectTriple[0].posN --;   //So we compute rank until position First - 1

    uchar foundSymbol = '\0';
    LetterNumber toRead = 0;
    LetterNumber *counters = new LetterNumber[sizeAlpha];  //it counts the number of each symbol into the i-Pile-BWT
    SequenceNumber j = 0;
    while ( j < 2 )
    {
        //The symbol for the sequences seqN in F[posN]  is the symbol
        //symbol = alphaInverse[vectTriple[k].pileN];
        //Now, I have to find the new symbol, it is in B[pileN] in position posN and so I can update pileN and posN
        for ( AlphabetSymbol i = 0 ; i < sizeAlpha; i++ )
            counters[i] = 0;

        AlphabetSymbol currentPile = vectTriple[j].pileN;
        //if (verboseDecode == 1)
        // std::cerr << "\n===Current BWT-partial= " << (int)currentPile << "(computeNewPositionForBackSearchByVector)\n";

        Filename newfilename( file1, "-B0", currentPile, "" );
        FILE *InFileBWT = fopen( newfilename, "rb" );
        if ( InFileBWT == NULL )
        {
            std::cerr << "computeNewPositionForBackSearchByVector: BWT file " << ( int )j << ": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        SequenceNumber k = j;
        LetterNumber cont = 0;   //number of the read symbols
        LetterNumber numberRead = 0;
        //uchar symbol;
        //SequenceLength lenCheck=0;
        LetterNumber numBlock = 0;
        while ( ( k < 2 ) && ( vectTriple[k].pileN == currentPile ) )
        {
            //cont is the number of symbols already read!
            //toRead = vectTriple[k].posN - cont;
            toRead = vectTriple[k].posN;
            for ( AlphabetSymbol i = 0 ; i < sizeAlpha; i++ )
                counters[i] = 0;

            if ( toRead > 0 )
            {
                //we need to know how many occurrences of each symbol there are up to the position toRead.
                //if ToRead > dimBlock, we can use vectorOcc in order to find the occorrences in the blocks precede the block where the position toRead is.
                //Before, we need to find the block where toRead position is.
                int result = findBlockToRead( counters, currentPile, &toRead, &numBlock );
                checkIfEqual( result, 1 );
            }

            if ( toRead <= DIMBLOCK )   //If toRead == DIMBLOCK, because I can need to known foundSymbol character
            {
                fseek ( InFileBWT, numBlock * DIMBLOCK, SEEK_SET );
                numberRead = rankManySymbolsByVector( *InFileBWT, counters, toRead, &foundSymbol );
                checkIfEqual( toRead, numberRead );
                cont += numberRead;
            }

            //I have to update the value in vectTriple[k].posN, it must contain the position of the symbol in F
            //Symbol is
            //newSymb[vectTriple[k].seqN] = symbol;   //it is not useful here
            //PosN is
            vectTriple[k].posN = counters[alpha[( int )symbol]];
            for ( AlphabetSymbol g = 0 ; g < currentPile; g++ )  //I have to count in each pile g= 0... (currentPile-1)-pile
            {
                vectTriple[k].posN = vectTriple[k].posN + tableOcc[g][alpha[( int )symbol]];
            }
            //pileN is
            vectTriple[k].pileN = alpha[( int )symbol];
            k++;
        }
        fclose( InFileBWT );

        j = k;
    }


    delete [] counters;

    //First = c[c] + rank (c, First - 1) + 1
    vectTriple[0].posN ++;  //We must to sum 1 to first

    return 1;
}

int BCRexternalBWT::findBlockToRead( LetterNumber *counters, AlphabetSymbol currentPile, LetterNumber *toRead, LetterNumber *numBlock )
{
    //Find the block numblock, where the position toRead is
    //numBlock = 0;
    *numBlock = ( LetterNumber )floor( ( long double )( ( *toRead - 1 ) / DIMBLOCK ) ) ; //The smallest integral value NOT less than x.
    //if (*numBlock >= numBlocksInPartialBWT[currentPile])
    //  std::cerr << "Error findBlockToRead: numBlock " << *numBlock << " and numBlocksInPartialBWT["<<(int)currentPile<<"]" << numBlocksInPartialBWT[currentPile] << "\n";
    //assert(*numBlock < numBlocksInPartialBWT[currentPile]);

    if ( *numBlock >= numBlocksInPartialBWT[currentPile] )
    {
        cerr << "Numblock size mismatch: " << *numBlock << " < "
             << numBlocksInPartialBWT[currentPile]
             << ". Aborting." << endl;
    }

    if ( *numBlock > 0 )
    {
        for ( AlphabetSymbol r = 0; r < sizeAlpha; r++ )
            counters[r] =  vectorOcc[currentPile][r][( *numBlock ) - 1]; //vectorOcc is indexed by 0, so we have numBlock-1
        *toRead = *toRead - ( *numBlock * DIMBLOCK ); //Number of symbols that we must read yet. it could be = DIMBLOCK
    }

    return 1;
}

int BCRexternalBWT::computeManyNewPositionForBackSearchByVector( char const *file1, char const *fileOutBwt, uchar *symbols, SequenceNumber nKmers )
{
    //Last = C[c] + rank (c, Last)    --> (vectTriple[1].pileN, vectTriple[1].posN)

    //First = C[c] + rank (c, First - 1) + 1    --> (vectTriple[0].pileN, vectTriple[0].posN)
    //So we write:
    for ( SequenceNumber i = 0; i < nKmers; i++ ) //For each kmer
        if ( LastVector[i].posN >= FirstVector[i].posN ) //if not, the kmer is not in the collection
            FirstVector[i].posN --;   //So we compute rank until position First - 1

    uchar foundSymbol = '\0';  //here, it is not useful
    LetterNumber toRead = 0;
    LetterNumber *counters = new LetterNumber[sizeAlpha];  //it counts the number of each symbol into the i-Pile-BWT
    SequenceNumber j = 0;
    while ( j < nKmers )
    {
        //The symbol for the sequences seqN in F[posN]  is the symbol
        //symbol = alphaInverse[vectTriple[k].pileN];
        //Now, I have to find the new symbol, it is in B[pileN] in position posN and so I can update pileN and posN
        for ( AlphabetSymbol i = 0 ; i < sizeAlpha; i++ )
            counters[i] = 0;

        AlphabetSymbol currentPile = FirstVector[j].pileN;
        //if (verboseDecode == 1)
        // std::cerr << "===Current BWT-partial= " << (int)currentPile << " (computeManyNewPositionForBackSearchByVector)\n";

        Filename newfilename( file1, "-B0", currentPile, "" );
        FILE *InFileBWT = fopen( newfilename, "rb" );
        if ( InFileBWT == NULL )
        {
            std::cerr << "computeManyNewPositionForBackSearchByVector: BWT file " << ( int )j << ": Error opening " << newfilename << std::endl;
            exit ( EXIT_FAILURE );
        }
        SequenceNumber k = j;
        //LetterNumber cont = 0;   //number of the read symbols
        LetterNumber numberRead = 0;
        //uchar symbol;
        //SequenceLength lenCheck=0;
        LetterNumber numBlock = 0;
        while ( ( k < nKmers ) && ( FirstVector[k].pileN == currentPile ) )
        {
            if ( FirstVector[k].posN <= LastVector[k].posN )
            {
                //FIRST
                //cont is the number of symbols already read!
                //toRead = vectTriple[k].posN - cont;
                toRead = FirstVector[k].posN;
                for ( AlphabetSymbol i = 0 ; i < sizeAlpha; i++ )
                    counters[i] = 0;

                if ( toRead > 0 )
                {
                    //we need to know how many occurrences of each symbol there are up to the position toRead.
                    //if ToRead > dimBlock, we can use vectorOcc in order to find the occurrences in the blocks that precede the block where the position toRead is.
                    //Before, we need to find the block where toRead position is.
                    int result = findBlockToRead( counters, currentPile, &toRead, &numBlock );
                    checkIfEqual( result, 1 );

                    if ( toRead <= DIMBLOCK )   //If toRead == DIMBLOCK, because I can need to known foundSymbol character
                    {
                        fseek ( InFileBWT, numBlock * DIMBLOCK, SEEK_SET );
                        numberRead = rankManySymbolsByVector( *InFileBWT, counters, toRead, &foundSymbol );
                        checkIfEqual( toRead, numberRead );
                        //cont += numberRead;
                    }
                }
                //I have to update the value in vectTriple[k].posN, it must contain the position of the symbol in F
                //Symbol is
                //newSymb[vectTriple[k].seqN] = symbol;   //it is not useful here
                //PosN is
                FirstVector[k].posN = counters[alpha[( int )symbols[FirstVector[k].seqN]]];
                for ( AlphabetSymbol g = 0 ; g < currentPile; g++ )  //I have to count in each pile g= 0... (currentPile-1)-pile
                {
                    FirstVector[k].posN = FirstVector[k].posN + tableOcc[g][alpha[( int )symbols[FirstVector[k].seqN]]];
                }
                //pileN is
                FirstVector[k].pileN = alpha[( int )symbols[FirstVector[k].seqN]];
                //First = c[c] + rank (c, First - 1) + 1
                FirstVector[k].posN ++;  //We must add 1 to first

                //LAST
                toRead = LastVector[k].posN;
                numBlock = 0;
                for ( AlphabetSymbol i = 0 ; i < sizeAlpha; i++ )
                    counters[i] = 0;

                if ( toRead > 0 )
                {
                    //we need to know how many occurrences of each symbol there are up to the position toRead.
                    //if ToRead > dimBlock, we can use vectorOcc in order to find the occorrences in the blocks precede the block where the position toRead is.
                    //Before, we need to find the block where toRead position is.
                    int result = findBlockToRead( counters, currentPile, &toRead, &numBlock );
                    checkIfEqual( result, 1 );

                    if ( toRead <= DIMBLOCK )   //If toRead == DIMBLOCK, because I can need to known foundSymbol character
                    {
                        fseek ( InFileBWT, numBlock * DIMBLOCK, SEEK_SET );
                        numberRead = rankManySymbolsByVector( *InFileBWT, counters, toRead, &foundSymbol );
                        checkIfEqual( toRead, numberRead );
                        //cont += numberRead;
                    }
                }
                //I have to update the value in vectTriple[k].posN, it must contain the position of the symbol in F
                //Symbol is
                //newSymb[vectTriple[k].seqN] = symbol;   //it is not useful here
                //PosN is
                LastVector[k].posN = counters[alpha[( int )symbols[FirstVector[k].seqN]]];
                for ( AlphabetSymbol g = 0 ; g < currentPile; g++ )  //I have to count in each pile g= 0... (currentPile-1)-pile
                {
                    LastVector[k].posN = LastVector[k].posN + tableOcc[g][alpha[( int )symbols[FirstVector[k].seqN]]];
                }
                //pileN is
                LastVector[k].pileN = alpha[( int )symbols[FirstVector[k].seqN]];
            }
            k++;
        }
        fclose( InFileBWT );

        j = k;
    }


    delete [] counters;

    return 1;
}

int BCRexternalBWT::backwardSearchManyBCR( char const *file1, char const *fileOutBwt, char const *fileOut, vector<string> kmers, SequenceLength lenKmer )
{
    assert( unbwtParams_ || searchParams_ );
    if ( ( unbwtParams_ && unbwtParams_->getValue( BeetlUnbwtParameters::PARAMETER_USE_VECTOR ) == BeetlUnbwtParameters::USE_VECTOR_ON )
         || searchParams_ )
    {
        std::cerr << "For the computation of the new position useful for BackSearch, it uses a sampling of the occurrences for each segment: " << DIMBLOCK << " size." << std::endl;
    }
    else
    {
        std::cerr << "backwardSearchManyBCR is only implemented by using the sampling." << std::endl;
        exit( EXIT_FAILURE );
    }

    //Initialization
    uchar *symbols = new uchar[kmers.size()];
    FirstVector.resize( kmers.size() );
    LastVector.resize( kmers.size() );
    for ( SequenceNumber i = 0; i < kmers.size(); i++ )
        std::cerr << kmers[i] << "\n";

    for ( SequenceNumber i = 0; i < kmers.size(); i++ ) //For each kmer
    {
        symbols[i] = kmers[i][lenKmer - 1];
        //FIRST
        FirstVector[i].seqN = i;  //It is not useful
        FirstVector[i].pileN = alpha[int( symbols[FirstVector[i].seqN] )];
        FirstVector[i].posN = 1;  //The first occurrence of symbol in F is in the first position in the pile Symbol

        //LAST
        LastVector[i].seqN = i;  //It is not useful
        LastVector[i].pileN = alpha[int( symbols[LastVector[i].seqN] )];
        //The last occurrence of the symbol prevSymbol in F is in the last position in the pile prevSymbol
        //It also corresponds to C[int(symbol) + 1]
        LastVector[i].posN = 0;
        for ( AlphabetSymbol mm = 0 ; mm < sizeAlpha; mm++ )
            LastVector[i].posN += tableOcc[LastVector[i].pileN][mm];
    }
    /*
    if (verboseDecode==1) {
            std::cerr << "Init triples: "  <<  "\n";
            std::cerr << "Symbols in positions " << lenKmer << "\n";
            for (SequenceNumber g = 0 ; g < kmers.size(); g++) {
                std::cerr << symbols[g]  << "\t";
            }
            std::cerr << std::endl;
            std::cerr << "Q  ";
            for (SequenceNumber g = 0 ; g < kmers.size(); g++) {
                std::cerr << (int)FirstVector[g].pileN << " " << (int)LastVector[g].pileN << "\t";
            }
            std::cerr << std::endl;
            std::cerr << "P  ";
            for (SequenceNumber g = 0 ; g < kmers.size(); g++) {
                    std::cerr << FirstVector[g].posN  << " " << LastVector[g].posN  << "\t";
            }
            std::cerr << std::endl;
            std::cerr << "N  ";
            for (SequenceNumber g = 0 ; g < kmers.size(); g++) {
                std::cerr << FirstVector[g].seqN  << " " << LastVector[g].seqN  << "\t" ;
            }
            std::cerr << std::endl;
    }
    */

    for ( SequenceLength posSymb = lenKmer - 1; posSymb > 0; posSymb-- ) //For each symbol of the kmer
    {

        for ( SequenceNumber i = 0; i < kmers.size(); i++ ) //For each kmer in accord to the order in the triples
            if ( LastVector[i].posN >= FirstVector[i].posN ) //if not, the kmer is not in the collection
                symbols[FirstVector[i].seqN] = kmers[FirstVector[i].seqN][posSymb - 1];

        quickSort( FirstVector );
        quickSort( LastVector );

        //For each symbol in the kmer we have to update First and Last
        int resultCompute = computeManyNewPositionForBackSearchByVector ( file1, fileOutBwt, symbols, kmers.size() );
        checkIfEqual( resultCompute, 1 );
        /*
        if (verboseDecode==1) {
            std::cerr << "After The computation of the new positions: "  <<  "\n";
            std::cerr << "Symbols in positions " << posSymb << "\n";
            for (SequenceNumber g = 0 ; g < kmers.size(); g++) {
                std::cerr << symbols[FirstVector[g].seqN]  << "\t\t";
            }
            std::cerr << std::endl;
            std::cerr << "Q  ";
            for (SequenceNumber g = 0 ; g < kmers.size(); g++) {
                std::cerr << (int)FirstVector[g].pileN << " " << (int)LastVector[g].pileN << "\t";
            }
            std::cerr << std::endl;
            std::cerr << "P  ";
            for (SequenceNumber g = 0 ; g < kmers.size(); g++) {
                    std::cerr << FirstVector[g].posN  << " " << LastVector[g].posN  << "\t";
            }
            std::cerr << std::endl;
            std::cerr << "N  ";
            for (SequenceNumber g = 0 ; g < kmers.size(); g++) {
                std::cerr << FirstVector[g].seqN  << " " << LastVector[g].seqN  << "\t" ;
            }
            std::cerr << std::endl;
        }
        */

    }

    delete[] symbols;

    return 1;
}

//Reconstruct 1 factor backwards by threading through the LF-mapping.
int BCRexternalBWT::backwardSearchBCR( char const *file1, char const *fileOutBwt, char const *fileOut, char const *kmer )
{
    // LetterNumber freq[256];  //contains the distribution of the symbols.
    // int resultInit = initializeUnbuildBCR(file1, fileOutBwt, freq);
    // assert (resultInit == 1);

    std::cerr << "Now: backward search\n";

    SequenceLength lenKmer = strlen( kmer );
    std::cerr << "kmer: " << kmer << " length " << lenKmer << "\n";

    LetterNumber posSymb = lenKmer - 1;
    uchar symbol = kmer[posSymb];
    vectTriple.resize( 2 );
    if ( verboseDecode == 1 )
        std::cerr << "\n>>>>>symbol is " <<  symbol << " in position " <<  posSymb + 1 <<  " of the pattern\n";

    //Initialize triplaFirst to find the first sequence
    //FIRST in position 0
    vectTriple[0].seqN = 0;  //It is not useful
    vectTriple[0].pileN = alpha[int( symbol )];
    vectTriple[0].posN = 1;  //The first occurrence of symbol in F is in the first position in the pile Symbol

    //Initialize triplaLast to find the last sequence
    //LAST in position 1
    vectTriple[1].seqN = 0;  //It is not useful
    vectTriple[1].pileN = alpha[int( symbol )];
    //The last occurrence of the symbol prevSymbol in F is in the last position in the pile prevSymbol
    //It also corresponds to C[int(symbol) + 1]
    vectTriple[1].posN = 0;
    for ( AlphabetSymbol j = 0 ; j < sizeAlpha; j++ )
        vectTriple[1].posN += tableOcc[vectTriple[1].pileN][j];

    if ( verboseDecode == 1 )
    {
        std::cerr << "Init triples: "  <<  "\n";
        std::cerr << "Q  ";
        for ( SequenceNumber g = 0 ; g < 2; g++ )
        {
            std::cerr << ( int )vectTriple[g].pileN << " ";
        }
        std::cerr << std::endl;
        std::cerr << "P  ";
        for ( SequenceNumber g = 0 ; g < 2; g++ )
        {
            std::cerr << vectTriple[g].posN  << " ";
        }
        std::cerr << std::endl;
        std::cerr << "N  ";
        for ( SequenceNumber g = 0 ; g < 2; g++ )
        {
            std::cerr << vectTriple[g].seqN  << " ";
        }
        std::cerr << std::endl;
    }

    //The new positions of symbol followed by kmer[posSymb] in F is computed by following function

    assert( unbwtParams_ || searchParams_ );
    if ( ( unbwtParams_ && unbwtParams_->getValue( BeetlUnbwtParameters::PARAMETER_USE_VECTOR ) == BeetlUnbwtParameters::USE_VECTOR_ON )
         || searchParams_ )
    {
        std::cerr << "For the computation of the new position useful for BackSearch, it uses a sampling of the occurrences for each segment: " << DIMBLOCK << " size." << std::endl;
    }
    else
    {
        std::cerr << "For the computation of the new position useful for BackSearch you don't use the vector of the occurrences. You read the file" << std::endl;
    }


    while ( ( ( vectTriple[0].pileN == vectTriple[1].pileN ) && ( vectTriple[0].posN < vectTriple[1].posN ) ) && ( posSymb >= 1 ) )
    {
        symbol = kmer[posSymb - 1];
        if ( verboseDecode == 1 )
            std::cerr << "\n>>>>>symbol is " <<  symbol << " in position " <<  posSymb <<  " of the pattern\n";

        //The new positions of symbol followed by kmer[posSymb] in F is computed by following function

        int resultCompute = 0;
        assert( unbwtParams_ || searchParams_ );
        if ( ( unbwtParams_ && unbwtParams_->getValue( BeetlUnbwtParameters::PARAMETER_USE_VECTOR ) == BeetlUnbwtParameters::USE_VECTOR_ON )
             || searchParams_ )
        {
            resultCompute = computeNewPositionForBackSearchByVector ( file1, fileOutBwt, symbol );
        }
        else
        {
            resultCompute = computeNewPositionForBackSearch ( file1, fileOutBwt, symbol );
        }
        checkIfEqual( resultCompute, 1 );

        if ( verboseDecode == 1 )
        {
            std::cerr << "New triples: "  <<  "\n";
            std::cerr << "Q  ";
            for ( SequenceNumber g = 0 ; g < 2; g++ )
            {
                std::cerr << ( int )vectTriple[g].pileN << " ";
            }
            std::cerr << std::endl;
            std::cerr << "P  ";
            for ( SequenceNumber g = 0 ; g < 2; g++ )
            {
                std::cerr << vectTriple[g].posN  << " ";
            }
            std::cerr << std::endl;
            std::cerr << "N  ";
            for ( SequenceNumber g = 0 ; g < 2; g++ )
            {
                std::cerr << vectTriple[g].seqN  << " ";
            }
            std::cerr << std::endl;
        }
        posSymb--;
        if ( verboseDecode == 1 )
            std::cerr << ">>>>>Next symbol in position " << posSymb << "\n";
    }

    return vectTriple[1].posN - vectTriple[0].posN + 1;
}

int BCRexternalBWT::computeVectorUnbuildBCR( char const *file1, char const *fileOutBwt, LetterNumber freq[] )
{
    assert( string( fileOutBwt ) == "" && "todo: remove this parameter if it is always null" );
    numBlocksInPartialBWT.resize( sizeAlpha );
    for ( AlphabetSymbol x = 0 ; x < sizeAlpha; x++ )
    {
        numBlocksInPartialBWT[x] = ( LetterNumber )ceil( ( long double )freq[alphaInverse[x]] / DIMBLOCK );
        if ( verboseDecode == 1 )
            std::cerr << "numBlocksInPartialBWT[ " << ( int )x << " ]= " << numBlocksInPartialBWT[x] << "\n";
    }


    // Start by allocating an array for array of arrays
    vectorOcc.resize( sizeAlpha );  //For each BWT-partial
    // Allocate an array for each element of the first array
    for ( AlphabetSymbol x = 0 ; x < sizeAlpha; x++ )       //For each block of BWT-partial
    {
        vectorOcc[x].resize( sizeAlpha );   //SumCumulative for each symbol and each block
        // Allocate an array of integers for each element of this symbol
        for ( AlphabetSymbol y = 0 ; y < sizeAlpha; y++ )       //For each block
            vectorOcc[x][y].resize( numBlocksInPartialBWT[x], 0 );
    }



    #pragma omp parallel for
    for ( AlphabetSymbol x = 0 ; x < sizeAlpha; x++ )   //For each BWT-partial
    {
        vector< uchar > bufBlock( DIMBLOCK );
        Filename newfilename( file1, "-B0", x, "" );

        FILE *InFileBWT = fopen( newfilename, "rb" );
        if ( InFileBWT == NULL )
        {
            std::cerr << "computeVectorUnbuildBCR: could not open file " << newfilename << " !" << std::endl;
            exit ( EXIT_FAILURE );
        }

        LetterNumber numBlock = 0;
        while ( !feof( InFileBWT ) && ( numBlock < numBlocksInPartialBWT[x] ) ) //Added check on numBlocks
        {
            SequenceLength num_read = fread( bufBlock.data(), sizeof( uchar ), DIMBLOCK, InFileBWT );
            for ( SequenceLength i = 0; i < num_read; i++ )
            {
                vectorOcc[x][alpha[( int )( bufBlock[i] )]][numBlock]++;
            }
            numBlock++;
        }
        fclose( InFileBWT );
        //Compute the sum cumulative for each BWT-partial
        for ( AlphabetSymbol z = 0 ; z < sizeAlpha; z++ )      //For each symbol z
            for ( LetterNumber y = 1; y < numBlocksInPartialBWT[x] ; y++ )      //For each block y>1 of partial-BWT x
                vectorOcc[x][z][y] = vectorOcc[x][z][y - 1] + vectorOcc[x][z][y]; //Sum the previous one: ie Blcok y and block y-1
    }

    /*
    #ifdef DEBUG
        for (AlphabetSymbol x = 0 ; x < sizeAlpha; x++) {
            std::cerr << "x = " << (int)x << " For the " << alphaInverse[x] << "-BWT-partial: the #symbols is " << freq[alphaInverse[x]] << " Number of block of the symbol " << numBlocksInPartialBWT[x] << "\n";
            for(AlphabetSymbol z = 0; z < sizeAlpha; ++z) {
                std::cerr << "Symbol: " << (int)z << ":\t";
                for(LetterNumber y = 0; y < numBlocksInPartialBWT[x]; ++y)
                    std::cerr << vectorOcc[x][z][y] << "\t";
                }
                std::cerr << "\n";
            }
            std::cerr << "\n";
        }
    #endif
    */
    return 1;
}


int BCRexternalBWT::initializeUnbuildBCR( char const *file1, char const *fileOutBwt, LetterNumber freq[] )
{
    assert( string( fileOutBwt ) == "" && "todo: remove this parameter if it is always null" );

    //We supposed that the symbols in the input file are the following
    //TODO
    for ( AlphabetSymbol i = 0; i < 255; i++ )
        freq[i] = 0;
    freq[int( terminatorChar )] = 1;
    freq[int( 'A' )] = 1;
    freq[int( 'C' )] = 1;
    freq[int( 'G' )] = 1;
    freq[int( 'N' )] = 1;
    freq[int( 'T' )] = 1;
    //GIOVANNA: ADDED THE SYMBOL Z IN THE ALPHABET, SO sizeAlpha = alphabetSize
#ifdef USE_EXTRA_CHARACTER_Z
    freq[int( 'Z' )] = 1;
#endif

    //Compute size of alphabet
    sizeAlpha = 0;
    for ( AlphabetSymbol i = 0; i < 255; i++ )
        if ( freq[i] > 0 )
            sizeAlpha++;

    //Compute alpha and alphaInverse
    alphaInverse = new AlphabetSymbol[sizeAlpha];
    AlphabetSymbol mmm = 0;
    for ( AlphabetSymbol i = 0; i < 255; i++ )
        if ( freq[i] > 0 )
        {
            alpha[i] = mmm;
            alphaInverse[mmm] = i;
            std::cerr << i << "\t" << freq[i] << "\t" << ( int )alpha[i] << "\t" << ( int )alphaInverse[mmm] << "\n";
            mmm++;
        }


    std::cerr << "sizeof(type size of alpha): " << sizeof( AlphabetSymbol ) << "\n";
    std::cerr << "sizeof(type of #sequences): " << sizeof( SequenceNumber ) << "\n";
    std::cerr << "sizeof(type of #characters): " << sizeof( LetterNumber ) << "\n";

    lengthTot = 0;  //Counts the number of symbols
    nText = 0;
    lengthRead = 0;
    lengthTot_plus_eof = 0;

    tableOcc = new LetterNumber*[sizeAlpha];
    for ( AlphabetSymbol j = 0 ; j < sizeAlpha; j++ )   //Counting for each pile: $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
    {
        tableOcc[j] = new LetterNumber[sizeAlpha];
    }
    for ( AlphabetSymbol j = 0 ; j < sizeAlpha; j++ )
        for ( AlphabetSymbol h = 0 ; h < sizeAlpha; h++ )
            tableOcc[j][h] = 0;


    LetterNumber lengthTotPlusEof = 0;
    #pragma omp parallel for reduction(+:lengthTotPlusEof)
    for ( AlphabetSymbol g = 0 ; g < sizeAlpha; g++ )
    {
        vector< uchar > buf( SIZEBUFFER );
        Filename newfilename( file1, "-B0", g, "" );

        FILE *InFileBWT = fopen( newfilename, "rb" );
        if ( InFileBWT == NULL )
        {
            std::cerr << "initializeUnbuildBCR: could not open file " << newfilename << " !" << std::endl;
            exit ( EXIT_FAILURE );
        }

        while ( !feof( InFileBWT ) )
        {
            SequenceLength num_read = fread( buf.data(), sizeof( uchar ), SIZEBUFFER, InFileBWT );
            for ( SequenceLength i = 0; i < num_read; i++ )
            {
                tableOcc[g][alpha[( int )( buf[i] )]]++;
            }
            lengthTotPlusEof += num_read;
        }
        fclose( InFileBWT );
    }
    lengthTot_plus_eof = lengthTotPlusEof;

    nText = 0;
    for ( AlphabetSymbol g = 0 ; g < sizeAlpha; g++ )
        nText += tableOcc[alpha[int( terminatorChar )]][g];
    lengthTot = lengthTot_plus_eof - nText;
    lengthRead = lengthTot / nText;
    std::cerr << "\nNumber of sequences: " << nText << "\n";
    std::cerr << "Length of each sequence: " << lengthRead << "\n\n";
    std::cerr << "Total length (without $): " << lengthTot << "\n";
    std::cerr << "Total length (with $): " << lengthTot_plus_eof << "\n";
    //if (verboseDecode == 1) {
    std::cerr << "TableOcc: "  << "\n";
    for ( AlphabetSymbol g = 0 ; g < sizeAlpha; g++ )
    {
        std::cerr << int( g )  << ":\t";
        for ( AlphabetSymbol j = 0 ; j < sizeAlpha; j++ )
            std::cerr << tableOcc[g][j]  << "\t";
        std::cerr << "\n";
    }
    //}

    for ( AlphabetSymbol j = 0 ; j < 255; j++ )
        freq[j] = 0;
    //Compute of the frequency of each symbol
    for ( AlphabetSymbol j = 0 ; j < sizeAlpha; j++ )
        for ( AlphabetSymbol h = 0 ; h < sizeAlpha; h++ )
            freq[( int )alphaInverse[j]] += tableOcc[j][h];

    return 1;
}

int BCRexternalBWT::unbuildBCR( char const *file1, char const *fileOutBwt, char const *fileOut, char const *fileOutput )
{
    bool processQualities = hasSuffix( fileOutput, ".fastq" );
    LetterNumber freq[256];  //contains the distribution of the symbols.
    int resultInit = initializeUnbuildBCR( file1, fileOutBwt, freq );
    checkIfEqual ( resultInit, 1 );

    if ( unbwtParams_->getValue( BeetlUnbwtParameters::PARAMETER_USE_VECTOR ) == BeetlUnbwtParameters::USE_VECTOR_ON )
    {
        resultInit = computeVectorUnbuildBCR( file1, fileOutBwt, freq );
        checkIfEqual( resultInit, 1 );
    }

    if ( unbwtParams_->getValue( BeetlUnbwtParameters::PARAMETER_DECODE_DIRECTION ) == BeetlUnbwtParameters::DECODE_DIRECTION_BACKWARD )
    {
        std::cerr << "Inverse BWT by Backward direction." << std::endl;
        decodeBCRmultipleReverse( file1, fileOutBwt, fileOut, processQualities );
        std::cerr << "The cyc files have been built. Building the sequences." << std::endl;
        TransposeFasta trasp;
        TmpFilename cycFilesInTmp( "cyc." );
        int res = trasp.convertFromCycFileToFastaOrFastq( cycFilesInTmp, fileOutput );
        checkIfEqual ( res, 1 );
        if ( deleteCycFile == 1 )
        {
            Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Removing auxiliary input files (cyc files)\n";

            // delete output files
            for ( SequenceLength i = 0; i < lengthRead; i++ )
            {
                TmpFilename filename1( fileOut, i, "" );
                if ( remove( filename1 ) != 0 )
                    std::cerr << filename1 << " BCRexternalBWT: Error deleting file" << std::endl;

                if ( processQualities )
                {
                    TmpFilename filename2( fileOut, i, ".qual" );
                    if ( remove( filename2 ) != 0 )
                        std::cerr << filename2 << " BCRexternalBWT: Error deleting file" << std::endl;
                }
            }
        }
    }
    else
    {
        std::cerr << "Inverse BWT by Forward direction."  << std::endl;
        decodeBCRnaiveForward( file1, fileOutBwt, fileOutput );
    }

    //Free the memory
    for ( AlphabetSymbol j = 0 ; j < sizeAlpha; j++ )
    {
        delete [] tableOcc[j];
        tableOcc[j] = NULL;
    }
    delete [] tableOcc;

    delete[] alphaInverse;

    return true;
}

int BCRexternalBWT::update_Pos_Pile( sortElement *tripla )
{
    //I have to find the position of toFindSymbol in corrected partial BWT
    //To find the pile, where the posN occurrences is, we use tableOcc.
    LetterNumber sumOccPileN = 0;
    AlphabetSymbol currentPile = 0;
    while ( ( sumOccPileN < tripla->posN ) && ( currentPile < sizeAlpha ) )
    {
        sumOccPileN += tableOcc[currentPile][tripla->pileN];
        currentPile++;
    }
    if ( sumOccPileN >= tripla->posN ) //it means that the pile, where toFindSymbol is, is currentPile-1 and the position
    {
        currentPile--;
        sumOccPileN = sumOccPileN - tableOcc[currentPile][tripla->pileN];
        tripla->posN = tripla->posN - sumOccPileN;
        tripla->pileN = currentPile;
    }
    else
        std::cerr << "update_Pos_Pile: symbol " << ( int )tripla->pileN << " not found: " << tripla->posN << "occurrence.\n";
    return 1;
}

int BCRexternalBWT::update_Pos_Pile_Blocks( LetterNumber *toRead, LetterNumber *numBlock, AlphabetSymbol currentPile, uchar toFindSymbol )
{
    //I have to find the position of toFindSymbol in corrected blocks in the partial BWT
    //To find the block in the pile, where the posN occurrences is, we use vectorOcc.
    /*
    //Linear scanning
    *numBlock = 0;
    while ((vectorOcc[currentPile][alpha[(int)(toFindSymbol)]][*numBlock] < *toRead) && (*numBlock <=numBlocksInPartialBWT[currentPile] ))  {  //Added checks on numBlocks
        (*numBlock)++;
    }
    assert (*numBlock <=numBlocksInPartialBWT[currentPile] );
    */
    //Binary search for scanning
    vector<LetterNumber>::iterator low;

    low = lower_bound ( vectorOcc[currentPile][alpha[( int )( toFindSymbol )]].begin(), vectorOcc[currentPile][alpha[( int )( toFindSymbol )]].end(), *toRead ); //          ^
    *numBlock = ( LetterNumber )( low - vectorOcc[currentPile][alpha[( int )( toFindSymbol )]].begin() );
    //assert (*numBlock <=numBlocksInPartialBWT[currentPile] );
    if ( *numBlock > numBlocksInPartialBWT[currentPile] )
    {
        cerr << "Numblock size mismatch: " << *numBlock << " < "
             << numBlocksInPartialBWT[currentPile]
             << ". Aborting." << endl;
    }

    if ( ( *numBlock ) > 0 )
    {
        *toRead = *toRead - vectorOcc[currentPile][alpha[( int )( toFindSymbol )]][*numBlock - 1]; //vectorOcc is indexed by 0
    }

    return 1;
}

vector <int> BCRexternalBWT::recoverNSequenceForward( char const *file1, char const *fileOutBwt, SequenceNumber numKmersInput )
{
    sortElement tripla;
    SequenceNumber numTotKmers = 0;
    for ( SequenceNumber g = 0 ; g < numKmersInput; g++ )
    {
        std::cerr << "Initialization for the k-mer of index " << LastVector[g].seqN << ". Number of occurrences " << LastVector[g].posN - FirstVector[g].posN + 1 << std::endl;
        //Initialize triple
        for ( LetterNumber  j = FirstVector[g].posN ; j <= LastVector[g].posN; j++ ) //For each position between first and last
        {
            tripla.seqN = numTotKmers;
            tripla.posN = j;
            tripla.pileN = FirstVector[g].pileN;
            vectTriple.push_back( tripla );
            numTotKmers++;
        }
    }
    std::cerr << "We want to compute the seqID of " << numTotKmers  << " sequences." << std::endl;

    quickSort( vectTriple );
    uchar *toFindSymbols = new uchar[numTotKmers];  //Symbol to find for each kmers
    SequenceNumber h;
    for ( h = 0 ; h < numTotKmers; h++ )
        toFindSymbols[h] = alphaInverse[vectTriple[h].pileN];  //The first time vectTriple[h].seqN == h

    h = 0 ;
    bool existDollars = false;  //We suppose that  there are not terminatorChar
    while ( ( h < numTotKmers ) && ( existDollars != true ) )
    {
        if ( toFindSymbols[h] != terminatorChar )
            existDollars = true;  //There are at least 1 symbol terminatorChar
        h++;
    }

    int result = 0;
    SequenceNumber countDollars = 0;
    while ( existDollars == true )   //If there is at least one terminatorChar symbol
    {
        //cerr << "another round existDollars" << endl;
        //Update the PileN where the number of occurrences is, for each h
        //posN is the number of occurrences (absolute value) in the entire BWT
        for ( h = 0 ; h < numTotKmers; h++ )
        {
            if ( toFindSymbols[vectTriple[h].seqN] != terminatorChar )   //if =terminatorChar, it means that we have already obtained the seqID
            {
                //cerr << "calling update_Pos_Pile" << endl;
                result = update_Pos_Pile( &vectTriple[vectTriple[h].seqN] ); //posN is the number of occurrences (relative value) in the pileN-BWT-partial
                checkIfEqual( result, 1 );
            }
        }
        //Compute the rank inverse and inserts the number of read symbols into posN. Update posN
        //cerr << "calling rankInverseManyByVector()";
        result = rankInverseManyByVector ( file1, fileOutBwt, numTotKmers, toFindSymbols );
        //cerr << "done." << endl;
        checkIfEqual( result, 1 );
        /*
        LetterNumber readChar = 0;
        for ( h = 0 ; h < numTotKmers; h++) {
            if (toFindSymbols[h] != terminatorChar) {
                readChar=findRankInBWTbyVector (file1, fileOutBwt, vectTriple[h].pileN, vectTriple[h].posN, toFindSymbols[h]);
                assert(readChar!=0);
                vectTriple[h].posN = readChar;
            }
        }
        */
        quickSort( vectTriple );
        //Update toFindSymbol and count the dollars
        countDollars = 0;
        for ( h = 0 ; h < numTotKmers; h++ )
        {
            if ( toFindSymbols[h] != terminatorChar )
            {
                toFindSymbols[h] = alphaInverse[vectTriple[h].pileN];
            }
            else
                countDollars++;
        }
        //std::cerr << "countDollars " << countDollars  << " ." << std::endl;

        if ( countDollars >= numTotKmers )
        {
            existDollars = false; //The end!
        }
    }

    vector <int> resultSeqId;
    resultSeqId.resize( numTotKmers );
    //The position is indexed by 1, the number of sequence by 0
    for ( h = 0 ; h < numTotKmers; h++ )
    {
        resultSeqId[vectTriple[h].seqN] = vectTriple[h].posN - 1;
    }
    vectTriple.clear();  //Erase all elements of vector.
    delete[] toFindSymbols;
    return resultSeqId;
}

int BCRexternalBWT::recoverNSequenceForwardSequentially( char const *file1, char const *fileOutBwt, SequenceNumber numKmersInput )
{
    //Compute the seqID sequentially

    //Now, you must find seqN for each position between vectTriple[0].posN and vectTriple[1].posN of the BWT-partial vectTriple[0].pileN=vectTriple[1].posN

    for ( SequenceNumber g = 0 ; g < numKmersInput; g++ )
    {
        //Recover the number of the sequence seqN of the kmer one symbol at time in reverse order
        uchar *sequence = new uchar[lengthRead + 2];
        sortElement tripla;
        std::cerr << "List of the seqID containing the k-mer with index" << g << ":" << std::endl;

        for ( LetterNumber  j = FirstVector[g].posN ; j <= LastVector[g].posN; j++ ) //For each position between first and last
        {
            for ( SequenceLength mmm = lengthRead + 2; mmm > 0; mmm-- )
                sequence[mmm - 1] = '\0';
            //Initialize tripla to find the sequence
            tripla.seqN = g;
            tripla.posN = j;
            tripla.pileN = FirstVector[g].pileN;
            SequenceLength lenSeq = 0;

            //#ifdef DEBUG
            // std::cerr << "Starting Tripla for the suffix: \tQ= " << (int)tripla.pileN << " P= " << tripla.posN  << " N= " << tripla.seqN  << std::endl;
            //#endif

            SequenceNumber numberOfSeq = recover1SequenceForward( file1, fileOutBwt, tripla, sequence, &lenSeq );
            //#ifdef DEBUG
            // std::cerr << " Computed suffix is " << sequence << "! It is long  " << lenSeq << ". It belongs to " << numberOfSeq << " sequence of the collection" << std::endl;
            //#endif

            std::cerr << "pos in the SA=\t" << j << "\t SeqId=\t" << numberOfSeq <<  std::endl;

        }
        delete [] sequence;
    }
    return 1;
}

//Reconstruct 1 sequence backwards by threading through the LF-mapping and reading the characters off of F column.
SequenceNumber BCRexternalBWT::recover1SequenceForward( char const *file1, char const *fileOutBwt, sortElement tripla, uchar *sequence, SequenceLength *lenCheck )
{
    //The toFindSymbol is into F column, it is in pileN-BWT in the position posN. So, it is the posN occurrences of alphaInverse[pileN] in F.
    //So, toFindSymbol is the alphaInverse[pileN]

    *lenCheck = 0;
    uchar toFindSymbol = alphaInverse[tripla.pileN];
    //LetterNumber rankFoundSymbol;

    // if (verboseDecode == 1) {
    //  std::cerr << "The symbol is: " << toFindSymbol << "\n";
    //  std::cerr << "\nI have to find the position of the " << tripla.posN << " " << toFindSymbol << " in the whole BWT\n";
    // }
    sequence[0] = toFindSymbol;
    //LetterNumber numcharWrite = fwrite (&toFindSymbol, sizeof(uchar), 1 , InfileOutDecode);
    //assert( numcharWrite == 1); // we should always read the same number of characters
    ( *lenCheck )++;
    while ( ( toFindSymbol != terminatorChar ) && ( *lenCheck <= lengthRead ) )
    {

        //posN is the number of occurrences (absolute value) in the entire BWT
        int result = update_Pos_Pile( &tripla ); //posN is the number of occurrences (relative value) in the pileN-BWT-partial
        checkIfEqual( result, 1 );

        //  if (verboseDecode == 1)
        //  std::cerr << "I have to find the position of the " << rankFoundSymbol << " occurrences of the symbol " <<  toFindSymbol << " in " << (int)tripla.pileN << " pile \n";
        //I have to read the pileN until I find rankFoundSymbol symbols. The found value is posN, i.e. the position of the next symbol

        LetterNumber readChar = 0;
        if ( unbwtParams_->getValue( BeetlUnbwtParameters::PARAMETER_USE_VECTOR ) != BeetlUnbwtParameters::USE_VECTOR_ON )
        {
            readChar = findRankInBWT ( file1, fileOutBwt, tripla.pileN, tripla.posN, toFindSymbol );
        }
        else
        {
            readChar = findRankInBWTbyVector ( file1, fileOutBwt, tripla.pileN, tripla.posN, toFindSymbol );
        }
        checkIfNotEqual( readChar, 0 );

        tripla.posN = readChar;
        //  if (verboseDecode == 1)
        //    std::cerr << "The occurrence " << rankFoundSymbol << " of the symbol " << toFindSymbol << " is in position " << tripla.posN << "\n\n";

        toFindSymbol = alphaInverse[tripla.pileN];
        sequence[*lenCheck] = toFindSymbol;
        //  if (verboseDecode == 1) {
        //    std::cerr << "The symbol is: " << toFindSymbol << "\n";
        //    std::cerr << "I have to find the position of the " << tripla.posN << "  " << toFindSymbol << " in the whole BWT\n";
        //  }
        ( *lenCheck )++;
    }
    //if (verboseDecode == 1)
    // std::cerr << lenCheck << " " << lengthRead << "\n";

    //if (verboseDecode==1) {
    // std::cerr << "***********Found the $-sign in First column \t";
    // std::cerr << "Q= " << (int)tripla.pileN << " P= " << tripla.posN  << " N= " << tripla.seqN  << std::endl;
    //}

    //The position is indexed by 1, the number of sequence by 0
    return tripla.posN - 1;
}

//Inverse BWT by Forward direction of nText sequences, one sequence at a time, in lexicographic order.
//Reconstruct the sequences one at a time in forward order
//file1 is the input file
//fileOutBwt is the suffix of the auxiliary files for the partial BWTs
//fileOutDecode is the output, that is the texts
int BCRexternalBWT::decodeBCRnaiveForward( char const *file1, char const *fileOutBwt, char const *fileOutDecode )
{
    LetterNumber numchar;

    Filename fileEndPos( file1, "-end-pos" );
    FILE *InFileEndPos;                  // input file of the end positions;
    InFileEndPos = fopen( fileEndPos, "rb" );
    if ( InFileEndPos == NULL )
    {
        std::cerr << "decodeBCRnaiveForward: could not open file " << fileEndPos << " !" << std::endl;
        exit ( EXIT_FAILURE );
    }

    FILE *InfileOutDecode = fopen( fileOutDecode, "wb" );
    if ( InfileOutDecode == NULL )
    {
        std::cerr << "decodeBCRnaiveForward: could not open file " << fileOutDecode << " !" << std::endl;
        exit ( EXIT_FAILURE );
    }
    SequenceNumber numText = 0;
    numchar = fread ( &numText, sizeof( SequenceNumber ), 1 , InFileEndPos );
    checkIfEqual ( numchar, 1 );
    checkIfEqual ( nText, numText ); // we should always read the same number of Texts of the bwt
    uint8_t subSequenceCount = 0;
    numchar = fread ( &subSequenceCount, sizeof( uint8_t ), 1 , InFileEndPos );
    checkIfEqual ( numchar, 1 );
    uint8_t hasRevComp = 0;
    numchar = fread ( &hasRevComp, sizeof( uint8_t ), 1 , InFileEndPos );
    checkIfEqual ( numchar, 1 );

    sortElement triple;
    std::cerr << "Recover the sequences of the collection in lexicographic order. A sequence at a time!" << std::endl;
    assert( unbwtParams_ || searchParams_ );
    if ( ( unbwtParams_ && unbwtParams_->getValue( BeetlUnbwtParameters::PARAMETER_USE_VECTOR ) == BeetlUnbwtParameters::USE_VECTOR_ON )
         || searchParams_ )
    {
        std::cerr << "It is using the sampling of the BWT. It requires more memory!" << std::endl;
        std::cerr << "In order to do this, it uses a sampling of the occurrences for each segment: " << DIMBLOCK << " size." << std::endl;
    }
    else
    {
        std::cerr << "It is not using the sampling of the BWT. It requires more time!" << std::endl;
    }
    for ( SequenceNumber i = 0; i < nText; i++ )
    {
        Logger::out() << "Decoding sequence " << i << endl;

        numchar = fread ( &triple.seqN, sizeof( SequenceNumber ), 1 , InFileEndPos );
        checkIfEqual( numchar, 1 ); // we should always read the same number of characters
        //        numchar = fread ( &triple.posN, sizeof( LetterNumber ), 1 , InFileEndPos ); //it is the relative position of the $ in the partial BWT
        //        checkIfEqual( numchar, 1 ); // we should always read the same number of characters
        //        numchar = fread ( &triple.pileN, sizeof( AlphabetSymbol ), 1 , InFileEndPos );
        //        checkIfEqual( numchar, 1 ); // we should always read the same number of characters
        assert( false && "todo: recalculate posN and pileN which are not stored in end-pos file anymore" );
        uint8_t subSequenceNum;
        numchar = fread ( &subSequenceNum, sizeof( uint8_t ), 1 , InFileEndPos );
        checkIfEqual( numchar, 1 ); // we should always read the same number of characters

        //  if (verboseDecode == 1)
        //   std::cerr << std::endl << "Starting Tripla: " << triple.seqN << " " << triple.posN << " " << (int)triple.pileN << "\n" << std::endl;

        uchar *sequence = new uchar[lengthRead + 2];
        for ( SequenceLength j = lengthRead + 2; j > 0; j-- )
            sequence[j - 1] = '\0';

        SequenceLength lenSeq = 0;
        SequenceNumber numberOfSeq = recover1SequenceForward( file1, fileOutBwt, triple, sequence, &lenSeq );
        checkIfEqual( numberOfSeq, triple.seqN );

        //  std::cerr << "The " << i+1 <<"-th/" << nText <<" computed sequence is " << sequence << "! It is long  " << lenSeq << ". It belongs to " << numberOfSeq << " sequence of the collection" << std::endl;
        if ( verboseDecode == 1 )
            cerr << numberOfSeq << "\t" << sequence << endl;

        SequenceLength numcharWrite = 0;
        numcharWrite = fwrite ( sequence, sizeof( uchar ), lenSeq , InfileOutDecode );
        checkIfEqual( numcharWrite , lenSeq ); // we should always read the same number of characters
        numcharWrite = fwrite ( "\n", sizeof( char ), 1 , InfileOutDecode );
        checkIfEqual( numcharWrite, 1 ); // we should always read the same number of characters
        delete [] sequence;
    }

    fclose( InFileEndPos );
    fclose( InfileOutDecode );
    return true;
}

//Multiple Decoding the sequences (Build reverse sequence)
//Reconstruct m sequences backwards by threading through the FL-mapping and reading the characters off of L.
//file1 is the input file
//fileOutBWT is the suffix of the filename of the partial BWTs
//fileOut is the prefix of the lengthRead-filename (transpose texts: cyc.i)
//Inverse BWT by Backward direction of nText sequences at the same time by lengthRead iterations.
int BCRexternalBWT::decodeBCRmultipleReverse( char const *file1, char const *fileOutBwt, char const *fileOut, bool processQualities )
{
    vectTriple.resize( nText );

    //As I want to compute the reverse sequences, I need the position of $ in F
    for ( SequenceNumber g = 0 ; g < nText; g++ )
    {
        vectTriple[g].pileN = alpha[int( terminatorChar )]; //So the 0-pile
        vectTriple[g].posN = g + 1;
        vectTriple[g].seqN = g;
    }

    if ( verboseDecode == 1 )
    {
        std::cerr << "The Initial triples of $ in first column are!" << std::endl;
        std::cerr << "Q  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << ( int )vectTriple[g].pileN << " ";
        }
        std::cerr << std::endl;
        std::cerr << "P  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].posN  << " ";
        }
        std::cerr << std::endl;
        std::cerr << "N  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].seqN  << " ";
        }
        std::cerr << std::endl;
    }

    uchar *newSymb = new uchar[nText];
    uchar *newQual = processQualities ? ( new uchar[nText] ) : NULL;

    //As we recover the symbol in reverse order, I store the first found symbol in cyc.(length-1) file
    //and the last found symbol in cyc.0 file
    assert( unbwtParams_ || searchParams_ );
    if ( ( unbwtParams_ && unbwtParams_->getValue( BeetlUnbwtParameters::PARAMETER_USE_VECTOR ) == BeetlUnbwtParameters::USE_VECTOR_ON )
         || searchParams_ )
    {
        std::cerr << "It is using the sampling of the BWT. It requires more memory!" << std::endl;
        std::cerr << "In order to do this, it uses a sampling of the occurrences for each segment: " << DIMBLOCK << " size." << std::endl;
    }
    else
    {
        std::cerr << "It is not using the sampling of the BWT. It requires more time!" << std::endl;
    }
    for ( SequenceLength m = lengthRead ; m > 0 ; m-- )
    {
        Logger::out() << "Decoding cycle " << m << endl;

        int resultNsymbol = -1;
        assert( unbwtParams_ || searchParams_ );
        if ( ( unbwtParams_ && unbwtParams_->getValue( BeetlUnbwtParameters::PARAMETER_USE_VECTOR ) == BeetlUnbwtParameters::USE_VECTOR_ON )
             || searchParams_ )
        {
            resultNsymbol = RecoverNsymbolsReverseByVector( file1, fileOutBwt, newSymb, newQual );
        }
        else
        {
            resultNsymbol = RecoverNsymbolsReverse ( file1, fileOutBwt, newSymb, newQual );
        }
        checkIfEqual ( resultNsymbol , 1 );

        TmpFilename filename( fileOut, m - 1, "" );
        FILE *InfileOutDecodeCyc = fopen( filename, "wb" );
        if ( InfileOutDecodeCyc == NULL )
        {
            std::cerr << "decodeBCRmultipleReverse: could not open file " << filename << " !" << std::endl;
            exit ( EXIT_FAILURE );
        }

        LetterNumber numcharWrite = fwrite ( newSymb, sizeof( uchar ), nText , InfileOutDecodeCyc );
        checkIfEqual( numcharWrite, nText ); // we should always read the same number of characters

        fclose( InfileOutDecodeCyc );

        if ( processQualities )
        {
            TmpFilename qualFilename( fileOut, m - 1, ".qual" );
            FILE *InfileOutDecodeCycQual = fopen( qualFilename, "wb" );
            if ( InfileOutDecodeCycQual == NULL )
            {
                std::cerr << "decodeBCRmultipleReverse: could not open file " << qualFilename << " !" << std::endl;
                exit ( EXIT_FAILURE );
            }
            LetterNumber numcharWriteQual = fwrite ( newQual, sizeof( uchar ), nText , InfileOutDecodeCycQual );
            checkIfEqual( numcharWriteQual, nText ); // we should always read the same number of characters
            fclose( InfileOutDecodeCycQual );
        }
    }
    delete [] newSymb;

    return true;
}

//It is used to reconstruct m sequences backwards by threading through the FL-mapping and reading the characters off of L.
int BCRexternalBWT::RecoverNsymbolsReverse( char const *file1, char const *fileOutBwt, uchar *newSymb, uchar *newQual )
{
    assert( string( fileOutBwt ) == "" && "todo: remove this parameter if it is always null" );

    if ( newQual != 0 )
    {
        cerr << "TODO: Quality decompression not implemented for RecoverNsymbolsReverse" << endl;
        exit( -1 );
    }

    LetterNumber toRead = 0;
    LetterNumber *counters = new LetterNumber[sizeAlpha];  //it counts the number of each symbol into the i-Pile-BWT
    SequenceNumber j = 0;
    while ( j < nText )
    {
        for ( AlphabetSymbol i = 0 ; i < sizeAlpha; i++ )
        {
            counters[i] = 0;
        }
        AlphabetSymbol currentPile = vectTriple[j].pileN;
        //if (verboseDecode == 1)
        // std::cerr << "===Current BWT-partial= " << (int)currentPile << "\n";

        Filename newfilename( file1, "-B0", currentPile, "" );
        FILE *InFileBWT = fopen( newfilename, "rb" );
        if ( InFileBWT == NULL )
        {
            std::cerr << "RecoverNsymbolsReverse: BWT file " << ( int )j << ": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        SequenceNumber k = j;
        LetterNumber cont = 0;   //number of the read symbols
        uchar foundSymbol;
        //SequenceLength lenCheck=0;
        LetterNumber numberRead = 0;
        while ( ( k < nText ) && ( vectTriple[k].pileN == currentPile ) )
        {
            if ( verboseDecode == 1 )
            {
                std::cerr << "Sequence number " << k << "\n";
                std::cerr << "j-1: Q[" << k << "]=" << ( int )vectTriple[k].pileN << " P[" << k << "]=" << ( LetterNumber )vectTriple[k].posN << " N[" << k << "]=" << ( SequenceNumber )vectTriple[k].seqN << "\n";
            }
            //The symbol for the sequences seqN in F[posN]  is the symbol
            //symbol = alphaInverse[vectTriple[k].pileN];
            //Now, I have to find the new symbol, it is in B[pileN] in position posN and so I can update pileN and posN

            //For any character (of differents sequences) in the same pile
            foundSymbol = '\0';
            //cont is the number of symbols already read!
            toRead = vectTriple[k].posN - cont;

            numberRead = rankManySymbols( *InFileBWT, counters, toRead, &foundSymbol );

            if ( verboseDecode == 1 )
            {
                std::cerr << "toRead " << toRead << " Found Symbol is " << foundSymbol << "\n";
            }
            checkIfEqual ( toRead, numberRead );
            cont += numberRead;

            //I have to update the value in vectTriple[k].posN, it must contain the position of the symbol in F
            //Symbol is
            if ( verboseDecode == 1 )
                std::cerr << "vectTriple[k].seqN = " << vectTriple[k].seqN << " Symbol = " << foundSymbol << "\n";
            newSymb[vectTriple[k].seqN] = foundSymbol;
            //PosN is
            vectTriple[k].posN = counters[alpha[( int )foundSymbol]];
            //if (verboseDecode == 1)
            // std::cerr << "\nCompute PosN\nInit New P["<< k <<"]= " << vectTriple[k].posN <<std::endl;
            for ( AlphabetSymbol g = 0 ; g < currentPile; g++ )  //I have to count in each pile g= 0... (currentPile-1)-pile
            {
                vectTriple[k].posN = vectTriple[k].posN + tableOcc[g][alpha[( int )foundSymbol]];
                //if (verboseDecode == 1) {
                // std::cerr << "g= " << (int)g << " symbol= " << (int)symbol << " alpha[symbol]= "<< (int)alpha[(int)symbol] <<std::endl;
                // std::cerr << "Add New posN[k]=" << vectTriple[k].posN << " tableOcc[g][alpha[(int)symbol]] " << tableOcc[g][alpha[(int)symbol]] <<std::endl;
                //}
            }
            //pileN is
            //std::cerr << "\nCompute Pile\n";
            vectTriple[k].pileN = alpha[( int )foundSymbol];
            if ( verboseDecode == 1 )
                std::cerr << "Result: j  : Q[q]=" << ( int )vectTriple[k].pileN << " P[q]=" << ( LetterNumber )vectTriple[k].posN <<  " N[q]=" << ( SequenceNumber )vectTriple[k].seqN << std::endl << std::endl;

            k++;
        }
        fclose( InFileBWT );
        j = k;
    }
    delete [] counters;

    if ( verboseDecode == 1 )
    {
        std::cerr << "NewSymbols " ;
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << newSymb[g] << " ";
        }
        std::cerr << std::endl;
        std::cerr << "Before Sorting" << std::endl;
        std::cerr << "Q  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << ( int )vectTriple[g].pileN << " ";
        }
        std::cerr << std::endl;
        std::cerr << "P  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].posN  << " ";
        }
        std::cerr << std::endl;
        std::cerr << "N  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].seqN  << " ";
        }
        std::cerr << std::endl;
    }

    quickSort( vectTriple );

    if ( verboseDecode == 1 )
    {
        std::cerr << "After Sorting" << std::endl;
        std::cerr << "Q  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << ( int )vectTriple[g].pileN << " ";
        }
        std::cerr << std::endl;
        std::cerr << "P  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].posN  << " ";
        }
        std::cerr << std::endl;
        std::cerr << "N  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].seqN  << " ";
        }
        std::cerr << std::endl;
    }

    return 1;
}

//It is used to reconstruct m sequences backwards by threading through the FL-mapping and reading the characters off of L.
int BCRexternalBWT::RecoverNsymbolsReverseByVector( char const *file1, char const *fileOutBwt, uchar *newSymb, uchar *newQual )
{
    assert( string( fileOutBwt ) == "" && "todo: remove this parameter if it is always null" );

    LetterNumber toRead = 0;
    LetterNumber *counters = new LetterNumber[sizeAlpha];  //it counts the number of each symbol into the i-Pile-BWT
    for ( AlphabetSymbol i = 0 ; i < sizeAlpha; i++ )
        counters[i] = 0;
    SequenceNumber j = 0;
    while ( j < nText )
    {
        for ( AlphabetSymbol i = 0 ; i < sizeAlpha; i++ )
            counters[i] = 0;

        AlphabetSymbol currentPile = vectTriple[j].pileN;
        //if (verboseDecode == 1)
        // std::cerr << "===Current BWT-partial= " << (int)currentPile << "\n";

        Filename newfilename( file1, "-B0", currentPile, "" );
        FILE *InFileBWT = fopen( newfilename, "rb" );
        if ( InFileBWT == NULL )
        {
            std::cerr << "RecoverNsymbolsReverseByVector: BWT file " << ( int )j << ": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        FILE *InFileBWTQual = 0;
        if ( newQual )
        {
            Filename qualFilename( file1, "-Q0", currentPile );
            InFileBWTQual = fopen( qualFilename, "rb" );
            if ( InFileBWT == NULL )
            {
                std::cerr << "RecoverNsymbolsReverseByVector: BWT Quality file " << ( int )j << ": Error opening " << std::endl;
                exit ( EXIT_FAILURE );
            }
        }
        LetterNumber currentReadPos = 0;
        LetterNumber nextReadPos = 0;
        SequenceNumber k = j;
        //LetterNumber cont = 0;   //number of the read symbols
        uchar foundSymbol, foundQual;
        //SequenceLength lenCheck=0;
        while ( ( k < nText ) && ( vectTriple[k].pileN == currentPile ) )
        {
            if ( verboseDecode == 1 )
            {
                std::cerr << "Sequence number " << k << "\n";
                std::cerr << "j-1: Q[" << k << "]=" << ( int )vectTriple[k].pileN << " P[" << k << "]=" << ( LetterNumber )vectTriple[k].posN << " N[" << k << "]=" << ( SequenceNumber )vectTriple[k].seqN << "\n";
            }
            //The symbol for the sequences seqN in F[posN]  is the symbol
            //symbol = alphaInverse[vectTriple[k].pileN];
            //Now, I have to find the new symbol, it is in B[pileN] in position posN and so I can update pileN and posN

            //For any character (of differents sequences) in the same pile
            foundSymbol = '\0';
            foundQual = '\0';
            //cont is the number of symbols already read!
            //toRead = vectTriple[k].posN - cont;
            nextReadPos = toRead = vectTriple[k].posN;
            // always true: assert( currentReadPos == ftell( InFileBWT ) );

            if ( toRead > currentReadPos && toRead - currentReadPos < DIMBLOCK )
            {
                // Next position is close enough to the last one, so that we can avoid to fseek
                toRead -= currentReadPos;
            }
            else
            {
                LetterNumber numBlock = 0;
                for ( AlphabetSymbol i = 0 ; i < sizeAlpha; i++ )
                    counters[i] = 0;
                //std::cerr << "toRead is " << toRead << "\n";
                if ( toRead > 0 )
                {
                    //we need to know how many occurrences of each symbol there are up to the position toRead.
                    //if ToRead > dimBlock, we can use vectorOcc in order to find the occorrences in the blocks precede the block where the position toRead is.
                    //Before, we need to find the block where toRead position is.
                    int result = findBlockToRead( counters, currentPile, &toRead, &numBlock );
                    checkIfEqual ( result , 1 );
                    //std::cerr << "numBlock: " << numBlock << " toRead " << toRead << "\n";
                }

                assert ( toRead <= DIMBLOCK );   //If toRead == DIMBLOCK, because I can need to known foundSymbol character
                //std::cerr << "Move file to the position " << numBlock*DIMBLOCK <<  "\n";
                fseek ( InFileBWT, numBlock * DIMBLOCK, SEEK_SET );
                if ( newQual )
                    fseek ( InFileBWTQual, numBlock * DIMBLOCK, SEEK_SET );
            }

            LetterNumber numberRead = rankManySymbolsByVector( *InFileBWT, counters, toRead, &foundSymbol, newQual ? &foundQual : NULL, InFileBWTQual );
            checkIfEqual ( toRead , numberRead );
            //std::cerr << "foundSymbol " << (int)foundSymbol <<  "\n";
            //cont += numberRead;
            /*
                std::cerr << "counters  after FirstVector:\t";
                for (AlphabetSymbol i = 0 ; i < sizeAlpha; i++)
                   std::cerr << " " << counters[i];
                std::cerr << "\n";
            */
            //numberRead = rankManySymbols(*InFileBWT, counters, toRead, &foundSymbol);

            //std::cerr << "toRead " << toRead << " Found Symbol is " << foundSymbol << "\n";
            //assert (toRead == numberRead);
            //cont += numberRead;

            //I have to update the value in vectTriple[k].posN, it must contain the position of the symbol in F
            //Symbol is
            if ( verboseDecode == 1 )
                std::cerr << "vectTriple[k].seqN = " << vectTriple[k].seqN << " Symbol = " << foundSymbol << "\n";

            newSymb[vectTriple[k].seqN] = foundSymbol;
            if ( newQual )
            {
                newQual[vectTriple[k].seqN] = foundQual;
            }

            //PosN is
            vectTriple[k].posN = counters[alpha[( int )foundSymbol]];

            //if (verboseDecode == 1)
            // std::cerr << "\nCompute PosN\nInit New P["<< k <<"]= " << vectTriple[k].posN <<std::endl;
            for ( AlphabetSymbol g = 0 ; g < currentPile; g++ )  //I have to count in each pile g= 0... (currentPile-1)-pile
            {
                vectTriple[k].posN = vectTriple[k].posN + tableOcc[g][alpha[( int )foundSymbol]];
                //if (verboseDecode == 1) {
                // std::cerr << "g= " << (int)g << " symbol= " << (int)symbol << " alpha[symbol]= "<< (int)alpha[(int)symbol] <<std::endl;
                // std::cerr << "Add New posN[k]=" << vectTriple[k].posN << " tableOcc[g][alpha[(int)symbol]] " << tableOcc[g][alpha[(int)symbol]] <<std::endl;
                //}
            }
            //pileN is
            //std::cerr << "\nCompute Pile\n";
            vectTriple[k].pileN = alpha[( int )foundSymbol];
            if ( verboseDecode == 1 )
                std::cerr << "Result: j  : Q[q]=" << ( int )vectTriple[k].pileN << " P[q]=" << ( LetterNumber )vectTriple[k].posN <<  " N[q]=" << ( SequenceNumber )vectTriple[k].seqN << std::endl << std::endl;

            currentReadPos = nextReadPos;
            k++;
        }
        fclose( InFileBWT );
        if ( InFileBWTQual )
            fclose( InFileBWTQual );
        j = k;
    }
    delete [] counters;

    if ( verboseDecode == 1 )
    {
        std::cerr << "NewSymbols " ;
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << ( char )newSymb[g] << " ";
        }
        std::cerr << std::endl;
        std::cerr << "Before Sorting" << std::endl;
        std::cerr << "Q  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << ( int )vectTriple[g].pileN << " ";
        }
        std::cerr << std::endl;
        std::cerr << "P  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].posN  << " ";
        }
        std::cerr << std::endl;
        std::cerr << "N  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].seqN  << " ";
        }
        std::cerr << std::endl;
    }

    quickSort( vectTriple );

    if ( verboseDecode == 1 )
    {
        std::cerr << "After Sorting" << std::endl;
        std::cerr << "Q  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << ( int )vectTriple[g].pileN << " ";
        }
        std::cerr << std::endl;
        std::cerr << "P  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].posN  << " ";
        }
        std::cerr << std::endl;
        std::cerr << "N  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].seqN  << " ";
        }
        std::cerr << std::endl;
    }

    return 1;
}


//It is used to reconstruct 1 sequences backwards by threading through the FL-mapping and reading the characters off of L.
int BCRexternalBWT::Recover1symbolReverse( char const *file1, char const *fileOutBwt, uchar *newSymbol, sortElement *tripla )
{
    assert( string( fileOutBwt ) == "" && "todo: remove this parameter if it is always null" );

    LetterNumber toRead = 0;
    LetterNumber *counters = new LetterNumber[sizeAlpha];  //it counts the number of each symbol into the i-Pile-BWT
    for ( AlphabetSymbol i = 0 ; i < sizeAlpha; i++ )
        counters[i] = 0;

    AlphabetSymbol currentPile = tripla->pileN;
    //if (verboseDecode == 1)
    // std::cerr << "===Current BWT-partial= " << (int)currentPile << "\n";

    Filename newfilename( file1, "-B0", currentPile, "" );
    FILE *InFileBWT = fopen( newfilename, "rb" );
    if ( InFileBWT == NULL )
    {
        std::cerr << "Recover1symbolReverse: BWT file " << ( int )currentPile << ": Error opening " << std::endl;
        exit ( EXIT_FAILURE );
    }
    LetterNumber cont = 0;   //number of the read symbols
    uchar foundSymbol;
    //if (verboseDecode == 1) {
    // std::cerr << "j-1: Q["<<(int)currentPile<<"]=" << (int)tripla->pileN << " P["<<(int)currentPile<<"]=" << (LetterNumber)tripla->posN << " N["<<(int)currentPile<<"]=" << (SequenceNumber)tripla->seqN << "\n";
    //}
    //The symbol for the sequences seqN in F[posN]  is the symbol
    //symbol = alphaInverse[vectTriple[k].pileN];
    //Now, I have to find the new symbol, it is in B[pileN] in position posN and so I can update pileN and posN

    //For any character (of differents sequences) in the same pile
    foundSymbol = '\0';
    //cont is the number of symbols already read!
    toRead = tripla->posN - cont;

    LetterNumber numberRead = 0;
    numberRead = rankManySymbols( *InFileBWT, counters, toRead, &foundSymbol );
    //std::cerr << "toRead " << toRead << "Found Symbol is " << foundSymbol << "\n";
    checkIfEqual ( toRead , numberRead );
    cont += numberRead;

    //I have to update the value in tripla.posN, it must contain the position of the symbol in F
    //Symbol is
    //if (verboseDecode == 1)
    // std::cerr << "tripla.seqN = " << tripla->seqN << " Symbol = " << foundSymbol << "\n";
    *newSymbol = foundSymbol;
    //PosN is
    tripla->posN = counters[alpha[( int )foundSymbol]];

    //if (verboseDecode == 1)
    // std::cerr << "\nCompute PosN\nInit New P= " << (LetterNumber)tripla->posN <<std::endl;
    for ( AlphabetSymbol g = 0 ; g < currentPile; g++ )  //I have to count in each pile g= 0... (currentPile-1)-pile
    {
        tripla->posN = tripla->posN + tableOcc[g][alpha[( int )foundSymbol]];
        //if (verboseDecode == 1) {
        // std::cerr << "g= " << (int)g << " symbol= " << (int)foundSymbol << " alpha[symbol]= "<< (int)alpha[(int)symbol] <<std::endl;
        // std::cerr << "Add New posN[k]=" << vectTriple[k].posN << " tableOcc[g][alpha[(int)foundSymbol]] " << tableOcc[g][alpha[(int)foundSymbol]] <<std::endl;
        //}
    }
    //pileN is
    //std::cerr << "\nCompute Pile\n";
    tripla->pileN = alpha[( int )foundSymbol];

    //if (verboseDecode == 1)
    // std::cerr << "Result: j  : Q[q]=" << (int)tripla->pileN << " P[q]=" << (LetterNumber)tripla->posN <<  " N[q]=" << (SequenceNumber)tripla->seqN << std::endl << std::endl;

    fclose( InFileBWT );

    delete [] counters;
    /*
    if (verboseDecode==1) {
        std::cerr << "NewSymbols " ;
        std::cerr << *newSymbol << " ";
        std::cerr << std::endl;
        std::cerr << "Q  ";
        std::cerr << (int)tripla->pileN << " ";
        std::cerr << std::endl;
        std::cerr << "P  ";
        std::cerr << tripla->posN  << " ";
        std::cerr << std::endl;
        std::cerr << "N  ";
        std::cerr << tripla->seqN  << " ";
        std::cerr << std::endl;
    }
    */
    return true;
}


BCRexternalBWT::~BCRexternalBWT()
{
}
