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

#include "BCRexternalBWT.hh"
#include "BWTCollection.hh"
#include "BwtReader.hh"
#include "BwtWriter.hh"
#include "Filename.hh"
#include "LetterCount.hh"
#include "Logger.hh"
#include "SeqReader.hh"
#include "Timer.hh"
#include "TransposeFasta.hh"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>

//#define USE_OPENMP
#ifdef USE_OPENMP
#include <omp.h>
#endif //ifdef USE_OPENMP

using namespace std;
using namespace BeetlBwtParameters;


//#ifdef COMPRESS_BWT
//typedef BwtReaderRunLength BwtReader;
//typedef BwtWriterRunLength BwtWriter;
//#else
//typedef BwtReaderASCII BwtReader;
//typedef BwtWriterASCII BwtWriter;
//#endif


#define SIZEBUFFER 1024 // TBD get rid of this
unsigned int debugCycle = 0;
unsigned int lastDefragCycle = 0;
Timer timer;

void dumpRamFiles()
{
    extern vector< vector<unsigned char> > ramFiles;

    for ( unsigned int i = 0; i < ramFiles.size(); ++i )
    {
        Filename filename( "ramdump", i );
        ofstream os( filename );
        for ( unsigned int j = 0; j < ramFiles[i].size(); j += 2 )
        {
            os << ( ( ( unsigned int )ramFiles[i][j] ) & 0xFF ) << ", " << ( ( unsigned int )( ramFiles[i][j + 1] ) & 0xFF ) << endl;
        }
    }
}

void debugRamFile( char *filenameIn, size_t n, char *filenameOut = "tmp.debug" )
{
    BwtReaderIncrementalRunLength *pReader;
    pReader = new BwtReaderIncrementalRunLength( filenameIn );
    assert( pReader != NULL );

    Logger::out( LOG_FOR_DEBUGGING ) << "Writing " << filenameOut << endl;
    BwtWriterASCII *pWriter = new BwtWriterASCII( filenameOut );
    for ( size_t i = 0; i < n; ++i )
    {
        pReader->readAndSend( *pWriter, 1 );
    }
    fflush( 0 );

    while ( 1 )
    {
        pReader->readAndSend( *pWriter, 1 );
    }

    delete pWriter;
    delete pReader;
}

void BCRexternalBWT::dumpRamFileToFile( const char *filenameIn, const char *filenameOut )
{
    BwtReaderBase *pReader = instantiateBwtReaderForIntermediateCycle( filenameIn );
    assert( pReader != NULL );

    BwtWriterBase *pWriter = instantiateBwtWriterForLastCycle( filenameOut );
    assert( pWriter != NULL );

    Logger::out( LOG_FOR_DEBUGGING ) << "Writing " << filenameOut << endl;
    while ( pReader->readAndSend( *pWriter, 1 ) ) {}
    Logger::out( LOG_FOR_DEBUGGING ) << " ... done" << endl;

    delete pReader;
    delete pWriter;
}

void BCRexternalBWT::ReadFilesForCycle( const char *prefix, const dataTypelenSeq cycle, const dataTypeNSeq nText, uchar *newSymb, const bool processQualities, uchar *newQual )
{
    Filename filename( prefix, cycle, "" );
    FILE *InFileInputText = fopen( filename, "rb" );
    if ( InFileInputText == NULL )
    {
        std::cerr << filename << " : Error opening " << std::endl;
        exit ( EXIT_FAILURE );
    }
    size_t num = fread( newSymb, sizeof( uchar ), nText, InFileInputText );
    checkIfEqual( num, nText ); // we should always read the same number of characters
    fclose( InFileInputText );
    if ( processQualities )
    {
        Filename qualFilename( prefix, "qual.", cycle, "" );
        FILE *InQualFileInputText = fopen( qualFilename, "rb" );
        if ( InQualFileInputText == NULL )
        {
            std::cerr << "buildBCR: " << qualFilename << " : Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        size_t numQual = fread( newQual, sizeof( uchar ), nText, InQualFileInputText );
        checkIfEqual( numQual, nText );
        fclose( InQualFileInputText );
    }
}

BwtReaderBase *BCRexternalBWT::instantiateBwtReaderForFirstCycle( const char *filenameIn )
{
    BwtReaderBase *pReader = new BwtReaderASCII( filenameIn );
    assert( pReader );
    return pReader;
}

BwtWriterBase *BCRexternalBWT::instantiateBwtWriterForFirstCycle( const char *filenameIn )
{
    BwtWriterBase *pWriter = new BwtWriterASCII( filenameIn );
    assert( pWriter );
    return pWriter;
}

BwtReaderBase *BCRexternalBWT::instantiateBwtReaderForIntermediateCycle( const char *filenameIn, bool allowDefrag )
{
    BwtReaderBase *pReader = NULL;
    int intermediateFormat = bwtParams_->getValue( BWT_OPTION_INTERMEDIATE_FORMAT );
    switch ( intermediateFormat )
    {
        case INTERMEDIATE_FORMAT_ASCII:
            pReader = new BwtReaderASCII( filenameIn );
            break;
        case INTERMEDIATE_FORMAT_RLE:
            pReader = new BwtReaderRunLength( filenameIn );
            break;
        case INTERMEDIATE_FORMAT_MULTIRLE:
            pReader = new BwtReaderIncrementalRunLength( filenameIn );
            break;
        case INTERMEDIATE_FORMAT_HUFFMAN:
            pReader = new BwtReaderHuffman( filenameIn );
            break;
        default:
            cerr << "Error in BCRexternalBWT::instantiateBwtReaderForIntermediateCycle: unknown intermediate format: " << intermediateFormat << endl;
            exit ( EXIT_FAILURE );
    }
    assert( pReader );

    if ( allowDefrag && intermediateFormat == INTERMEDIATE_FORMAT_MULTIRLE )
    {
#define RAM_FILES_DEFRAG_FREQUENCY 10
        if ( debugCycle > 1
             && ( debugCycle % RAM_FILES_DEFRAG_FREQUENCY ) == 0
             && ( debugCycle + RAM_FILES_DEFRAG_FREQUENCY / 2 ) < lengthRead )
        {
            BwtReaderIncrementalRunLength *p = static_cast<BwtReaderIncrementalRunLength *>( pReader );
            p->defragment();

            #pragma omp critical
            {
                cout << "After defrag, time now: " << timer.timeNow();
            }

            delete( pReader );
            pReader = new BwtReaderIncrementalRunLength( filenameIn );


#ifdef DUMP_EACH_CYCLE
            {
                Filename debugFilename( "", debugCycle, string( filenameIn ) + ".afterDefrag.debug" );
                dumpRamFileToFile( filenameIn, debugFilename );
            }
#endif //ifdef DUMP_EACH_CYCLE

            #pragma omp critical
            lastDefragCycle = debugCycle;
        }
    }

    return pReader;
}

BwtWriterBase *BCRexternalBWT::instantiateBwtWriterForIntermediateCycle( const char *filenameOut )
{
    BwtWriterBase *pWriter = NULL;
    int intermediateFormat = bwtParams_->getValue( BWT_OPTION_INTERMEDIATE_FORMAT );
    switch ( intermediateFormat )
    {
        case INTERMEDIATE_FORMAT_ASCII:
            pWriter = new BwtWriterASCII( filenameOut );
            break;
        case INTERMEDIATE_FORMAT_RLE:
            pWriter = new BwtWriterRunLength( filenameOut );
            break;
        case INTERMEDIATE_FORMAT_MULTIRLE:
            pWriter = new BwtWriterIncrementalRunLength( filenameOut );
            break;
        case INTERMEDIATE_FORMAT_HUFFMAN:
            pWriter = new BwtWriterHuffman( filenameOut );
            break;
        default:
            cerr << "Error in BCRexternalBWT::instantiateBwtWriterForIntermediateCycle: unknown intermediate format: " << intermediateFormat << endl;
            exit ( EXIT_FAILURE );
    }
    assert( pWriter );
    return pWriter;
}

BwtWriterBase *BCRexternalBWT::instantiateBwtWriterForLastCycle( const char *filenameOut )
{
    BwtWriterBase *pWriter = NULL;
    int outputFormat = bwtParams_->getValue( BWT_OPTION_OUTPUT_FORMAT );
    switch ( outputFormat )
    {
        case OUTPUT_FORMAT_ASCII:
            pWriter = new BwtWriterASCII( filenameOut );
            break;
        case OUTPUT_FORMAT_RLE:
            pWriter = new BwtWriterRunLength( filenameOut );
            break;
        case OUTPUT_FORMAT_HUFFMAN:
            pWriter = new BwtWriterHuffman( filenameOut );
            break;
        default:
            cerr << "Error in BCRexternalBWT::instantiateBwtWriterForLastCycle: unknown output format: " << outputFormat << endl;
            exit ( EXIT_FAILURE );
    }
    assert( pWriter );
    return pWriter;
}


int BCRexternalBWT::buildBCR( char const *file1, char const *fileOut, const BwtParameters *bwtParams )
{
#ifdef USE_OPENMP
    //    if ( bwtParams->getValue( BWT_OPTION_PARALLEL_PROCESSING ) != PARALLEL_PROCESSING_OFF )
    {
        // Use nested openmp parallelisation
        omp_set_nested( 1 );
    }
#endif //ifdef USE_OPENMP

    const bool permuteQualities = ( bwtParams_->getValue( BWT_OPTION_PROCESS_QUALITIES ) == PROCESS_QUALITIES_PERMUTE );
    const bool readQualities = permuteQualities;

    string cycFilesPrefix;
    TransposeFasta transp;
    if ( bwtParams->getValue( BWT_OPTION_INPUT_FORMAT ) == INPUT_FORMAT_CYC )
    {
        cycFilesPrefix = string( file1 );
        transp.inputCycFile( cycFilesPrefix );
    }
    else
    {
        cycFilesPrefix = string( fileOut );
        FILE *f = fopen( file1, "rb" );
        SeqReaderFile *pReader( SeqReaderFile::getReader( f ) );
        transp.init( pReader, readQualities );
        transp.convert( file1, cycFilesPrefix );
        delete pReader;
        fclose( f );
    }

    nText = transp.nSeq;
    lengthRead = transp.lengthRead;
    lengthTot = transp.lengthTexts;
    bool processQualities = transp.hasProcessedQualities();

    dataTypelenSeq currentIteration = 0; // iteration counter of symbols insertion
    dataTypelenSeq currentCycleFileNum; // file num processed in the current iteration (may not be equal to file num being read, in case of prefetch)
    int cycleFileNumIncrement;
    if ( bwtParams->getValue( BWT_OPTION_REVERSE ) == REVERSE_ON )
    {
        currentCycleFileNum = 0;
        cycleFileNumIncrement = 1;
    }
    else
    {
        currentCycleFileNum = lengthRead - 1;
        cycleFileNumIncrement = -1;
    }

    //#ifdef REPLACE_TABLEOCC
    sizeAlpha = 0;
    for ( dataTypedimAlpha i = 0; i < 255; ++i )
        if ( transp.freq[i] > 0 )
        {
            alpha[i] = sizeAlpha;
            sizeAlpha++;
        }
    lengthTot_plus_eof = lengthTot + nText;

    if ( Logger::isActive( LOG_FOR_DEBUGGING ) )
    {
        Logger::out( LOG_FOR_DEBUGGING ) << "We supposed that the symbols in the input file are:\n";
        for ( dataTypedimAlpha i = 0; i < 255; ++i )
            if ( transp.freq[i] > 0 )
                Logger::out( LOG_FOR_DEBUGGING ) << i << " " << transp.freq[i] << " " << ( int )alpha[i] << "\n";
        //#endif

        Logger::out( LOG_FOR_DEBUGGING ) << "sizeof(type size of alpha): " << sizeof( dataTypedimAlpha ) << "\n";
        Logger::out( LOG_FOR_DEBUGGING ) << "sizeof(type of #sequences): " << sizeof( dataTypeNSeq ) << "\n";
        Logger::out( LOG_FOR_DEBUGGING ) << "sizeof(type of #characters): " << sizeof( dataTypeNChar ) << "\n";

        Logger::out( LOG_FOR_DEBUGGING ) << "\nalphabetSize: " << ( int )alphabetSize << "\n";
        Logger::out( LOG_FOR_DEBUGGING ) << "Number of sequences: " << nText << "\n";
        Logger::out( LOG_FOR_DEBUGGING ) << "Length of each sequence: " << lengthRead << "\n\n";
        Logger::out( LOG_FOR_DEBUGGING ) << "Total length (without $): " << lengthTot << "\n";
        Logger::out( LOG_FOR_DEBUGGING ) << "Total length (with $): " << lengthTot_plus_eof << "\n";
    }

    dataTypeNChar numchar;

    Logger::out( LOG_SHOW_IF_VERBOSE ) << "Partial File name for input: " << cycFilesPrefix << " \n\n";

    uchar *newSymb = new uchar[nText];
    uchar *newQual = processQualities ? ( new uchar[nText] ) : NULL;
    uchar *nextSymb = new uchar[nText];
    uchar *nextQual = processQualities ? ( new uchar[nText] ) : NULL;
    vectTriple.resize( nText );

#ifdef REPLACE_TABLEOCC
    tableOcc = new dataTypeNChar*[sizeAlpha];
    for ( dataTypedimAlpha j = 0 ; j < sizeAlpha; j++ )   //Counting for each pile: $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
    {
        tableOcc[j] = new dataTypeNChar[sizeAlpha];
    }
    for ( dataTypedimAlpha j = 0 ; j < sizeAlpha; j++ )
        for ( dataTypedimAlpha h = 0 ; h < sizeAlpha; h++ )
            tableOcc[j][h] = 0;
#endif
    tableOcc_.clear();

    Logger::out( LOG_SHOW_IF_VERBOSE ) << "\nFirst symbols: " << "Iteration " << 0 << " - symbols in (zero-based) position " << currentCycleFileNum << "\n";
    Logger::out( LOG_ALWAYS_SHOW ) << "Starting iteration " << currentIteration << ", time now: " << timer.timeNow();
    Logger::out( LOG_ALWAYS_SHOW ) << "Starting iteration " << currentIteration << ", usage: " << timer << endl;

    ReadFilesForCycle( cycFilesPrefix.c_str(), currentCycleFileNum, nText, newSymb, processQualities, newQual );
    InsertFirstsymbols( newSymb, newQual );
    // Update iteration counters
    ++currentIteration;
    currentCycleFileNum += cycleFileNumIncrement;
    debugCycle = currentIteration;

    if ( lengthRead >= 2 )
    {
        Logger::out( LOG_SHOW_IF_VERBOSE ) << "Reading next cycle files, time now: " << timer.timeNow();
        Logger::out( LOG_SHOW_IF_VERBOSE ) << "Reading next cycle files, usage: " << timer << endl;
        ReadFilesForCycle( cycFilesPrefix.c_str(), currentCycleFileNum, nText, newSymb, processQualities, newQual );
    }

    while ( currentIteration <= lengthRead - 2 )
        //    for ( dataTypelenSeq t = lengthRead - 2 ; t > 0; t-- ) //dataTypelenSeq is unsigned
    {
        pauseBetweenCyclesIfNeeded ();

        Logger::out( LOG_SHOW_IF_VERBOSE ) << "Iteration " << ( int ) currentIteration << " - symbols in position " << ( int ) currentCycleFileNum << std::endl;
        Logger::out( LOG_ALWAYS_SHOW ) << "Starting iteration " << currentIteration << ", time now: " << timer.timeNow();
        Logger::out( LOG_ALWAYS_SHOW ) << "Starting iteration " << currentIteration << ", usage: " << timer << endl;


        //To insert the symbol from position m-3 to position 1
        //The last inserted symbol is in position i+1 (or it is newSymb[j]),
        //the next symbol (to insert) is in position i


        #pragma omp parallel sections
        {
            #pragma omp section
            {
                InsertNsymbols( newSymb, currentIteration, newQual );
            }
            #pragma omp section
            {
                ReadFilesForCycle( cycFilesPrefix.c_str(), currentCycleFileNum + cycleFileNumIncrement, nText, nextSymb, processQualities, nextQual );
                Logger::out( LOG_SHOW_IF_VERBOSE ) << "Reading input file, time now: " << timer.timeNow();
            }
        }

        // swap data pointers
        uchar *tmp;
        tmp = newSymb;
        newSymb = nextSymb;
        nextSymb = tmp;
        if ( processQualities )
        {
            tmp = newQual;
            newQual = nextQual;
            nextQual = tmp;
        }

        // Update iteration counters
        ++currentIteration;
        currentCycleFileNum += cycleFileNumIncrement;
        debugCycle = currentIteration;
    }

    pauseBetweenCyclesIfNeeded ();

    //The last inserted symbol was in position 1 (or it is newSymb[j]),
    //the next symbol (to insert) is in position 0
    Logger::out( LOG_SHOW_IF_VERBOSE ) << "Iteration " << ( int ) currentIteration << " - symbols in position " << ( int ) currentCycleFileNum << std::endl;
    Logger::out( LOG_ALWAYS_SHOW ) << "Starting iteration " << currentIteration << ", time now: " << timer.timeNow();
    Logger::out( LOG_ALWAYS_SHOW ) << "Starting iteration " << currentIteration << ", usage: " << timer << endl;
    assert( currentIteration == lengthRead - 1 );
    assert( currentCycleFileNum == 0 || currentCycleFileNum == lengthRead - 1 ); // depending on the --reverse flag
    InsertNsymbols( newSymb, currentIteration, newQual );
    // Update iteration counters
    ++currentIteration;
    debugCycle = currentIteration;

    pauseBetweenCyclesIfNeeded ();

    //The last inserted symbol was in position 0 (or it is newSymb[j]),
    //the next symbol (to insert) is in position m-1, that is, I have to inserted the symbols $
    Logger::out( LOG_SHOW_IF_VERBOSE ) << "Final iteration " << ( int ) currentIteration << " - Inserting $=" << ( int )TERMINATE_CHAR << "=" << TERMINATE_CHAR << " symbols" << std::endl;
    Logger::out( LOG_ALWAYS_SHOW ) << "Starting iteration " << currentIteration << ", time now: " << timer.timeNow();
    Logger::out( LOG_ALWAYS_SHOW ) << "Starting iteration " << currentIteration << ", usage: " << timer << endl;
    assert( currentIteration == lengthRead );
    for ( dataTypeNSeq j = 0 ; j < nText; j++ )
    {
        newSymb[j] = '$';
        if ( readQualities )
            newQual[j] = 0;
    }
    InsertNsymbols( newSymb, currentIteration, newQual );

    Logger::out( LOG_ALWAYS_SHOW ) << "Final iteration complete, time now: " << timer.timeNow();
    Logger::out( LOG_ALWAYS_SHOW ) << "Final iteration complete, usage: " << timer << endl;


    if ( bwtParams->getValue( BWT_OPTION_GENERATE_ENDPOSFILE ) || BUILD_SA )
    {
        Logger::out( LOG_SHOW_IF_VERBOSE ) << "Storing 'outFileEndPos', time now: " << timer.timeNow();
        Logger::out( LOG_SHOW_IF_VERBOSE ) << "Storing 'outFileEndPos', usage: " << timer << endl;

        Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "Stores the 'end positions' of the $!" << std::endl;
        const char *fileEndPos = "outFileEndPos.bwt";
        static FILE *OutFileEndPos;                  // output file of the end positions;
        OutFileEndPos = fopen( fileEndPos, "wb" );
        if ( OutFileEndPos == NULL )
        {
            std::cerr << "Error opening \"" << fileEndPos << "\" file" << std::endl;
            exit ( EXIT_FAILURE );
        }

        /*
        //Each symbol newSymb[seqN[i]] has been in position posN[i] into the pile pileN[i]
        //We have to store the absolute positions in the entire BWT
        //So we need to use the tableOcc_.
        //The symbol $ of the sequence i is in the position endPos[SeqN[i]]
        for (dataTypeNSeq i = 0; i < nText; i++) {
        for (dataTypedimAlpha r = 0; r < vectTriple[i].pileN; r++) {
        for (dataTypedimAlpha t = 0; t < alphabetSize; t++) {
        vectTriple[i].posN += tableOcc_[r][t];
        }
        }
        }

        std::cerr << "Positions of the EOF into BWT" << std::endl;
        for (dataTypeNSeq i = 0; i < nText; i++) {
        std::cerr << posN[i] << " ";
        }
        std::cerr << std::endl;
        */

        numchar = fwrite ( &nText, sizeof( dataTypeNChar ), 1 , OutFileEndPos );
        assert( numchar == 1 ); // we should always read the same number of characters

        for ( dataTypeNSeq i = 0; i < nText; i++ )
        {
            if ( verboseEncode == 1 )
                std::cerr << "Triple: " << vectTriple[i].seqN << " " << vectTriple[i].posN << " " << ( int )vectTriple[i].pileN << std::endl;
            numchar = fwrite ( &vectTriple[i].seqN, sizeof( dataTypeNSeq ), 1 , OutFileEndPos );
            assert( numchar == 1 ); // we should always read the same number of characters
            numchar = fwrite ( &vectTriple[i].posN, sizeof( dataTypeNChar ), 1 , OutFileEndPos ); //here vectTriple[i].posN is the relative position of $ in the partial BWT
            assert( numchar == 1 ); // we should always read the same number of characters
            numchar = fwrite ( &vectTriple[i].pileN, sizeof( dataTypedimAlpha ), 1 , OutFileEndPos );
            assert( numchar == 1 ); // we should always read the same number of characters
        }

        fclose( OutFileEndPos );
        Logger::out( LOG_ALWAYS_SHOW ) << "'end positions' stored!" << std::endl;
    }

    /* We shouldn't need this anymore as long as we never send the last cycle to ram
      if (bwtParams->getValue( BWT_OPTION_INTERMEDIATE_STORAGE_MEDIUM ) == INTERMEDIATE_STORAGE_MEDIUM_RAM)
      {
    #pragma omp parallel for
          for (dataTypedimAlpha i = 1; i < alphabetSize; i++)
          {
              char inputFilename[1000];
              char outputFilename[1000];
              Filename inputFilename( i );
              Filename outputFilename( "final.", i );
              dumpRamFileToFile( outputCompression_, inputFilename, compressionRunLength, outputFilename);
          }
          cout << "Disk output complete, usage: " << timer << endl;
      }
    */

    // to delete those:
    delete [] newSymb;
    // vectTriple.~vector<sortElement>();
    /*
      delete [] seqN;
      delete [] pileN;
      delete [] posN;
    */
    //  std::cerr << std::endl;
    //std::cerr << "The input text is long " << lengthTot << std::endl;

    dataTypeNChar numCharInTable = 0;
    for ( dataTypedimAlpha r = 0; r < alphabetSize; r++ )
    {
        for ( dataTypedimAlpha t = 0; t < alphabetSize; t++ )
        {
            numCharInTable += tableOcc_[r].count_[t];
        }
    }
    Logger::out( LOG_FOR_DEBUGGING ) << "In tableOcc_, there are " << numCharInTable << " letters" << std::endl;

#ifdef REPLACING_TABLEOCC
    for ( dataTypedimAlpha j = 0 ; j < sizeAlpha; j++ )
    {
        delete [] tableOcc[j];
        tableOcc[j] = NULL;
    }
    delete [] tableOcc;
#endif

    if ( bwtParams_->getValue( BWT_OPTION_GENERATE_LCP ) == GENERATE_LCP_ON
         && bwtParams_->getValue( BWT_OPTION_OUTPUT_FORMAT ) != OUTPUT_FORMAT_ASCII )
    {
        // LCP source code only generate ASCII BWT, so we convert it here if needed
        Logger::out( LOG_SHOW_IF_VERBOSE ) << "Converting ASCII BWT to " << bwtParams_->getValueAsString( BWT_OPTION_OUTPUT_FORMAT ) << std::endl;
        for ( dataTypedimAlpha g = 1 ; g < alphabetSize; g++ )
        {
            Filename filename1( "", g, "" );
            Filename filename2( "new_", g, "" );

            BwtReaderBase *pReader = new BwtReaderASCII( filename1 );
            BwtWriterBase *pWriter = instantiateBwtWriterForLastCycle( filename2 );

            while ( pReader->readAndSend( *pWriter, 1000000000 ) > 0 ) {}

            delete pReader;
            delete pWriter;

            if ( remove( filename1 ) != 0 )
                std::cerr << "Error deleting file " << filename1 << std::endl;
            else if ( rename( filename2, filename1 ) )
                std::cerr << "Error renaming " << filename2 << " to " << filename1 << std::endl;
        }
    }

    return permuteQualities ? 2 : 1;
} // ~buildBCR

void BCRexternalBWT::InsertFirstsymbols( uchar const *newSymb, uchar const *newSymbQual )
{
    for ( dataTypeNSeq j = 0 ; j < nText; j++ )
    {
        vectTriple[j] = sortElement( 0, j + 1, j );
    }
    if ( verboseEncode == 1 )
    {
        std::cerr << "First step" << std::endl;
        std::cerr << "Q  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            std::cerr << ( int )vectTriple[g].pileN << " ";
        }
        std::cerr << std::endl;
        std::cerr << "P  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].posN  << " ";
        }
        std::cerr << std::endl;
        std::cerr << "N  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].seqN  << " ";
        }
        std::cerr << std::endl;

        if ( bwtParams_->getValue( BWT_OPTION_GENERATE_LCP ) == GENERATE_LCP_ON )
        {
            std::cerr << "C  ";             //LCP current
            for ( dataTypeNSeq g = 0 ; g < nText; g++ )
            {
                std::cerr << ( int )vectTriple[g].getLcpCurN()  << " ";
            }
            std::cerr << std::endl;
            std::cerr << "S  ";                     //LCP successive
            for ( dataTypeNSeq g = 0 ; g < nText; g++ )
            {
                std::cerr << ( int )vectTriple[g].getLcpSucN()  << " ";
            }
            std::cerr << std::endl;
        }
    }


    static FILE *OutFileBWT;                  // output file BWT;

    OutFileBWT = fopen( "0", "wb" );
    if ( OutFileBWT == NULL )
    {
        std::cerr << "BWT file $: Error opening " << std::endl;
        exit ( EXIT_FAILURE );
    }
    for ( dataTypeNSeq j = 0 ; j < nText; j++ )
    {
        tableOcc_[0].count_[whichPile[( int )newSymb[j]]]++;     //counting the number of occurrences in BWT of the $-pile
        /*
        if (newSymb[j] < 50) {
            std::cerr << (int)newSymb[j] <<" " << std::endl;
        }
        else {
            std::cerr << newSymb[j] <<" " << std::endl;
        }
        */
    }
    //Store newSymb into $-pile BWT
    dataTypeNChar num = fwrite ( newSymb, sizeof( uchar ), nText , OutFileBWT );
    checkIfEqual( num, nText ); // we should always read the same number of characters

    //if (num != nText)
    // std::cerr << "the written characters is not equal to number of the texts" << num << " and "<< nText <<"\n";
    fclose( OutFileBWT );

    const bool permuteQualities = ( bwtParams_->getValue( BWT_OPTION_PROCESS_QUALITIES ) == PROCESS_QUALITIES_PERMUTE );
    if ( permuteQualities )
    {
        assert ( newSymbQual );
        Filename filenameQualOut( "0.qual" );
        FILE *OutFileBWTQual = fopen( filenameQualOut, "wb" );
        if ( OutFileBWTQual == NULL )
        {
            std::cerr << "BWT file $: Error opening " << filenameQualOut << std::endl;
            exit ( EXIT_FAILURE );
        }
        dataTypeNChar numQual = fwrite ( newSymbQual, sizeof( uchar ), nText, OutFileBWTQual );
        checkIfEqual( numQual, nText );
        fclose( OutFileBWTQual );
    }


    //Creates one file for each letter in the alphabet. From 1 to alphabetSize-1
    //GIOVANNA: In the For, it is need the ''='' symbol. The maximal value must be alphabetSize-1
    for ( dataTypedimAlpha i = 1; i <= alphabetSize - 1; i++ )
    {
        Filename filenameOut( i );
        OutFileBWT = fopen( filenameOut, "wb" );
        if ( OutFileBWT == NULL )
        {
            std::cerr << "BWT file " << ( int )i << " : Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        fclose( OutFileBWT );

        if ( permuteQualities )
        {
            Filename filenameQualOut( filenameOut.str() + ".qual" );
            FILE *OutFileBWTQual = fopen( filenameQualOut, "wb" );
            if ( OutFileBWTQual == NULL )
            {
                std::cerr << "BWT file $: Error opening " << filenameQualOut << std::endl;
                exit ( EXIT_FAILURE );
            }
            fclose( OutFileBWTQual );
        }
    }

    if ( bwtParams_->getValue( BWT_OPTION_GENERATE_LCP ) == GENERATE_LCP_ON )
    {
        dataTypelenSeq *vectLCP = new dataTypelenSeq[nText];
        for ( dataTypeNSeq j = 0 ; j < nText; j++ )
        {
            vectLCP[j] = 0;
        }
        FILE *OutFileLCP;                  // output and input file LCP;
        Filename filenameOut( "0.lcp" );
        OutFileLCP = fopen( filenameOut, "wb" );
        if ( OutFileLCP == NULL )
        {
            std::cerr << "LCP file $: " << filenameOut << " Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }

        num = fwrite ( vectLCP, sizeof( dataTypelenSeq ), nText , OutFileLCP );
        checkIfEqual( num , nText ); // we should always read the same number of integers
        //vectLCP.clear();
        fclose( OutFileLCP );
        //Creates one file for each letter in the alphabet for LCP. From 1 to alphabetSize-1
        for ( dataTypedimAlpha i = 1; i < alphabetSize; i++ )
        {
            Filename filenameOut( "", i, ".lcp" );
            OutFileLCP = fopen( filenameOut, "wb" );
            if ( OutFileLCP == NULL )
            {
                std::cerr << "LCP file: " << filenameOut << " : Error opening " << std::endl;
                exit ( EXIT_FAILURE );
            }
        }
        fclose( OutFileLCP );
        delete [] vectLCP;
    }

    //Do we want compute the extended suffix array (position and number of sequence)?
    if ( BUILD_SA == 1 )  //To store the SA
    {
        Filename filenameOut( "sa_0" );
        static FILE  *OutFileSA = fopen( filenameOut, "wb" );    // output file SA;
        if ( OutFileSA == NULL )
        {
            std::cerr << "SA file: Error opening: " << filenameOut << " (SA file $)" << std::endl;
            exit ( EXIT_FAILURE );
        }

        ElementType *newEle = new ElementType[nText];
        for ( dataTypeNSeq j = 0 ; j < nText; j++ )
        {
            //newEle[j].sa=(iterationNum + 1) % (lengthRead + 1);
            newEle[j].sa = lengthRead;
            newEle[j].numSeq = j;
            //std::cerr << "(" << (int)newEle[j].sa << ", " << newEle[j].numSeq << ")\n";
        }
        //Store into $-pile SA
        dataTypeNChar num = fwrite ( newEle, sizeof( ElementType ), nText , OutFileSA );
        if ( num != nText )
            std::cerr << "Error: The written characters is not equal to number of the texts in SA" << num << " and " << nText << "\n";
        assert( num == nText );

        fclose( OutFileSA );

        //Creates one file for each letter in the alphabet. From 1 to alphabetSize-1
        for ( dataTypedimAlpha i = 1; i < alphabetSize; i++ )
        {
            Filename filenameOut( "sa_", i );

            OutFileSA = fopen( filenameOut, "wb" );
            if ( OutFileBWT == NULL )
            {
                std::cerr << "SA file " << ( int )i << " : Error opening " << std::endl;
                exit ( EXIT_FAILURE );
            }
            fclose( OutFileSA );
        }
    }
}

void BCRexternalBWT::InsertNsymbols( uchar const *newSymb, dataTypelenSeq iterationNum, uchar const *newQual )
{
    FILE *InFileBWT;                  // output and input file BWT;
    dataTypeNChar numchar = 0;

    // We first calculate at which index each pile starts
    vector<dataTypeNSeq> pileStarts( alphabetSize );
    pileStarts[0] = 0;
    dataTypeNSeq index = 0;
    for ( int pile = 1; pile < alphabetSize; ++pile )
    {
        while ( index < nText && vectTriple[index].pileN < pile )
            ++index;
        pileStarts[pile] = index;
    }
    /*
      for (unsigned int j = 0; j < alphabetSize-1; ++j)
      {
      clog << " => pile " << j << ": " << pileStarts[j] << "-" << pileStarts[j+1] <<endl;
      }
    */
    vector< vector< vector< sortElement > > > parallelVectTriplePerNewPile( alphabetSize, vector< vector< sortElement > >( alphabetSize ) );

    int parallelPile;
    //    for (int parallelPile = alphabetSize-2; parallelPile >= 0; --parallelPile)
    #pragma omp parallel for
    for ( parallelPile = 0; parallelPile < alphabetSize - 1; ++parallelPile )
    {
        InsertNsymbols_parallelPile( newSymb, iterationNum, newQual, parallelPile, pileStarts[parallelPile], pileStarts[parallelPile + 1], parallelVectTriplePerNewPile[parallelPile] );
    }


    if ( verboseEncode == 1 )
    {
        std::cerr << "The segments before inserting are:\n";
        uchar *buffer = new uchar[SIZEBUFFER];
        dataTypedimAlpha mmm = 0;
        while ( mmm < alphabetSize )
        {
            Filename filenameIn( mmm );
            //printf("===currentPile= %d\n",mmm);
            InFileBWT = fopen( filenameIn, "r" );
            for ( dataTypeNSeq g = 0 ; g < SIZEBUFFER; g++ )
                buffer[g] = '\0';
            numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, InFileBWT );
            std::cerr << "B[" << mmm << "]:\t";
            if ( numchar == 0 )
                std::cerr  << "empty\n";
            else
                std::cerr  << buffer << "\n";
            fclose( InFileBWT );
            mmm++;
        }
        delete [] buffer;
    } // ~if verboseEncode

    if ( verboseEncode == 1 )
    {
        if ( bwtParams_->getValue( BWT_OPTION_GENERATE_LCP ) == GENERATE_LCP_ON )
        {
            FILE *InFileLCP;                  // output and input file LCP;
            dataTypelenSeq *bufferLCP = new dataTypelenSeq[SIZEBUFFER];
            dataTypedimAlpha mmm = 0;
            while ( mmm < alphabetSize )
            {
                Filename filenameInLCP( "", mmm, ".lcp" );
                InFileLCP = fopen( filenameInLCP, "rb" );
                for ( dataTypeNChar g = 0 ; g < SIZEBUFFER; g++ )
                    bufferLCP[g] = 0;
                numchar = fread( bufferLCP, sizeof( dataTypelenSeq ), SIZEBUFFER, InFileLCP );
                std::cerr << "L[" << ( int )mmm << "]:\t";
                if ( numchar == 0 )
                    std::cerr  << "empty";
                else
                    for ( dataTypeNSeq g = 0 ; g < numchar; g++ )
                        std::cerr  << ( int )bufferLCP[g] << " ";
                std::cerr  << "\n";
                fclose( InFileLCP );
                mmm++;
            }
            delete [] bufferLCP;
        }

        std::cerr << "NewSymbols " ;
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            std::cerr << newSymb[g] << " ";
        }
        std::cerr << std::endl;
        std::cerr << "Before Sorting" << std::endl;
        std::cerr << "Q  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            std::cerr << ( int )vectTriple[g].pileN << " ";
        }
        std::cerr << std::endl;
        std::cerr << "P  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].posN  << " ";
        }
        std::cerr << std::endl;
        std::cerr << "N  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].seqN  << " ";
        }
        std::cerr << std::endl;
    } // ~if verboseEncode

#ifdef XXX
    delete [] counters;
#endif

    Logger::out( LOG_SHOW_IF_VERBOSE ) << "Before quicksort, time now: " << timer.timeNow();
    Logger::out( LOG_SHOW_IF_VERBOSE ) << "Before quicksort, usage: " << timer << endl;

    //    quickSort(vectTriple);
    vectTriple.clear();
    for ( uint newPile = 0; newPile < alphabetSize - 1; ++newPile )
    {
        for ( uint prevPile = 0; prevPile < alphabetSize - 1; ++prevPile )
        {
            if ( !parallelVectTriplePerNewPile[prevPile][newPile].empty() )
            {
                vectTriple.insert( vectTriple.end(), parallelVectTriplePerNewPile[prevPile][newPile].begin(), parallelVectTriplePerNewPile[prevPile][newPile].end() );
            }
        }
    }

    if ( Logger::isActive( LOG_FOR_DEBUGGING ) )
    {
        Logger::out( LOG_FOR_DEBUGGING ) << "After Sorting" << std::endl;
        Logger::out( LOG_FOR_DEBUGGING ) << "U  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            Logger::out( LOG_FOR_DEBUGGING ) << newSymb[g] << " ";
        }
        Logger::out( LOG_FOR_DEBUGGING ) << std::endl;
        Logger::out( LOG_FOR_DEBUGGING ) << "Q  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            Logger::out( LOG_FOR_DEBUGGING ) << ( int )vectTriple[g].pileN << " ";
        }
        Logger::out( LOG_FOR_DEBUGGING ) << std::endl;
        Logger::out( LOG_FOR_DEBUGGING ) << "P  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            Logger::out( LOG_FOR_DEBUGGING ) << vectTriple[g].posN  << " ";
        }
        Logger::out( LOG_FOR_DEBUGGING ) << std::endl;
        Logger::out( LOG_FOR_DEBUGGING ) << "N  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            Logger::out( LOG_FOR_DEBUGGING ) << vectTriple[g].seqN  << " ";
        }
        Logger::out( LOG_FOR_DEBUGGING ) << std::endl;

        if ( bwtParams_->getValue( BWT_OPTION_GENERATE_LCP ) == GENERATE_LCP_ON )
        {
            Logger::out( LOG_FOR_DEBUGGING ) << "C  ";
            for ( dataTypeNSeq g = 0 ; g < nText; g++ )
            {
                Logger::out( LOG_FOR_DEBUGGING ) << ( int )vectTriple[g].getLcpCurN()  << " ";
            }
            Logger::out( LOG_FOR_DEBUGGING ) << std::endl;
            Logger::out( LOG_FOR_DEBUGGING ) << "S  ";
            for ( dataTypeNSeq g = 0 ; g < nText; g++ )
            {
                Logger::out( LOG_FOR_DEBUGGING ) << ( int )vectTriple[g].getLcpSucN()  << " ";
            }
            Logger::out( LOG_FOR_DEBUGGING ) << std::endl;
        }
    } // ~if verboseEncode

    Logger::out( LOG_SHOW_IF_VERBOSE ) << "After quicksort, time now: " << timer.timeNow();
    Logger::out( LOG_SHOW_IF_VERBOSE ) << "After quicksort, usage: " << timer << endl;

    if ( bwtParams_->getValue( BWT_OPTION_GENERATE_LCP ) == GENERATE_LCP_OFF )
    {
        storeBWT( newSymb, newQual );
    }
    else
    {
        storeBWTandLCP( newSymb );
    }

    //std::cerr << "End storing BWT" << std::endl;

    //Do we want to compute the generalized suffix array (position and number of sequence)?
    if ( BUILD_SA == 1 )
    {
        storeSA( iterationNum );
    }

    //  delete pReader;
}

void BCRexternalBWT::InsertNsymbols_parallelPile( uchar const *newSymb, dataTypelenSeq iterationNum, uchar const *newQual, unsigned int parallelPile, dataTypeNSeq startIndex, dataTypeNSeq endIndex, vector< vector< sortElement > > &newVectTriplePerNewPile )
{
    if ( Logger::isActive( LOG_FOR_DEBUGGING ) )
    {
        #pragma omp critical
        {
            Logger::out( LOG_FOR_DEBUGGING ) << "InsertNsymbols_parallelPile: pile=" << parallelPile << " from " << startIndex << " to " << endIndex << endl;
        }
    }
#ifdef DEBUG
    ulong zz( 0 ); // TEMP
#endif

    //<<<<<<< BuildBCR.cpp
    //they are not first symbols
    //std::cerr << "Compute new posN" << std::endl;

    LetterCount counters;

    dataTypeNChar toRead = 0;
    //Find the positions of the new symbols
    dataTypeNSeq j = startIndex;

    BwtReaderBase *pReader( NULL );
    sortElement newVectTripleItem;

    counters.clear();

    if ( parallelPile == 6 )
    {
        assert( startIndex == endIndex );
        return;
    }
    if ( parallelPile == 0 && startIndex == endIndex )
    {
        return;
    }
    //    clog << "---------- " << (int)vectTriple[j].pileN << "|" << parallelPile << endl;
    dataTypedimAlpha currentPile = parallelPile; //vectTriple[j].pileN;
    assert ( currentPile < alphabetSize - 1 );
    Filename filename( "", currentPile, "" );
    //printf("===Current BWT-partial= %d\n",currentPile);

    //#define DUMP_EACH_CYCLE
#ifdef DUMP_EACH_CYCLE
    if ( bwtParams_->getValue( BWT_OPTION_INTERMEDIATE_STORAGE_MEDIUM ) == INTERMEDIATE_STORAGE_MEDIUM_RAM )
    {
        Filename debugFilename( "", debugCycle, "." + filename.str() + ".debug" );
        dumpRamFileToFile( filename, debugFilename );
    }
#endif //ifdef DUMP_EACH_CYCLE

    if ( iterationNum == 1 )
        pReader = instantiateBwtReaderForFirstCycle( filename );
    else
        pReader = instantiateBwtReaderForIntermediateCycle( filename );
    assert( pReader != NULL );

    if ( j < endIndex )
    {
        //    BwtReader reader(filename.str().c_str());

        dataTypeNSeq k = j;
        //For each pile, we have a different counter of characters
        for ( dataTypedimAlpha i = 0 ; i < alphabetSize; i++ )
            counters.count_[i] = 0;
        dataTypeNChar cont = 0;   //number of the read symbols
        uchar foundSymbol;
        dataTypeNChar numberRead = 0;
        while ( ( k < endIndex ) && ( vectTriple[k].pileN == currentPile ) )
        {
            if ( verboseEncode == 1 )
                std::cerr << "j-1: Q[" << k << "]=" << ( int )vectTriple[k].pileN << " P[" << k << "]=" << ( dataTypeNChar )vectTriple[k].posN << " N[" << k << "]=" << ( dataTypeNSeq )vectTriple[k].seqN << "\t";

            //std::cerr << "--k= " << k << " pileN[k]= " << pileN[k] << " posN[k]= " << posN[k] <<  " seqN[k]= " << seqN[k] << std::endl;
            //For any character (of differents sequences) in the same pile
            foundSymbol = '\0';

            //cont is the number of symbols already read!
            //      toRead = vectTriple[k].posN - cont -1;



            toRead = vectTriple[k].posN - cont;
            //      cout << "toRead: " << toRead << endl;
            if ( toRead > 0 )
            {
                if ( toRead > 1 )
                {
                    numberRead = ( *pReader ).readAndCount( counters, toRead - 1 );
                    assert ( toRead - 1 == numberRead );
                }
                assert( ( *pReader )( ( char * )&foundSymbol, 1 ) == 1 );
                if ( whichPile[( int )foundSymbol] < alphabetSize - 1 ) {}
                else
                {
                    cout << ( int )foundSymbol << " " << foundSymbol << endl;
                    assert( 1 == 0 );
                }

                counters.count_[whichPile[( int )foundSymbol]]++;
                cont += toRead;
            }

#ifdef DEBUG
            cout << zz << " " << toRead << " " << foundSymbol  << " ";
            counters.print();
#endif

            //std::cerr << "toRead " << toRead << "Found Symbol is " << foundSymbol << "\n";


            //I have to update the value in vectTriple[k].posN, it must contain the position of the new symbol
            //#ifdef XXX
            //vectTriple[k].posN = counters.count_[whichPile[(int)foundSymbol]];
            newVectTripleItem.posN = counters.count_[whichPile[( int )foundSymbol]];
            //#endif
            //std::cerr << "--New posN[k]=" << (int)posN[k] <<std::endl;
            if ( verboseEncode == 1 )
                std::cerr << "\nInit New P[" << k << "]= " << vectTriple[k].posN << std::endl; //TODO: update this to newVectTripleItem

            for ( dataTypedimAlpha g = 0 ; g < currentPile; g++ )  //I have to count in each pile g= 0... (currentPile-1)-pile
            {
                //                vectTriple[k].posN = vectTriple[k].posN + tableOcc_[g].count_[whichPile[(int)foundSymbol]];
                newVectTripleItem.posN += tableOcc_[g].count_[whichPile[( int )foundSymbol]];
                //std::cerr << "--New posN[k]=" << (int)posN[k] << " tableOcc[g][whichPile[(int)symbol]] " << tableOcc[g][whichPile[(int)symbol]] <<std::endl;
                if ( verboseEncode == 1 )
                {
                    std::cerr << "g= " << ( int )g << " symbol= " << ( int )foundSymbol << " whichPile[symbol]= "
                              << ( int )whichPile[( int )foundSymbol] << std::endl;
                    std::cerr << "Add New posN[k]=" << vectTriple[k].posN << " tableOcc[g][whichPile[(int)symbol]] "
                              << tableOcc_[g].count_[whichPile[( int )foundSymbol]] << std::endl;
                }

            }
            //I have to insert the new symbol in the symbol-pile
            assert( whichPile[( int )foundSymbol] < alphabetSize - 1 );
            //            vectTriple[k].pileN=whichPile[(int)foundSymbol];
            newVectTripleItem.pileN = whichPile[( int )foundSymbol];
            //std::cerr << "New posN[k]=" << (int)posN[k] << " New pileN[k]=" << (int)pileN[k] << std::endl;
            if ( verboseEncode == 1 )
                std::cerr << "j  : Q[q]=" << ( int )vectTriple[k].pileN << " P[q]=" << ( dataTypeNChar )vectTriple[k].posN <<  " N[q]=" << ( dataTypeNSeq )vectTriple[k].seqN << std::endl;

            newVectTripleItem.seqN = vectTriple[k].seqN;
            newVectTripleItem.setLcpCurN( vectTriple[k].getLcpCurN() );
            newVectTripleItem.setLcpSucN( vectTriple[k].getLcpSucN() );
            //            vectTriple[k] = newVectTripleItem;
            newVectTriplePerNewPile[ newVectTripleItem.pileN ].push_back( newVectTripleItem );

            k++;
        }
        j = k;

        assert ( j == endIndex ); // while loop removed, as we now only process one pile here
    }

    delete pReader;
    pReader = NULL;
}

void BCRexternalBWT::storeBWT( uchar const *newSymb, uchar const *newQual )
{

    //I have found the position where I have to insert the chars in the position t of the each text
    //Now I have to update the BWT in each file.

    // We first calculate at which index each pile starts
    vector<dataTypeNSeq> pileStarts( alphabetSize );
    pileStarts[0] = 0;
    dataTypeNSeq index = 0;
    for ( int pile = 1; pile < alphabetSize; ++pile )
    {
        while ( index < nText && vectTriple[index].pileN < pile )
            ++index;
        pileStarts[pile] = index;
    }


    int parallelPile;
    //    for (int parallelPile = alphabetSize-2; parallelPile >= 0; --parallelPile)
    #pragma omp parallel for
    for ( parallelPile = 0; parallelPile < alphabetSize - 1; ++parallelPile )
    {
        storeBWT_parallelPile( newSymb, newQual, parallelPile, pileStarts[parallelPile], pileStarts[parallelPile + 1] );
    }


    //static FILE *tmpFile;
    //tmpFile = fopen("sizeFileBwt.txt", "a");
    static FILE *OutFileBWT;                  // output and input file BWT;

    //Renaming new to old
    for ( dataTypedimAlpha g = 0 ; g < alphabetSize - 1; g++ )
    {
        Filename filenameIn( g );
        Filename filenameOut( "new_", g );
        //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
        OutFileBWT = fopen( filenameOut, "rb" );

        if ( OutFileBWT != NULL ) //If it exists
        {
            fclose( OutFileBWT );
            if ( remove( filenameIn ) != 0 )
                std::cerr << filenameIn << ": Error deleting file" << std::endl;
            else if ( rename( filenameOut, filenameIn ) )
                std::cerr << filenameOut << ": Error renaming " << std::endl;
        }

        if ( verboseEncode == 1 )
        {
            struct stat results;
            if ( stat( filenameIn, &results ) == 0 )
                // The size of the file in bytes is in results.st_size
                //fprintf(tmpFile,"%s\t%u\n", filenameIn, results.st_size);
                std::cerr << filenameIn << "\t" << results.st_size << std::endl;
            else
                //fprintf(tmpFile,"An error occurred %s\n", filenameIn);
                std::cerr << "An error occurred" << std::endl;
        }

        const bool permuteQualities = ( bwtParams_->getValue( BWT_OPTION_PROCESS_QUALITIES ) == PROCESS_QUALITIES_PERMUTE );
        if ( permuteQualities )
        {
            Filename filenameQualIn( "", g, ".qual" );
            Filename filenameQualOut( "new_", g, ".qual" );
            FILE *OutQualFileBWT = fopen( filenameQualOut, "rb" );

            if ( OutQualFileBWT != NULL ) //If it exists
            {
                fclose( OutQualFileBWT );
                if ( remove( filenameQualIn ) != 0 )
                    std::cerr << filenameQualIn << ": Error deleting file" << std::endl;
                else if ( rename( filenameQualOut, filenameQualIn ) )
                    std::cerr << filenameQualOut << ": Error renaming " << std::endl;
            }
        }
    }
    //std::cerr <<  std::endl;
    //fprintf(tmpFile,"\n");
    //fclose(tmpFile);
    //  delete pReader;
    //  delete pWriter;

}

char getPredictionBasedEncoding( const char actualBase, char predictedBase, bool &isCorrectlyPredicted )
{
#define USE_PBE_ALGO1
    //#define USE_PBE_ALGO2
    char result;
#ifdef USE_PBE_ALGO1
    if ( predictedBase == notInAlphabet )
        predictedBase = 'A';

    switch ( predictedBase )
    {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            if ( actualBase == predictedBase )
            {
                result = 'A';
            }
            else if ( actualBase == 'A' )
            {
                result = predictedBase;
            }
            else
            {
                result = actualBase;
            }
            break;

        default:
            cerr << "Error: Predicted base = " << predictedBase << endl;
            assert( false && "predicted base should have been A, C, G or T" );
    }

    isCorrectlyPredicted = ( result == 'A' );
    return result;
#endif // USE_PBE_ALGO1


#ifdef USE_PBE_ALGO2
    // from observations: when the insertion is in the middle of a BWT run of the same letter, it is more frequent to have the following insertions: A<->G and C<->T
    char predictedBase2;
    // if (predictedBase == predictedBase2) <- encoding used when we use both the BWT bases before and after the inserting point
    switch ( predictedBase )
    {
        case 'A':
            predictedBase2 = 'G';
            break;
        case 'C':
            predictedBase2 = 'T';
            break;
        case 'G':
            predictedBase2 = 'A';
            break;
        case 'T':
            predictedBase2 = 'C';
            break;
    }

    if ( actualBase == predictedBase )
    {
        result = 'A';
    }
    else if ( actualBase == predictedBase2 )
    {
        result = 'C';
    }
    else if ( actualBase == 'A' )
    {
        if ( predictedBase != 'C' )
            result = predictedBase;
        else
            result = predictedBase2;
    }
    else if ( actualBase == 'C' )
    {
        result = predictedBase2;
    }
    else
    {
        result = actualBase;
    }

    isCorrectlyPredicted = ( result == 'A' );
    return result;
#endif // USE_PBE_ALGO2
}

void BCRexternalBWT::storeBWT_parallelPile( uchar const *newSymb, uchar const *newQual, unsigned int parallelPile, dataTypeNSeq startIndex, dataTypeNSeq endIndex )
{
    if ( 0 )
    {
        #pragma omp critical
        {
            clog << "storeBWT_parallelPile: pile=" << parallelPile << " from " << startIndex << " to " << endIndex << endl;
        }
    }
    //I have found the position where I have to insert the chars in the position t of the each text
    //Now I have to update the BWT in each file.
    dataTypeNChar toRead = 0;
    dataTypeNSeq j;
    dataTypedimAlpha currentPile = parallelPile;
    const bool permuteQualities = ( bwtParams_->getValue( BWT_OPTION_PROCESS_QUALITIES ) == PROCESS_QUALITIES_PERMUTE );

    BwtReaderBase *pReader( NULL );
    BwtWriterBase *pWriter( NULL );
    BwtReaderBase *pQualReader( NULL );
    BwtWriterBase *pQualWriter( NULL );
    BwtWriterBase *pWriterPredictionBasedEncoding( NULL );
    BwtWriterBase *pQualWriterPredictionBasedEncoding( NULL );

    // If there are no letters to add at the last cycle, still process it
    if ( startIndex >= endIndex && debugCycle >= lengthRead && currentPile > 0 )
    {
        Filename filenameIn( "", currentPile, "" );
        Filename filenameOut( "new_", currentPile, "" );
        dumpRamFileToFile( filenameIn, filenameOut );
    }

    j = startIndex;
    while ( j < endIndex )
    {
        assert( currentPile == vectTriple[j].pileN );
        if ( verboseEncode == 1 )
            std::cerr << "index j= " << j << " current BWT segment " << ( int )currentPile << std::endl;

        //std::cerr << "Pile " << (int)currentPile << std::endl;
        Filename filenameIn( "", currentPile, "" );
        Filename filenameOut( "new_", currentPile, "" );

        if ( pReader != NULL ) delete pReader;
        if ( pWriter != NULL ) delete pWriter;

        if ( currentPile == 0 )
        {
            pReader = instantiateBwtReaderForFirstCycle( filenameIn );
            pWriter = instantiateBwtWriterForFirstCycle( filenameOut );
        }
        else if ( debugCycle < lengthRead )
        {
            pReader = instantiateBwtReaderForIntermediateCycle( filenameIn );
            pWriter = instantiateBwtWriterForIntermediateCycle( filenameOut );
        }
        else
        {
            pReader = instantiateBwtReaderForIntermediateCycle( filenameIn );
            pWriter = instantiateBwtWriterForLastCycle( filenameOut );
        }
        assert( pReader != NULL );
        assert( pWriter != NULL );


        if ( permuteQualities )
        {
            if ( pQualReader != NULL ) delete pQualReader;
            if ( pQualWriter != NULL ) delete pQualWriter;
            Filename filenameQualIn( filenameIn, ".qual" );
            Filename filenameQualOut( filenameOut, ".qual" );
            pQualReader = new BwtReaderASCII( filenameQualIn );
            pQualWriter = new BwtWriterASCII( filenameQualOut );
        }

        //For each new symbol in the same pile
        dataTypeNSeq k = j;
        dataTypeNChar cont = 0;
        while ( ( k < nText ) && ( vectTriple[k].pileN == currentPile ) )
        {
            if ( verboseEncode == 1 )
                std::cerr << "k= " << k << " Q[k]= " << ( int )vectTriple[k].pileN << " P[k]= " << vectTriple[k].posN << " cont = " << cont << std::endl;
            //std::cerr << "k= " << k << " pileN[k]= " << pileN[k] << " posN[k]= " << posN[k] << std::endl;
            //So I have to read the k-BWT and I have to count the number of the symbols up to the position posN.
            //As PosN starts to the position 1 and I have to insert the new symbol in position posN[k]
            // I have to read posN[k]-1 symbols
            //cont is the number of symbols already read!
            toRead = ( vectTriple[k].posN - 1 ) - cont;
            if ( verboseEncode == 1 )
                std::cerr << "Start: to Read " << toRead << "\n";

            ( *pReader ).readAndSend( ( *pWriter ), toRead );
            if ( permuteQualities )
            {
                ( *pQualReader ).readAndSend( ( *pQualWriter ), toRead );
            }

            cont += toRead;

            //      toRead=0;


            //Now I have to insert the new symbol associated with the suffix of the sequence k
            //And I have to update the number of occurrences of each symbol
            //   if (toRead==0) {

            const char charToWrite = newSymb[vectTriple[k].seqN];
            ( *pWriter )( ( char * )&charToWrite, 1 );
            if ( permuteQualities )
            {
                assert( newQual );
                ( *pQualWriter )( ( char * )&newQual[vectTriple[k].seqN], 1 );
            }

            tableOcc_[currentPile].count_[whichPile[( int )newSymb[vectTriple[k].seqN]]]++;     //update the number of occurrences in BWT of the pileN[k]
            //std::cerr << "new number write " << numchar << "\n";
            cont++;    //number of read symbols
            //toRead--;
            //      }
            k++;   //  I changed the number of the sequence. New iteration.
        }

        //it means that posN[k]<>currentPile, so I have to change BWT-file
        //But before, I have to copy the remainder symbols from the old BWT to new BWT
        ( *pReader ).readAndSend( ( *pWriter ) );
        if ( permuteQualities )
        {
            ( *pQualReader ).readAndSend( ( *pQualWriter ) );
        }

        j = k;
    }

    delete pReader;
    pReader = NULL; // %%%
    delete pWriter;
    pWriter = NULL; // %%%
    if ( permuteQualities )
    {
        delete pQualReader;
        pQualReader = NULL;
        delete pQualWriter;
        pQualWriter = NULL;
    }

    if ( pWriterPredictionBasedEncoding )
    {
        delete pWriterPredictionBasedEncoding;
        pWriterPredictionBasedEncoding = NULL;
    }
    if ( pQualWriterPredictionBasedEncoding )
    {
        delete pQualWriterPredictionBasedEncoding;
        pWriterPredictionBasedEncoding = NULL;
    }
}

void BCRexternalBWT::storeEntireBWT( const char *fn )
{

    static FILE *OutFileBWT, *InFileBWT;                  // output and input file BWT;
    dataTypeNChar numchar = 0;
    dataTypeNChar numcharWrite = 0;

    uchar *buffer = new uchar[SIZEBUFFER];

    dataTypeNChar *freqOut = new dataTypeNChar [256];
    for ( unsigned i = 0; i < 255; ++i )
        freqOut[i] = 0;

    OutFileBWT = fopen( fn, "wb" );
    if ( OutFileBWT == NULL )
    {
        std::cerr << "storeEntireBWT: Error opening " << std::endl;
        exit ( EXIT_FAILURE );
    }

    if ( verboseEncode == 1 )
    {
        std::cerr << "\nThe last BWT-segment:" << std::endl;
        unsigned int mmm = 0;
        while ( mmm < alphabetSize )
        {
            Filename filenameIn( mmm );
            //printf("===Current BWT-partial= %d\n",mmm);
            InFileBWT = fopen( filenameIn, "rb" );
            for ( dataTypeNChar g = 0 ; g < SIZEBUFFER; g++ )
                buffer[g] = '\0';
            numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, InFileBWT );
            std::cerr << "B[" << ( int )mmm << "]:\t";
            if ( numchar == 0 )
                std::cerr  << "empty";
            else
                std::cerr  << buffer;
            while ( numchar != 0 )
            {
                for ( dataTypeNChar g = 0 ; g < SIZEBUFFER; g++ )
                    buffer[g] = '\0';
                numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, InFileBWT );
                if ( numchar != 0 )
                    std::cerr  << buffer;
            }
            std::cerr << std::endl;

            fclose( InFileBWT );
            mmm++;
        }
    }

    std::cerr << "Entire BWT file" << std::endl;
    std::cerr << "Concatenation of " << ( int )alphabetSize << " segments \n";

    std::cerr << "Compute the distribution of chars \n";

    for ( dataTypedimAlpha g = 0 ; g < alphabetSize; g++ )
    {
        Filename filenameIn( g );
        InFileBWT = fopen( filenameIn, "rb" );
        if ( InFileBWT == NULL )
        {
            std::cerr << "storeEntireBWT: " << "BWT file " << ( int )g << ": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        //std::cerr << "BWT file " << (int)g << "= ";
        while ( numchar != 0 )
        {
            numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, InFileBWT );
            //std::cerr << "number read " << numchar << "\n";
            numcharWrite = fwrite ( buffer, sizeof( uchar ), numchar , OutFileBWT );
            checkIfEqual( numchar, numcharWrite ); // we should always read/write the same number of characters

            for ( unsigned j = 0 ; j < numchar; j++ )
                freqOut[( int )( buffer[j] )]++;
        }

        fclose( InFileBWT );
    }
    fclose( OutFileBWT );

    if ( verboseEncode == 1 )
    {
        std::cerr << "\nThe Entire BWT:" << std::endl;
        OutFileBWT = fopen( fn, "rb" );
        if ( OutFileBWT == NULL )
        {
            std::cerr << "storeEntireBWT: Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }

        for ( dataTypeNSeq g = 0 ; g < SIZEBUFFER; g++ )
            buffer[g] = '\0';
        numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, OutFileBWT );
        if ( numchar == 0 )
            std::cerr  << "empty\n";
        else
            std::cerr  << buffer << "\n";
        fclose( OutFileBWT );
    }
    delete [] buffer;


    std::cerr << "Distribution in BWT\n";
    for ( dataTypedimAlpha i = 0; i < 255; ++i )
        if ( freqOut[i] > 0 )
            std::cerr << i << " " << freqOut[i] << "\n";
    delete [] freqOut;
}

void BCRexternalBWT::storeSA( dataTypelenSeq iterationNum )
{

    //I have found the position where I have to insert the chars in the position t of the each text
    //Now I have to update the SA in each file.
    static FILE *OutFileSA, *InFileSA;                  // output and input file SA;

    dataTypeNChar numchar = 0;
    dataTypeNChar numcharWrite = 0;
    ElementType *buffer = new ElementType[SIZEBUFFER];
    dataTypeNChar toRead = 0;

    dataTypeNSeq j;
    dataTypedimAlpha currentPile;
    j = 0;
    while ( j < nText )
    {
        currentPile = vectTriple[j].pileN;
        //if (verboseEncode==1)
        // std::cerr << "\nNew Segment; index text j= " << j << " current SA segment is " << (int)currentPile << std::endl;
        //std::cerr << "Pile " << (int)currentPile << std::endl;
        Filename filenameIn( "sa_", currentPile );
        InFileSA = fopen( filenameIn, "rb" );
        if ( InFileSA == NULL )
        {
            std::cerr << "In SA file " << ( int )j << ": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }

        Filename filenameOut( "new_sa_", currentPile );
        OutFileSA = fopen( filenameOut, "wb" );
        if ( OutFileSA == NULL )
        {
            std::cerr << "Out SA file " << ( int )j << ": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        //std::cerr << "In File " << filenameIn << std::endl;
        //std::cerr << "Out File " << filenameOut << std::endl;

        //For each new symbol in the same pile
        dataTypeNSeq k = j;
        dataTypeNChar cont = 0;
        while ( ( k < nText ) && ( vectTriple[k].pileN == currentPile ) )
        {

            //if (verboseEncode==1)
            // std::cerr << "k= " << k << " Q[k]= " << (int)vectTriple[k].pileN << " P[k]= " << vectTriple[k].posN << " cont = "<< cont << std::endl;
            //std::cerr << "k= " << k << " pileN[k]= " << pileN[k] << " posN[k]= " << posN[k] << std::endl;
            //So I have to read the k-SA and I have to count the number of the symbols up to the position posN.
            //As PosN starts to the position 1 and I have to insert the new symbol in position posN[k]
            // I have to read posN[k]-1 symbols
            //cont is the number of symbols already read!
            toRead = ( vectTriple[k].posN - 1 ) - cont;
            /*
            if (verboseEncode == 1)
                std::cerr << "Start: to Read " << toRead << "\n";
            */
            while ( toRead > 0 )            //((numchar!=0) && (toRead > 0)) {
            {
                if ( toRead < SIZEBUFFER ) //The last reading for this sequence
                {
                    numchar = fread( buffer, sizeof( ElementType ), toRead, InFileSA );
                    /*
                    if (verboseEncode == 1)
                        std::cerr << "number read " << numchar << " to Read " << toRead << "\n";
                    */
                    checkIfEqual( numchar, toRead ); // we should always read/write the same number of characters

                    numcharWrite = fwrite ( buffer, sizeof( ElementType ), numchar , OutFileSA );
                    checkIfEqual( numchar, numcharWrite ); // we should always read/write the same number of characters
                    //std::cerr << "toread number write " << numcharWrite << "\n";
                }
                else
                {
                    numchar = fread( buffer, sizeof( ElementType ), SIZEBUFFER, InFileSA );
                    //if (verboseEncode == 1)
                    // std::cerr << "number read " << numchar << "\n";
                    checkIfEqual( numchar, SIZEBUFFER ); // we should always read/write the same number of characters
                    numcharWrite = fwrite ( buffer, sizeof( ElementType ), numchar , OutFileSA );
                    checkIfEqual( numchar , numcharWrite ); // we should always read/write the same number of characters
                    //std::cerr << "sizebuffer number write " << numcharWrite << "\n";
                }

                cont   += numchar;  //number of read symbols
                toRead -= numchar;
                if ( ( numchar == 0 ) && ( toRead > 0 ) ) //it means that we have read 0 character, but there are still toRead characters to read
                {
                    std::cerr << "storeSA: sequence number" << ( int )k << " read 0 character, but there are still " << toRead << " characters to read  " << std::endl;
                    exit ( EXIT_FAILURE );
                }
            }
            //Now I have to insert the new symbol associated with the suffix of the sequence k
            //And I have to update the number of occurrences of each symbol
            if ( toRead == 0 )
            {
                ElementType newEle;
                newEle.sa = iterationNum; //( posSymb + 1 ) % ( lengthRead + 1 );
                newEle.numSeq = vectTriple[k].seqN;

                numchar = fwrite ( &newEle, sizeof( ElementType ), 1, OutFileSA );
                checkIfEqual( numchar, 1 ); // we should always read/write the same number of characters
                //it is not useful for the suffix array
                //tableOcc[currentPile][alpha[(int)newSymb[vectTriple[k].seqN]]]++;       //update the number of occurrences in SA of the pileN[k]
                //std::cerr << "new number write " << numchar << "\n";
                cont++;    //number of read symbols
                toRead--;
            }

            k++;   //  I changed the number of the sequence. New iteration.
        }

        //it means that posN[k]<>currentPile, so I have to change SA-file
        //But before, I have to copy the remainder symbols from the old SA to new SA
        while ( numchar != 0 )
        {
            numchar = fread( buffer, sizeof( ElementType ), SIZEBUFFER, InFileSA );
            //std::cerr << "After insert: " << numchar << "\n";
            numcharWrite = fwrite ( buffer, sizeof( ElementType ), numchar , OutFileSA );
            checkIfEqual( numchar, numcharWrite ); // we should always read/write the same number of characters
        }

        fclose( InFileSA );
        fclose( OutFileSA );
        j = k;
    }

    //Renaming new to old
    for ( dataTypedimAlpha g = 0 ; g < alphabetSize; g++ )
    {
        Filename filenameIn( "sa_", g );
        Filename filenameOut( "new_sa_", g );
        //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
        OutFileSA = fopen( filenameOut, "rb" );

        if ( OutFileSA != NULL ) //If it exists
        {
            fclose( OutFileSA );
            if ( remove( filenameIn ) != 0 )
                std::cerr << filenameIn << ": Error deleting file" << std::endl;
            else if ( rename( filenameOut, filenameIn ) )
                std::cerr << filenameOut << ": Error renaming " << std::endl;
        }
        /*
        if (verboseEncode == 1) {
            struct stat results;
            if (stat(filenameIn, &results) == 0)
                // The size of the file in bytes is in results.st_size
                //fprintf(tmpFile,"%s\t%u\n", filenameIn, results.st_size);
                std::cerr << filenameIn <<"\t" << results.st_size << std::endl;
            else
                //fprintf(tmpFile,"An error occurred %s\n", filenameIn);
                std::cerr << "An error occurred" << std::endl;
        }
        */
    }
    //std::cerr <<  std::endl;
    //fprintf(tmpFile,"\n");
    //fclose(tmpFile);

    delete [] buffer;
}

void BCRexternalBWT::storeEntirePairSA( const char *fn )
{

    std::cerr << "\nEntire Pairs SA file (position, number of sequence)" << std::endl;

    static FILE *OutFileSA, *InFileSA;                  // output and input file SA;
    dataTypeNChar numcharWrite, numcharRead;
    ElementType *buffer = new ElementType[SIZEBUFFER];

    Filename fnSA( fn, ".pairSA" );

    OutFileSA = fopen( fnSA, "wb" );
    if ( OutFileSA == NULL )
    {
        std::cerr << "Entire Pairs SA file: Error opening " << fnSA << std::endl;
        exit ( EXIT_FAILURE );
    }
    /* //it will be useful for varying length reads
        vector <dataTypelenSeq> vectLen;
        vectLen.resize(nText);

        char *fileLen="outFileLen";
        static FILE *InFileLen;                  // file of the lengths;
        InFileLen = fopen(fileLen, "rb");
        if (InFileLen==NULL) {
                std::cerr << "storeEntireSAfromPairSA: could not open file \"" << fileLen << "\"!"<< std::endl;
                exit (EXIT_FAILURE);
        }

        numcharRead = fread (&vectLen[0], sizeof(dataTypelenSeq), vectLen.size() , InFileLen);
        checkIfEqual(numcharRead , nText); // we should always read the same number of characters

        fclose(InFileLen);
        */
    for ( dataTypedimAlpha g = 0 ; g < alphabetSize; g++ )
    {
        Filename filenameIn( "sa_", g );
        InFileSA = fopen( filenameIn, "rb" );
        if ( InFileSA == NULL )
        {
            std::cerr << "SA file " << ( int )g << ": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        //std::cerr << "SA file " << (int)g << "= ";

        numcharRead = fread( buffer, sizeof( ElementType ), SIZEBUFFER, InFileSA );

        /* //it will be useful for varying length reads
        //Correction of the length of the sequences.
        for (dataTypeNSeq num = 0; num < numcharRead; num++) {
            buffer[num].sa = buffer[num].sa - (lengthRead - vectLen[buffer[num].numSeq]);
        }
        */

        numcharWrite = fwrite ( buffer, sizeof( ElementType ), numcharRead , OutFileSA );
        checkIfEqual ( numcharRead , numcharWrite );

        while ( numcharRead != 0 )
        {
            numcharRead = fread( buffer, sizeof( ElementType ), SIZEBUFFER, InFileSA );
            /* //it will be useful for varying length reads
            //Correction of the length of the sequences.
            for (dataTypeNSeq num = 0; num < numcharRead; num++) {
                buffer[num].sa = buffer[num].sa - (lengthRead - vectLen[buffer[num].numSeq]);
            }
            */
            numcharWrite = fwrite ( buffer, sizeof( ElementType ), numcharRead , OutFileSA );
            checkIfEqual ( numcharRead , numcharWrite );
        }

        fclose( InFileSA );
        if ( remove( filenameIn ) != 0 )
            std::cerr << filenameIn << ": Error deleting file" << std::endl;
    }

    fclose( OutFileSA );

    if ( verboseEncode == 1 )
    {
        OutFileSA = fopen( fnSA, "rb" );
        if ( OutFileSA == NULL )
        {
            std::cerr << "Entire SA file: Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }

        numcharRead = fread( buffer, sizeof( ElementType ), SIZEBUFFER, OutFileSA );
        if ( numcharRead == 0 )
            std::cerr  << "empty\n";
        else
            for ( dataTypeNChar g = 0 ; g < numcharRead; g++ )
            {
                std::cerr  << "(" << ( int )buffer[g].sa << "," << buffer[g].numSeq << ") ";
            }
        while ( numcharRead != 0 )
        {
            numcharRead = fread( buffer, sizeof( ElementType ), SIZEBUFFER, OutFileSA );
            for ( dataTypeNChar g = 0 ; g < numcharRead; g++ )
            {
                std::cerr  << "(" << buffer[g].sa << "," << buffer[g].numSeq << ") ";
            }
        }
        std::cerr << std::endl;

        fclose( OutFileSA );
    }

    delete [] buffer;
}


void BCRexternalBWT::storeEntireSAfromPairSA( const char *fn )
{
    static FILE *OutFileSA, *InFilePairSA;                  // output and input file SA;

    std::cerr << "\nSA file from pair SA file" << std::endl;

    dataTypeNChar numchar, numcharWrite;
    /* //it will be useful for varying length reads

        vector <dataTypelenSeq> vectSumCumLen;
        vectSumCumLen.resize(nText+1);
        char *fileLen="outFileLen";
        static FILE *InFileLen;                  // file of the lengths;
        InFileLen = fopen(fileLen, "rb");
        if (InFileLen==NULL) {
                std::cerr << "storeEntireSAfromPairSA: could not open file \"" << fileLen << "\"!"<< std::endl;
                exit (EXIT_FAILURE);
        }
        dataTypelenSeq lenSeq=0;
        numchar = fread (&lenSeq, sizeof(dataTypelenSeq), 1 , InFileLen);
        checkIfEqual( numchar , 1); // we should always read the same number of characters
        vectSumCumLen[0] = 0;
        vectSumCumLen[1] = lenSeq + 1;   //Plus $
        for (dataTypeNSeq num = 2; num < nText+1; num++) {
            numchar = fread (&lenSeq, sizeof(dataTypelenSeq), 1 , InFileLen);
            checkIfEqual(numchar , 1); // we should always read the same number of characters
            vectSumCumLen[num] = vectSumCumLen[num-1] + lenSeq + 1;  //Plus $
        }
        fclose(InFileLen);
    */
    Filename fnSA    ( fn, ".sa" );
    Filename fnPairSA( fn, ".pairSA" );

    InFilePairSA = fopen( fnPairSA, "rb" );
    if ( InFilePairSA == NULL )
    {
        std::cerr << "Entire Pairs SA file: Error opening " << fnPairSA << std::endl;
        exit ( EXIT_FAILURE );
    }

    OutFileSA = fopen( fnSA, "wb" );
    if ( OutFileSA == NULL )
    {
        std::cerr << "Entire SA file: Error opening " << fnSA << std::endl;
        exit ( EXIT_FAILURE );
    }

    ElementType *buffer = new ElementType[SIZEBUFFER];
    dataTypeNChar *bufferNChar = new dataTypeNChar[SIZEBUFFER];

    while ( !feof( InFilePairSA ) )
    {
        numchar = fread( buffer, sizeof( ElementType ), SIZEBUFFER, InFilePairSA );
        //std::cerr << "number read " << numchar << "\n";
        if ( numchar > 0 )
        {
            for ( dataTypeNChar i = 0; i < numchar; i++ )
            {
                bufferNChar[i] = ( dataTypeNChar )( buffer[i].numSeq * ( lengthRead + 1 ) + buffer[i].sa );
                //std::cerr << buffer[i].numSeq << " " << (int)lengthRead << " " << (int)buffer[i].sa << " --> " << (int)bufferNChar[i] << "\n";
                //bufferNChar[i] = (dataTypeNChar)(vectSumCumLen[buffer[i].numSeq] + buffer[i].sa);       //it will be useful for varying length reads
                //std::cerr << "vectSumCumLen["<< buffer[i].numSeq<< "]= " << (int)vectSumCumLen[buffer[i].numSeq] << " + " << (int)buffer[i].sa << " --> " << (int)bufferNChar[i] << "\n";
            }
            numcharWrite = fwrite ( bufferNChar, sizeof( dataTypeNChar ), numchar, OutFileSA );
            checkIfEqual( numchar , numcharWrite );
            //std::cerr << "number write " << numcharWrite << "\n";
        }
    }
    fclose( InFilePairSA );
    fclose( OutFileSA );

    if ( verboseEncode == 1 )
    {
        std::cerr << "\nThe Entire SA. The file is " << fnSA << std::endl;
        OutFileSA = fopen( fnSA, "rb" );
        if ( OutFileSA == NULL )
        {
            std::cerr << "Entire SA file: Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }

        numchar = fread( bufferNChar, sizeof( dataTypeNChar ), SIZEBUFFER, OutFileSA );
        if ( numchar == 0 )
            std::cerr  << "empty\n";
        else
            for ( dataTypeNChar g = 0 ; g < numchar; g++ )
            {
                std::cerr  << bufferNChar[g] << " ";
            }
        while ( numchar != 0 )
        {
            numchar = fread( buffer, sizeof( ElementType ), SIZEBUFFER, OutFileSA );
            for ( dataTypeNChar g = 0 ; g < numchar; g++ )
            {
                std::cerr  << bufferNChar[g] << " ";
            }
        }
        std::cerr << std::endl;

        fclose( OutFileSA );
    }

    delete [] buffer;
    delete [] bufferNChar;
}

void BCRexternalBWT::storeBWTandLCP( uchar const *newSymb )
{
    dataTypelenSeq maxValueLen = lengthRead + 2;
    vector <dataTypelenSeq> minLCPcur;
    vector <bool> minLCPcurFound;
    vector <dataTypelenSeq> minLCPsuc;
    vector <dataTypeNSeq> minLCPsucText;
    vector <bool> minLCPsucToFind;
    minLCPcur.resize( alphabetSize );     //for each symbol of the alphabet
    minLCPcurFound.resize( alphabetSize );
    minLCPsuc.resize( alphabetSize );
    minLCPsucText.resize( alphabetSize );
    minLCPsucToFind.resize( alphabetSize );

    //I have found the position where I have to insert the chars in the position t of the each text
    //Now I have to update the BWT in each file.
    static FILE *OutFileBWT, *InFileBWT;                  // output and input file BWT;
    static FILE *OutFileLCP, *InFileLCP;                  // output and input file LCP;

    dataTypeNChar numchar = 0;
    dataTypeNChar numcharWrite = 0;
    uchar *buffer = new uchar[SIZEBUFFER];
    dataTypelenSeq *bufferLCP = new dataTypelenSeq[SIZEBUFFER];
    dataTypeNChar toRead = 0;

    dataTypeNSeq j;
    dataTypedimAlpha currentPile;
    //uchar symbol='\0';

    j = 0;
    while ( j < nText )
    {
        currentPile = vectTriple[j].pileN;
        for ( dataTypedimAlpha g = 0 ; g < alphabetSize; g++ )
        {
            minLCPcur[g] = maxValueLen;
            minLCPcurFound[g] = 0;
            minLCPsuc[g] = maxValueLen;
            minLCPsucText[g] = 0;   //denotes the number of the text associated with the symbol g
            minLCPsucToFind[g] = 0;
        }
        assert( currentPile <= alphabetSize );
        Filename filenameIn( currentPile );
        InFileBWT = fopen( filenameIn, "rb" );
        if ( InFileBWT == NULL )
        {
            std::cerr << "(storeBWTandLCP) In BWT file " << filenameIn << ": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        Filename filenameOut( "new_", currentPile );
        OutFileBWT = fopen( filenameOut, "wb" );
        if ( OutFileBWT == NULL )
        {
            std::cerr << "(storeBWTandLCP) Out BWT file " << filenameIn << ": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        Filename filenameInLCP( "", currentPile, ".lcp" );
        InFileLCP = fopen( filenameInLCP, "rb" );
        if ( InFileLCP == NULL )
        {
            std::cerr << "In LCP file " << filenameInLCP << ": Error opening "  << std::endl;
            exit ( EXIT_FAILURE );
        }
        Filename filenameOutLCP( "new_", currentPile, ".lcp" );
        OutFileLCP = fopen( filenameOutLCP, "wb" );
        if ( OutFileLCP == NULL )
        {
            std::cerr << "Out LCP file " << filenameInLCP << ": Error opening " << filenameOutLCP << std::endl;
            exit ( EXIT_FAILURE );
        }

        //For each new symbol in the same pile
        dataTypeNSeq k = j;
        dataTypeNChar cont = 0;
        while ( ( k < nText ) && ( vectTriple[k].pileN == currentPile ) )
        {
            //So I have to read the k-BWT and I have to count the number of the symbols up to the position posN.
            //symbol = '\0';
            //PosN is indexed from the position 1 and I have to insert the new symbol in position posN[k], then I have to read posN[k]-1 symbols
            //cont is the number of symbols already read!
            toRead = ( vectTriple[k].posN - 1 ) - cont;
            while ( toRead > 0 )            //((numchar!=0) && (toRead > 0)) {
            {
                if ( toRead < SIZEBUFFER ) //The last reading for this sequence
                {
                    numchar = fread( buffer, sizeof( uchar ), toRead, InFileBWT );
                    checkIfEqual( numchar , toRead );
                    numcharWrite = fwrite ( buffer, sizeof( uchar ), numchar , OutFileBWT );
                    checkIfEqual( numchar , numcharWrite );
                    numchar = fread( bufferLCP, sizeof( dataTypelenSeq ), toRead, InFileLCP );
                    checkIfEqual( numchar , toRead );
                    numcharWrite = fwrite ( bufferLCP, sizeof( dataTypelenSeq ), numchar , OutFileLCP );
                    checkIfEqual( numchar , numcharWrite );
                }
                else
                {
                    numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, InFileBWT );
                    checkIfEqual( numchar , SIZEBUFFER );
                    numcharWrite = fwrite ( buffer, sizeof( uchar ), numchar , OutFileBWT );
                    checkIfEqual( numchar , numcharWrite );
                    numchar = fread( bufferLCP, sizeof( dataTypelenSeq ), SIZEBUFFER, InFileLCP );
                    checkIfEqual( numchar , SIZEBUFFER );
                    numcharWrite = fwrite ( bufferLCP, sizeof( dataTypelenSeq ), numchar , OutFileLCP );
                    checkIfEqual( numchar , numcharWrite );
                }
                //I must to compute the minimum LCP. It needs to compute the lcpValue for the next iteration
                //std::cerr << "For each letter in the buffer before of the position where I have to insert the new symbol\n";
                for ( dataTypeNChar bb = 0 ; bb < numcharWrite; bb++ )
                {
                    //Update the min1 for each letter of the alphabet, for which I have already met the symbol
                    for ( dataTypedimAlpha gg = 0 ; gg < alphabetSize; gg++ )
                    {
                        if ( minLCPcurFound[gg] == 1 ) //I already met the symbol gg. So, I must compute the minimum
                            if ( bufferLCP[bb] < minLCPcur[gg] ) //comparison with the last inserted lcp
                                minLCPcur[gg] = bufferLCP[bb];
                    }

                    minLCPcur[whichPile[( int )buffer[bb]]] = maxValueLen; //For each occurrence of buffer[bb], I have to set the min1 (the interval starts from the next symbol)
                    minLCPcurFound[whichPile[( int )buffer[bb]]] = 1; //So I open the LCP interval for buffer[bb] (for the first occurrence of buffer[bb])

                    //First, it needs to check if the symbol buffer[bb] closes a previous interval or it is in the middle or no.
                    //In any case, for each symbol, we have to update the minimum
                    for ( dataTypedimAlpha gg = 0 ; gg < alphabetSize; gg++ )
                    {
                        if ( minLCPsucToFind[( int )gg] == 1 )  //We have to compute the minimum for the lcp of the symbol gg
                        {
                            if ( bufferLCP[bb] < minLCPsuc[( int )gg] ) //comparison with the last inserted lcp
                                minLCPsuc[( int )gg] = bufferLCP[bb];
                        }
                    }
                    //If the symbol buffer[bb] does not close an LCP interval, ok!
                    if ( minLCPsucToFind[whichPile[( int )buffer[bb]]] == 1 ) //The symbol buffer[bb] closed an interval LCP
                    {
                        //We have to compute the minimum for the lcp of the symbol buffer[bb]
                        //I have already computed the minimum (close the previous lcp interval).
                        //I can set lcpSucN of minLCPsucText[alpha[(int)[buffer[bb]]]
                        vectTriple[minLCPsucText[whichPile[( int )buffer[bb]]]].setLcpSucN( minLCPsuc[whichPile[( int )buffer[bb]]] );
                        //Since it closes the LCP interval, then
                        minLCPsuc[whichPile[( int )buffer[bb]]] = maxValueLen;
                        minLCPsucToFind[whichPile[( int )buffer[bb]]] = 0;
                        minLCPsucText[whichPile[( int )buffer[bb]]] = 0;
                    }
                }
                cont   += numchar;  //number of read symbols
                toRead -= numchar;
                if ( ( numchar == 0 ) && ( toRead > 0 ) ) //it means that we have read 0 character, but there are still toRead characters to read
                {
                    // NO, abort program
                    std::cerr << "Error storeBWT: sequence number" << ( int )k << " read 0 character, but there are still " << toRead << " characters to read  " << std::endl;
                    exit ( EXIT_FAILURE );
                }

            }
            //Now I have to insert the new symbol associated with the suffix of the sequence k
            //And I have to update the number of occurrences of each symbol
            //And I have to insert the valueLCP store in lcpCurN + 1 in the previous iteration
            if ( toRead == 0 )
            {
                //std::cerr << "\nNow I can insert the new symbol and lcp, indeed toRead= " << toRead << std::endl;
                numchar = fwrite ( &newSymb[vectTriple[k].seqN], sizeof( uchar ), 1, OutFileBWT );
                checkIfEqual( numchar , 1 ); // we should always read/write the same number of characters
                tableOcc_[currentPile].count_[whichPile[( int )newSymb[vectTriple[k].seqN]]]++;
                //tableOcc[currentPile][whichPile[(int)newSymb[vectTriple[k].seqN]]]++;       //update the number of occurrences in BWT of the pileN[k]
                dataTypelenSeq lcpValueNow;
                if ( vectTriple[k].posN == 1 )   //it is the first symbol of the segment. So, the lcp value is 0
                {
                    lcpValueNow = vectTriple[k].getLcpCurN();
                }
                else
                    lcpValueNow = vectTriple[k].getLcpCurN() + 1;
                numchar = fwrite ( &lcpValueNow, sizeof( dataTypelenSeq ), 1, OutFileLCP ); //Insert the lcp for the new symbol
                checkIfEqual( numchar , 1 );
                //std::cerr << "I insert the symbol= " << newSymb[vectTriple[k].seqN] <<  " and lcp " << lcpValueNow << std::endl;
                //Update the lcpCurN for the next iteration
                if ( minLCPcurFound[whichPile[( int )newSymb[vectTriple[k].seqN]]] == 0 )
                {
                    if ( minLCPcur[whichPile[( int )newSymb[vectTriple[k].seqN]]] == maxValueLen )
                    {
                        //it means that we have not met the symbol before of the position posN because minLCPcurFound is 0.
                        //if minLCPcurFound is 0, then minLCPcur is maxValueLen
                        vectTriple[k].setLcpCurN( 0 );         //The next suffix has suffix 0+1=1
                    }
                }
                else    //it means that (minLCPcurFound[alpha[(int)newSymb[vectTriple[k].seqN]]] == 1)
                {
                    if ( minLCPcur[whichPile[( int )newSymb[vectTriple[k].seqN]]] == maxValueLen )
                    {
                        //it means that we have met the symbol before of the position posN because minLCPcurFound is 1.
                        //But minLCPcur is maxValueLen, this means that the previous occurrences of new symbol is the previous position
                        vectTriple[k].setLcpCurN( lcpValueNow );         //The next suffix has suffix lcpValueNow+1
                    }
                    else
                    {
                        if ( lcpValueNow < minLCPcur[whichPile[( int )newSymb[vectTriple[k].seqN]]] )
                        {
                            //comparison with the last inserted lcp. It means that the previous occurrence of new symbol is not the previous symbol
                            vectTriple[k].setLcpCurN( lcpValueNow );
                        }
                        else
                        {
                            vectTriple[k].setLcpCurN( minLCPcur[whichPile[( int )newSymb[vectTriple[k].seqN]]] );
                        }
                    }
                }

                //it may happen that the new symbol closes a previous interval or it is in the middle or no.
                //In any case, for each symbol, we have to update the minimum
                for ( dataTypedimAlpha gg = 0 ; gg < alphabetSize; gg++ )
                {
                    if ( minLCPsucToFind[( int )gg] == 1 )  //We have to compute the minimum for the lcp of the symbol gg
                    {
                        //if (minLCPsuc[gg] != maxValueLen) {      //There are an LCP interval opened for the symbol c_g
                        if ( lcpValueNow < minLCPsuc[( int )gg] ) //comparison with the last inserted lcp
                            minLCPsuc[( int )gg] = lcpValueNow;
                        //}
                    }

                    if ( minLCPcurFound[( int )gg] == 1 )  //We have to compute the minimum for the lcp of the symbol gg
                    {
                        if ( lcpValueNow < minLCPcur[( int )gg] ) //comparison with the last inserted lcp
                            minLCPcur[( int )gg] = lcpValueNow;
                    }
                }

                //I have to re-set
                minLCPcur[whichPile[( int )newSymb[vectTriple[k].seqN]]] = maxValueLen;  //I initialized the minLCPcur
                minLCPcurFound[whichPile[( int )newSymb[vectTriple[k].seqN]]] = 1;  // I set the fact that I met the new symbol

                if ( minLCPsucToFind[whichPile[( int )newSymb[vectTriple[k].seqN]]] == 1 ) //If the new symbol closes a LCP interval
                {
                    //I have already computed the minimum in the previous FOR (close the previous lcp interval).
                    //I can set lcpSucN of minLCPsucText[alpha[(int)[vectTriple[k].seqN]]]
                    vectTriple[minLCPsucText[whichPile[( int )newSymb[vectTriple[k].seqN]]]].setLcpSucN( minLCPsuc[whichPile[( int )newSymb[vectTriple[k].seqN]]] );
                }
                //It set minLVP for the new LCP interval for the new symbol
                //					minLCPsuc[alpha[(int)newSymb[vectTriple[k].seqN]]] = vectTriple[k].lcpCurN+1;   //It sets the min_2 for successive symbol with the current of the new symbol (next text)
                minLCPsuc[whichPile[( int )newSymb[vectTriple[k].seqN]]] = maxValueLen;
                minLCPsucText[whichPile[( int )newSymb[vectTriple[k].seqN]]] = k;
                minLCPsucToFind[whichPile[( int )newSymb[vectTriple[k].seqN]]] = 1; //I have to open the lcp succ for this sequence
                //std::cerr << "minSuc was not the maxValue. The new value for minLCPsuc[" << newSymb[vectTriple[k].seqN] << "] is " << minLCPsuc[alpha[(int)newSymb[vectTriple[k].seqN]]] << "\n";

                //Since we have inserted, we must update them
                cont++;    //number of read symbols
                toRead--;

                //Now I have to update the lcp of the next symbol, if it exists.
                if ( ( k + 1 < nText ) && ( vectTriple[k + 1].pileN == currentPile ) && ( vectTriple[k + 1].posN == vectTriple[k].posN + 1 ) )
                {
                    //If the next symbol is the new symbol associated with another text (in the same segment), that I have to insert yet
                    //I can ignored this updating, because the lcp of the new symbol is already computed
                    //We set the minLCPsuc with the value LCP associated with the next symbol that we have to insert in it.
                    //it should be vectTriple[k].lcpSucN + 1 == vectTriple[k+1].lcpCurN
                    if ( vectTriple[k].getLcpSucN() + 1 != vectTriple[k + 1].getLcpCurN() + 1 )
                    {
                        std::cerr << "???? Warning!--Should be the same? triple[" << k << "].lcpSucN(=" << vectTriple[k].getLcpSucN() << ") + 1= " << vectTriple[k].getLcpSucN() + 1 << " == triple[" << k + 1 << "].lcpCurN+1= " << vectTriple[k + 1].getLcpCurN() + 1 << " ";
                        std::cerr << ", Seq k N. " << vectTriple[k].seqN << " and Seq k+1 N. " << vectTriple[k + 1].seqN << "\n";
                    }

                    //Hence, at the next step, I have to insert the symbol newSymb[vectTriple[k+1].seqN]
                    //I check if newSymb[vectTriple[k+1].seqN] is equal to the inserted symbol now, that is newSymb[vectTriple[k].seqN]
                    if ( newSymb[vectTriple[k].seqN] == newSymb[vectTriple[k + 1].seqN] )
                    {
                        //In this case, I can set the lcpSuc of newSymb[vectTriple[k].seqN] to vectTriple[k+1].lcpCurN + 1
                        vectTriple[k].setLcpSucN( vectTriple[k + 1].getLcpCurN() + 1 );    //I set the lcpSucN of the current symbol (in position k)
                        //I close the LCP interval for newSymb[vectTriple[k].seqN]], so
                        minLCPsuc[whichPile[( int )newSymb[vectTriple[k].seqN]]] = maxValueLen;
                        minLCPsucToFind[whichPile[( int )newSymb[vectTriple[k].seqN]]] = 0; //closes the LCP interval
                        minLCPsucText[whichPile[( int )newSymb[vectTriple[k].seqN]]] = 0;
                    }
                    else
                    {
                        //In this case, I cannot set the lcpSuc of newSymb[vectTriple[k].seqN], because the symbol corresponding to k+1 is different
                        //I set minLCPsuc of newSymb[vectTriple[k].seqN] to vectTriple[k+1].lcpCurN +1, and I search the symbol newSymb[vectTriple[k].seqN]
                        minLCPsuc[whichPile[( int )newSymb[vectTriple[k].seqN]]] = vectTriple[k + 1].getLcpCurN() + 1;	//set the lcp interval
                        minLCPsucToFind[whichPile[( int )newSymb[vectTriple[k].seqN]]] = 1;
                        minLCPsucText[whichPile[( int )newSymb[vectTriple[k].seqN]]] = k;
                    }
                }
                else
                {
                    //The next symbol in the current segment, if there exists, it is an old symbol
                    uchar sucSymbol = '\0';
                    //it there exists another symbol in the segment file then I have to update it.
                    numchar = fread( &sucSymbol, sizeof( uchar ), 1, InFileBWT );
                    if ( numchar == 1 )  //There exists at least a symbol in the current segment
                    {
                        numcharWrite = fwrite ( &sucSymbol, sizeof( uchar ), numchar, OutFileBWT );
                        assert( numchar == numcharWrite ); // we should always read/write the same number of characters
                        //I have to update the lcp of the next symbol
                        dataTypelenSeq sucLCP = 0;
                        numchar = fread( &sucLCP, sizeof( dataTypelenSeq ), 1, InFileLCP ); //I have to change it
                        checkIfEqual( numchar , 1 ); // we should always read/write the same number of characters

                        //I have to update the lcp of this symbol and I have to copy it into the new bwt segment
                        dataTypelenSeq lcpValueNow = vectTriple[k].getLcpSucN() + 1;
                        numcharWrite = fwrite ( &lcpValueNow , sizeof( dataTypelenSeq ), numchar , OutFileLCP ); //Updated the lcpSuc
                        checkIfEqual( numchar , numcharWrite ); // we should always read/write the same number of characters

                        //Now, I have to check if the symbol sucSymbol close the LCP interval the new symbol
                        if ( newSymb[vectTriple[k].seqN] == sucSymbol )
                        {
                            //If it is equal to newSymb[vectTriple[k].seqN] then I can set lcpSucN of newSymb[vectTriple[k].seqN]
                            vectTriple[k].setLcpSucN( lcpValueNow );
                            minLCPsucToFind[whichPile[( int )newSymb[vectTriple[k].seqN]]] = 0; //Close the LCP interval
                            minLCPsuc[whichPile[( int )newSymb[vectTriple[k].seqN]]] = maxValueLen;
                            minLCPsucText[whichPile[( int )newSymb[vectTriple[k].seqN]]] = 0;
                        }
                        else    			//std::cerr << "The succSymb is not equal to the new symbol\n";
                        {
                            minLCPsucToFind[whichPile[( int )newSymb[vectTriple[k].seqN]]] = 1; //I have to search the symbol newSymb[vectTriple[k].seqN]]
                            minLCPsuc[whichPile[( int )newSymb[vectTriple[k].seqN]]] = lcpValueNow; //I set the minLCPsuc for the new symbol
                            minLCPsucText[whichPile[( int )newSymb[vectTriple[k].seqN]]] = k;
                            //It could close the LCP interval for the symbol sucSymb if it is opened
                            //If the symbol sucSymbol does not close an LCP interval, ok!
                            if ( minLCPsucToFind[whichPile[( int )sucSymbol]] == 1 )  //We have to compute the minimum for the lcp of the symbol (int)sucSymbol
                            {
                                //it means that there is an interval lcp to close for the symbol (int)sucSymbol
                                //The symbol sucSymbol closes a LCP interval
                                if ( lcpValueNow < minLCPsuc[whichPile[( int )sucSymbol]] ) //comparison with the last inserted lcp
                                    minLCPsuc[whichPile[( int )sucSymbol]] = lcpValueNow;
                                vectTriple[minLCPsucText[whichPile[( int )sucSymbol]]].setLcpSucN( minLCPsuc[whichPile[( int )sucSymbol]] ); //I can set lcpSucN of minLCPsucText[alpha[(int)[sucSymbol]]
                                //It closes the LCP interval, so
                                minLCPsucToFind[whichPile[( int )sucSymbol]] = 0; //Close the LCP interval for the symbol sucSymbol
                                minLCPsuc[whichPile[( int )sucSymbol]] = maxValueLen;
                                minLCPsucText[whichPile[( int )sucSymbol]] = 0;
                            }
                        }

                        //Now, I have to update the lcpSucc of the opened interval lcp
                        for ( dataTypedimAlpha gg = 0 ; gg < alphabetSize; gg++ )
                        {
                            //Update the minLCPsuc
                            if ( minLCPsucToFind[( int )gg] == 1 )  //We have to compute the minimum for the lcp of the symbol gg
                            {
                                //if (minLCPsuc[gg] != maxValueLen) {      //There are an LCP interval opened for the symbol c_g
                                if ( lcpValueNow < minLCPsuc[( int )gg] ) //comparison with the last inserted lcp
                                {
                                    minLCPsuc[( int )gg] = lcpValueNow;
                                }
                            }
                            //Update the minLCPCur
                            if ( minLCPcurFound[gg] == 1 ) //I already met the symbol gg. So, I must compute the minimum
                                if ( lcpValueNow < minLCPcur[gg] ) //comparison with the last inserted lcp for the symbol gg
                                {
                                    minLCPcur[gg] = lcpValueNow;
                                }
                        }

                        //For the current LCP
                        minLCPcur[whichPile[( int )sucSymbol]] = maxValueLen; //For each occurrence of sucSymbol, I have to set the min1(the interval starts from the next symbol)
                        minLCPcurFound[whichPile[( int )sucSymbol]] = 1; //So I open the LCP interval for sucSymbol (for the first occurrence of sucSymbol)

                        //We have read another symbol of bwt and its associated lcp
                        cont++;    //number of read symbols
                        //toRead--;
                    }
                    else    //Then there are not other symbols.
                    {
                        //it means that the file does not contain other symbols and we have inserted the symbol in the last position
                        //all lcp intervals have to be close and initializate
                        for ( dataTypedimAlpha gg = 0 ; gg < alphabetSize; gg++ )
                        {
                            if ( minLCPsucToFind[( int )gg] == 1 )  //We have to close the lcp interval of the symbol gg
                            {
                                //if (minLCPsuc[gg] != maxValueLen) {      //LCP interval is apened for the symbol c_g

                                //I have to set the lcpSuc of the text minLCPsucText[(int)gg] to 0
                                vectTriple[minLCPsucText[( int )gg]].setLcpSucN( 0 );
                                minLCPsucToFind[( int )gg] = 0;
                                minLCPsuc[( int )gg] = maxValueLen;
                                minLCPsucText[( int )gg] = 0;
                            }
                        }
                    }
                }
            }
            //}
            k++;   //  I changed the number of the sequence. New iteration.
        }
        //it means that posN[k]<>currentPile, so I have to change BWT-file
        //But before, I have to copy the remainder symbols from the old BWT to new BWT
        //We could need to compute the minLCPsuc for some text
        //		numchar = 1;                   //***********************************
        //if we have inserted the new symbol and we don't read or read the successive symbol, numchar=1,
        //if we have inserted the new symbol and the successive symbol does not exists, numchar=0
        while ( numchar != 0 )
        {
            //For BWT
            numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, InFileBWT );
            numcharWrite = fwrite ( buffer, sizeof( uchar ), numchar , OutFileBWT );
            assert( numchar == numcharWrite ); // we should always read/write the same number of characters
            //For LCP
            numchar = fread( bufferLCP, sizeof( dataTypelenSeq ), SIZEBUFFER, InFileLCP );
            numcharWrite = fwrite ( bufferLCP, sizeof( dataTypelenSeq ), numchar , OutFileLCP );
            assert( numchar == numcharWrite ); // we should always read/write the same number of characters
            //Compute lcpSucN for the other texts
            //For each symbol in the buffer, we check it it close any interval, while each entry in minLcpSuc is maxValue

            //TBD: TO OPTIMIZE. IT CAN END BEFORE. IT DOES NOT NEED TO READ THE ENTIRE BUFFER

            for ( dataTypeNChar bb = 0 ; bb < numchar; bb++ )
            {
                //First, I check if the symbol bb closes a previous interval or it is in the middle or no.
                //In any case, for each symbol, we have to update the minimum
                for ( dataTypedimAlpha gg = 0 ; gg < alphabetSize; gg++ )
                {
                    if ( minLCPsucToFind[( int )gg] == 1 )  //We have to compute the minimum for the lcp of the symbol gg
                    {
                        //if (minLCPsuc[gg] != maxValueLen) {      //There are an LCP interval apened for the symbol c_g
                        if ( bufferLCP[bb] < minLCPsuc[( int )gg] ) //comparison with the last inserted lcp
                        {
                            minLCPsuc[( int )gg] = bufferLCP[bb];
                        }
                    }
                }
                //If the symbol buffer[bb] does not close an LCP interval, ok!
                if ( minLCPsucToFind[whichPile[( int )buffer[bb]]] == 1 )  //We have to compute the minimum for the lcp of the symbol gg
                {
                    //if (minLCPsuc[alpha[(int)buffer[bb]]] != maxValueLen) {       //it means that there is an interval lcp to close for this symbol
                    //The symbol bb closes a LCP interval
                    //I have already computed the minimum (close the previous lcp interval).
                    //I can set lcpSucN of minLCPsucText[alpha[(int)[bb]]
                    vectTriple[minLCPsucText[whichPile[( int )buffer[bb]]]].setLcpSucN( minLCPsuc[whichPile[( int )buffer[bb]]] );
                    //It close the LCP interval, so
                    minLCPsucToFind[whichPile[( int )buffer[bb]]] = 0;
                    minLCPsuc[whichPile[( int )buffer[bb]]] = maxValueLen;
                    minLCPsucText[whichPile[( int )buffer[bb]]] = 0;
                }
            }  //~For
            //           } //~If       //***********************************
        }   //~While
        //The file is finished! but some interval lcp could be opened
        //In this case, we have to set the lcpSuc to 0
        for ( dataTypedimAlpha gg = 0 ; gg < alphabetSize; gg++ )
        {
            if ( minLCPsucToFind[( int )gg] == 1 )  //We have to close the lcp interval of the symbol gg
            {
                //if (minLCPsuc[gg] != maxValueLen) {      //There are an LCP interval opened for the symbol c_g
                vectTriple[minLCPsucText[( int )gg]].setLcpSucN( 0 );
                minLCPsucToFind[( int )gg] = 0;
                minLCPsuc[( int )gg] = maxValueLen;
                minLCPsucText[( int )gg] = 0;
            }
        }

        fclose( InFileBWT );
        fclose( OutFileBWT );
        fclose( InFileLCP );
        fclose( OutFileLCP );

        //Rename files
        //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
        OutFileBWT = fopen( filenameOut, "r" );
        if ( OutFileBWT != NULL ) //If it exists
        {
            fclose( OutFileBWT );
            if ( remove( filenameIn ) != 0 )
                std::cerr << filenameIn << ": Error deleting file" << std::endl;
            else if ( rename( filenameOut, filenameIn ) )
                std::cerr << filenameOut << ": Error renaming " << std::endl;
        }
        //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
        OutFileLCP = fopen( filenameOutLCP, "r" );
        if ( OutFileLCP != NULL ) //If it exists
        {
            fclose( OutFileLCP );
            if ( remove( filenameInLCP ) != 0 )
                std::cerr << filenameInLCP << ": Error deleting file" << std::endl;
            else if ( rename( filenameOutLCP, filenameInLCP ) )
                std::cerr << filenameOutLCP << ": Error renaming " << std::endl;
        }

        j = k;
    }

    //static FILE *tmpFile;
    //tmpFile = fopen("sizeFileBwt.txt", "a");

    if ( verboseEncode == 1 )
    {
        std::cerr << "After the computation of LCP for the next iteration" << std::endl;
        std::cerr << "Q  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            std::cerr << ( int )vectTriple[g].pileN << " ";
        }
        std::cerr << std::endl;
        std::cerr << "P  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].posN  << " ";
        }
        std::cerr << std::endl;
        std::cerr << "N  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].seqN  << " ";
        }
        std::cerr << std::endl;
        std::cerr << "C  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].getLcpCurN()  << " ";
        }
        std::cerr << std::endl;
        std::cerr << "S  ";
        for ( dataTypeNSeq g = 0 ; g < nText; g++ )
        {
            std::cerr << vectTriple[g].getLcpSucN()  << " ";
        }
        std::cerr << std::endl;
    }

    //Renaming new to old
    for ( dataTypedimAlpha g = 0 ; g < alphabetSize; g++ )
    {
        Filename filenameOut( "new_", g );
        Filename filenameIn( g );
        //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
        OutFileBWT = fopen( filenameOut, "rb" );
        if ( OutFileBWT != NULL ) //If it exists
        {
            fclose( OutFileBWT );
            if ( remove( filenameIn ) != 0 )
                std::cerr << filenameIn << ": Error deleting file" << std::endl;
            else if ( rename( filenameOut, filenameIn ) )
                std::cerr << filenameOut << ": Error renaming " << std::endl;
        }

        Filename filenameOutLCP( "new_", g, ".lcp" );
        Filename filenameInLCP( "", g, ".lcp" );
        //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
        OutFileLCP = fopen( filenameOutLCP, "rb" );
        if ( OutFileLCP != NULL ) //If it exists
        {
            fclose( OutFileLCP );
            if ( remove( filenameInLCP ) != 0 )
                std::cerr << filenameInLCP << ": Error deleting file" << std::endl;
            else if ( rename( filenameOutLCP, filenameInLCP ) )
                std::cerr << filenameOutLCP << ": Error renaming " << std::endl;
        }
    }

    delete [] buffer;
    delete [] bufferLCP;
}

void BCRexternalBWT::storeEntireLCP( const char *fn )
{
    assert( false && "TODO" );
}

void BCRexternalBWT::pauseBetweenCyclesIfNeeded()
{
    if ( bwtParams_->getValue( BWT_OPTION_PAUSE_BETWEEN_CYCLES ) == PAUSE_BETWEEN_CYCLES_ON )
    {
        char c;
        fflush( 0 );
        clog << "Iteration complete" << endl;
        clog << " Press a letter + Return to continue" << endl;
        cin >> c;
    }
}
