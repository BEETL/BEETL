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

    for ( unsigned int i=0; i<ramFiles.size(); ++i )
    {
        char filename[100];
        sprintf( filename, "ramdump%d", i );
        ofstream os( filename );
        for ( unsigned int j=0; j<ramFiles[i].size(); j+=2 )
        {
            os << ( ( ( unsigned int )ramFiles[i][j] )&0xFF ) << ", " << ( ( unsigned int )( ramFiles[i][j+1] )&0xFF ) << endl;
        }
    }
}

void debugRamFile( char *filenameIn, size_t n, char *filenameOut = "tmp.debug" )
{
    BwtReaderIncrementalRunLength *pReader;
    pReader = new BwtReaderIncrementalRunLength( filenameIn );
    assert( pReader!=NULL );

    clog << "Writing " << filenameOut << endl;
    BwtWriterASCII *pWriter = new BwtWriterASCII( filenameOut );
    for ( size_t i=0; i<n; ++i )
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
    assert( pReader!=NULL );

    BwtWriterBase *pWriter = instantiateBwtWriterForLastCycle( filenameOut );
    assert( pWriter!=NULL );

    clog << "Writing " << filenameOut << endl;
    while ( pReader->readAndSend( *pWriter, 1 ) ) {}
    clog << " ... done" << endl;

    delete pReader;
    delete pWriter;
}

void BCRexternalBWT::ReadFilesForCycle( const char *prefix, const dataTypelenSeq cycle, const dataTypeNSeq nText, uchar *newSymb, const bool processQualities, uchar *newQual )
{
    char filename[100];
    sprintf ( filename, "%s%u.txt", prefix, cycle );
    FILE *InFileInputText = fopen( filename, "rb" );
    if ( InFileInputText==NULL )
    {
        std::cerr << filename <<" : Error opening " << std::endl;
        exit ( EXIT_FAILURE );
    }
    size_t num = fread( newSymb,sizeof( uchar ),nText,InFileInputText );
    checkIfEqual( num,nText ); // we should always read the same number of characters
    fclose( InFileInputText );
    if ( processQualities )
    {
        char qualFilename[100];
        sprintf ( qualFilename, "%squal.%u.txt", prefix, cycle );
        FILE *InQualFileInputText = fopen( qualFilename, "rb" );
        if ( InQualFileInputText==NULL )
        {
            std::cerr << "buildBCR: " << qualFilename <<" : Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        size_t numQual = fread( newQual,sizeof( uchar ),nText,InQualFileInputText );
        checkIfEqual( numQual,nText );
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
    int intermediateFormat = bwtParams_->getValue( OPTION_INTERMEDIATE_FORMAT );
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
                char debugFilename[1000];
                sprintf ( debugFilename, "%d.%s.afterDefrag.debug", debugCycle, filenameIn );
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
    int intermediateFormat = bwtParams_->getValue( OPTION_INTERMEDIATE_FORMAT );
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
    int outputFormat = bwtParams_->getValue( OPTION_OUTPUT_FORMAT );
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
    if ( bwtParams->getValue( OPTION_PARALLEL_PROCESSING ) != PARALLEL_PROCESSING_OFF )
    {
        // Use nested openmp parallelisation
        omp_set_nested( 1 );
    }
#endif //ifdef USE_OPENMP

    string cycFilesPrefix;
    TransposeFasta transp;
    if ( bwtParams->getValue( OPTION_INPUT_FORMAT ) == INPUT_FORMAT_CYC )
    {
        cycFilesPrefix = string( file1 );
        transp.inputCycFile( cycFilesPrefix );
    }
    else
    {
        cycFilesPrefix = string( fileOut );
        SeqReaderFile *pReader( SeqReaderFile::getReader( fopen( file1,"rb" ) ) );
        transp.init( pReader, bwtParams->getValue( OPTION_PERMUTE_QUALITIES ) == PERMUTE_QUALITIES_ON );
        transp.convert( file1, cycFilesPrefix );
        delete pReader;
    }

    nText = transp.nSeq;
    lengthRead = transp.lengthRead;
    lengthTot = transp.lengthTexts;
    bool processQualities = transp.hasProcessedQualities();

    //#ifdef REPLACE_TABLEOCC
    sizeAlpha=0;
    for ( dataTypedimAlpha i = 0; i < 255; ++i )
        if ( transp.freq[i] > 0 )
        {
            alpha[i] = sizeAlpha;
            sizeAlpha++;
        }
    lengthTot_plus_eof = lengthTot+nText;

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
    char *filename = new char[cycFilesPrefix.length()+sizeof( dataTypelenSeq )*8];
    char *qualFilename = new char[cycFilesPrefix.length()+sizeof( dataTypelenSeq )*8+5];

    Logger::out( LOG_SHOW_IF_VERBOSE ) << "Partial File name for input: " << cycFilesPrefix <<" \n\n";

    static FILE *InFileInputText;

    uchar *newSymb = new uchar[nText];
    uchar *newQual = processQualities?( new uchar[nText] ):NULL;
    uchar *nextSymb = new uchar[nText];
    uchar *nextQual = processQualities?( new uchar[nText] ):NULL;
    vectTriple.resize( nText );

#ifdef REPLACE_TABLEOCC
    tableOcc = new dataTypeNChar*[sizeAlpha];
    for ( dataTypedimAlpha j = 0 ; j < sizeAlpha; j++ )   //Counting for each pile: $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
    {
        tableOcc[j] = new dataTypeNChar[sizeAlpha];
    }
    for ( dataTypedimAlpha j = 0 ; j < sizeAlpha; j++ )
        for ( dataTypedimAlpha h = 0 ; h < sizeAlpha; h++ )
            tableOcc[j][h]=0;
#endif
    tableOcc_.clear();

    //lengthTot = 0;  //Counts the number of symbols

    Logger::out( LOG_SHOW_IF_VERBOSE ) << "\nFirst symbols: "<< "j= "<< 0 <<" - symbols in position " << lengthRead << "\n";

    ReadFilesForCycle( cycFilesPrefix.c_str(), lengthRead-1, nText, newSymb, processQualities, newQual );
    //lengthTot += num;   //Increment the number of chars
    /*
      std::cerr << filename  << std::endl;
      std::cerr << "number read " << num << "\n";
      for (dataTypeNSeq j = 0 ; j < nText; j++)
      std::cerr << newSymb[j];
      std::cerr << ".\n";
    */
    InsertFirstsymbols( newSymb, newQual );

    if ( lengthRead-2 > 0 )
    {
        Logger::out( LOG_SHOW_IF_VERBOSE ) << "Reading next cycle files, time now: " << timer.timeNow();
        Logger::out( LOG_SHOW_IF_VERBOSE ) << "Reading next cycle files, usage: " << timer << endl;
        ReadFilesForCycle( cycFilesPrefix.c_str(), lengthRead-2, nText, newSymb, processQualities, newQual );
    }


    //maxLengthRead-2
    for ( dataTypelenSeq t = lengthRead-2 ; t > 0; t-- )  //dataTypelenSeq is unsigned
    {
        //        if (verboseEncode == 1)
        Logger::out( LOG_SHOW_IF_VERBOSE ) << "j= "<< ( int )( lengthRead - t - 1 ) <<" - symbols in position " << ( int )t << std::endl;
        debugCycle = lengthRead - t - 1;
        Logger::out( LOG_ALWAYS_SHOW ) << "Starting iteration " << t << ", time now: " << timer.timeNow();
        Logger::out( LOG_ALWAYS_SHOW ) << "Starting iteration " << t << ", usage: " << timer << endl;


        //To insert the symbol from position m-3 to position 1
        //The last inserted symbol is in position i+1 (or it is newSymb[j]),
        //the next symbol (to insert) is in position i


        #pragma omp parallel sections
        {
            #pragma omp section
            {
                ReadFilesForCycle( cycFilesPrefix.c_str(), t-1, nText, nextSymb, processQualities, nextQual );
                Logger::out( LOG_SHOW_IF_VERBOSE ) << "Reading input file, time now: " << timer.timeNow();
            }
            #pragma omp section
            {
                InsertNsymbols( newSymb, t, newQual );
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
    }

    debugCycle = lengthRead - 1;
    //The last inserted symbol was in position 1 (or it is newSymb[j]),
    //the next symbol (to insert) is in position 0
    Logger::out( LOG_SHOW_IF_VERBOSE ) << "\n" << "j= "<< lengthRead-1 <<" - symbols in position " << 0 << "\n";

    //    ReadFilesForCycle( fileOut, 0, nText, newSymb, processQualities, newQual );
    InsertNsymbols( newSymb, 0, newQual );

    //The last inserted symbol was in position 0 (or it is newSymb[j]),
    //the next symbol (to insert) is in position m-1, that is, I have to inserted the symbols $
    Logger::out( LOG_SHOW_IF_VERBOSE ) << "\n" << "j= "<< lengthRead <<" - symbols in position " << lengthRead  << ". Inserting $=" << ( int )TERMINATE_CHAR << "=" << TERMINATE_CHAR << " symbol\n\n";
    Logger::out( LOG_ALWAYS_SHOW ) << "Starting iteration " << 0 << ", time now: " << timer.timeNow();
    Logger::out( LOG_ALWAYS_SHOW ) << "Starting iteration " << 0 << ", usage: " << timer << endl;
    for ( dataTypeNSeq j = 0 ; j < nText; j++ )
    {
        newSymb[j] = '$';
        if ( newQual )
            newQual[j] = 0;
        //std::cerr << newSymb[j];
    }
    //std::cerr << "\n";

    debugCycle = lengthRead;
    InsertNsymbols( newSymb, lengthRead, newQual );

    Logger::out( LOG_ALWAYS_SHOW ) << "Final iteration complete, usage: " << timer << endl;


    if ( bwtParams->getValue( OPTION_GENERATE_ENDPOSFILE ) || BUILD_SA )
    {
        Logger::out( LOG_SHOW_IF_VERBOSE ) << "Storing 'outFileEndPos', time now: " << timer.timeNow();
        Logger::out( LOG_SHOW_IF_VERBOSE ) << "Storing 'outFileEndPos', usage: " << timer << endl;

        Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "Stores the 'end positions' of the $!"<< std::endl;
        const char *fileEndPos="outFileEndPos.bwt";
        static FILE *OutFileEndPos;                  // output file of the end positions;
        OutFileEndPos = fopen( fileEndPos, "wb" );
        if ( OutFileEndPos==NULL )
        {
            std::cerr << "Error opening \"" << fileEndPos << "\" file"<< std::endl;
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
        Logger::out( LOG_ALWAYS_SHOW ) << "'end positions' stored!"<< std::endl;
    }

    /* We shouldn't need this anymore as long as we never send the last cycle to ram
      if (bwtParams->getValue( OPTION_INTERMEDIATE_STORAGE_MEDIUM ) == INTERMEDIATE_STORAGE_MEDIUM_RAM)
      {
    #pragma omp parallel for
          for (dataTypedimAlpha i = 1; i < alphabetSize; i++)
          {
              char inputFilename[1000];
              char outputFilename[1000];
              sprintf (inputFilename, "%u", i);
              sprintf (outputFilename, "final.%u", i);
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

    return processQualities?2:1;
} // ~buildBCR

void BCRexternalBWT::InsertFirstsymbols( uchar const *newSymb, uchar const *newSymbQual )
{
    for ( dataTypeNSeq j = 0 ; j < nText; j++ )
    {
        vectTriple[j].posN = j+1;  // position of the suffix (1-based)
        vectTriple[j].seqN = j;   // number of sequence
        vectTriple[j].pileN = 0;    //The first symbols are in $-pile
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
    }


    static FILE *OutFileBWT;                  // output file BWT;
    char *filenameOut = new char[12];
    char *filename = new char[8];
    dataTypeNChar numchar;
    const char *ext = "";
    numchar=sprintf ( filename, "%d",0 );

    numchar=sprintf ( filenameOut,"%s%s",filename,ext );

    OutFileBWT = fopen( filenameOut, "wb" );
    if ( OutFileBWT==NULL )
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
    checkIfEqual( num,nText ); // we should always read the same number of characters

    //if (num != nText)
    // std::cerr << "the written characters is not equal to number of the texts" << num << " and "<< nText <<"\n";
    fclose( OutFileBWT );

    if ( newSymbQual )
    {
        char filenameQualOut[17];
        numchar=sprintf ( filenameQualOut,"%s.qual",filename );
        FILE *OutFileBWTQual = fopen( filenameQualOut, "wb" );
        if ( OutFileBWTQual==NULL )
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
    for ( dataTypedimAlpha i = 1; i <= alphabetSize-1; i++ )
    {
        numchar=sprintf ( filename, "%d", i );
        numchar=sprintf ( filenameOut,"%s%s",filename,ext );
        OutFileBWT = fopen( filenameOut, "wb" );
        if ( OutFileBWT==NULL )
        {
            std::cerr << "BWT file " << ( int )i <<" : Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        fclose( OutFileBWT );

        if ( newSymbQual )
        {
            char filenameQualOut[17];
            numchar=sprintf ( filenameQualOut,"%s.qual",filename );
            FILE *OutFileBWTQual = fopen( filenameQualOut, "wb" );
            if ( OutFileBWTQual==NULL )
            {
                std::cerr << "BWT file $: Error opening " << filenameQualOut << std::endl;
                exit ( EXIT_FAILURE );
            }
            fclose( OutFileBWTQual );
        }
    }

    //Do we want compute the extended suffix array (position and number of sequence)?
    if ( BUILD_SA == 1 )  //To store the SA
    {
        numchar=sprintf ( filename, "sa_%d",0 );
        numchar=sprintf ( filenameOut,"%s%s",filename,ext );
        static FILE  *OutFileSA = fopen( filenameOut, "wb" );    // output file SA;
        if ( OutFileSA==NULL )
        {
            std::cerr << "SA file: Error opening: " << filenameOut << " (SA file $)" << std::endl;
            exit ( EXIT_FAILURE );
        }

        ElementType *newEle = new ElementType[nText];
        for ( dataTypeNSeq j = 0 ; j < nText; j++ )
        {
            //newEle[j].sa=(posSymb + 1) % (lengthRead + 1);
            newEle[j].sa= lengthRead;
            newEle[j].numSeq=j;
            //std::cerr << "(" << (int)newEle[j].sa << ", " << newEle[j].numSeq << ")\n";
        }
        //Store into $-pile SA
        dataTypeNChar num = fwrite ( newEle, sizeof( ElementType ), nText , OutFileSA );
        if ( num != nText )
            std::cerr << "Error: The written characters is not equal to number of the texts in SA" << num << " and "<< nText <<"\n";
        assert( num == nText );

        fclose( OutFileSA );

        //Creates one file for each letter in the alphabet. From 1 to sizeAlpha-1
        for ( dataTypedimAlpha i = 1; i < sizeAlpha; i++ )
        {
            numchar=sprintf ( filename, "sa_%d", i );
            numchar=sprintf ( filenameOut,"%s%s",filename,ext );

            OutFileSA = fopen( filenameOut, "wb" );
            if ( OutFileBWT==NULL )
            {
                std::cerr << "SA file " << ( int )i <<" : Error opening " << std::endl;
                exit ( EXIT_FAILURE );
            }
            fclose( OutFileSA );
        }
    }

    delete [] filenameOut;
    delete [] filename;
}

void BCRexternalBWT::InsertNsymbols( uchar const *newSymb, dataTypelenSeq posSymb, uchar const *newQual )
{
    FILE *InFileBWT;                  // output and input file BWT;
    char *filenameIn = new char[12];
    char *filename = new char[8];
    const char *ext = "";
    dataTypeNChar numchar=0;

    // We first calculate at which index each pile starts
    vector<dataTypeNSeq> pileStarts( alphabetSize,-1 );
    int lastPile = -1;
    for ( dataTypeNSeq j=0; j < nText; ++j )
    {
        int currentPile = vectTriple[j].pileN;
        if ( currentPile == lastPile )
        {
        }
        else
        {
            Logger::out( LOG_FOR_DEBUGGING ) << "pile " << currentPile << " starts at index " << j << endl;
            pileStarts[currentPile] = j;
            lastPile = currentPile;
        }
    }
    Logger::out( LOG_FOR_DEBUGGING ) << "piles finish at index " << nText << endl;
    pileStarts[alphabetSize-1] = nText;
    for ( unsigned int j = alphabetSize-2; j >= 1; --j )
    {
        if ( pileStarts[j] == -1 )
            pileStarts[j] = pileStarts[j+1];
    }
    pileStarts[0] = 0;
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
    for ( parallelPile = 0; parallelPile < alphabetSize-1; ++parallelPile )
    {
        InsertNsymbols_parallelPile( newSymb, posSymb, newQual, parallelPile, pileStarts[parallelPile], pileStarts[parallelPile+1], parallelVectTriplePerNewPile[parallelPile] );
    }


    if ( verboseEncode==1 )
    {
        uchar *buffer = new uchar[SIZEBUFFER];
        dataTypedimAlpha mmm=0;
        while ( mmm < alphabetSize )
        {
            numchar=sprintf ( filename, "%d", mmm );
            numchar=sprintf ( filenameIn,"%s%s",filename,ext );
            //printf("===currentPile= %d\n",mmm);
            InFileBWT = fopen( filenameIn, "r" );
            for ( dataTypeNSeq g = 0 ; g < SIZEBUFFER; g++ )
                buffer[g] = '\0';
            numchar = fread( buffer,sizeof( uchar ),SIZEBUFFER,InFileBWT );
            std::cerr << "B[" << mmm << "]:\t";
            if ( numchar==0 )
                std::cerr  << "empty\n";
            else
                std::cerr  << buffer << "\n";
            fclose( InFileBWT );
            mmm++;
        }
        delete [] buffer;
    } // ~if verboseEncode


#ifdef XXX
    delete [] counters;
#endif
    delete [] filenameIn;
    delete [] filename;

    if ( verboseEncode==1 )
    {
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

    Logger::out( LOG_SHOW_IF_VERBOSE ) << "Before quicksort, time now: " << timer.timeNow();
    Logger::out( LOG_SHOW_IF_VERBOSE ) << "Before quicksort, usage: " << timer << endl;

    //    quickSort(vectTriple);
    vectTriple.clear();
    for ( uint newPile = 0; newPile < alphabetSize-1; ++newPile )
    {
        for ( uint prevPile = 0; prevPile < alphabetSize-1; ++prevPile )
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
    } // ~if verboseEncode

    Logger::out( LOG_SHOW_IF_VERBOSE ) << "After quicksort, time now: " << timer.timeNow();
    Logger::out( LOG_SHOW_IF_VERBOSE ) << "After quicksort, usage: " << timer << endl;

    storeBWT( newSymb,newQual );
    //std::cerr << "End storing BWT" << std::endl;

    //Do we want to compute the generalized suffix array (position and number of sequence)?
    if ( BUILD_SA == 1 )
    {
        storeSA( posSymb );
    }

    //  delete pReader;
}

void BCRexternalBWT::InsertNsymbols_parallelPile( uchar const *newSymb, dataTypelenSeq posSymb, uchar const *newQual, unsigned int parallelPile, dataTypeNSeq startIndex, dataTypeNSeq endIndex, vector< vector< sortElement > > &newVectTriplePerNewPile )
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
    FILE *InFileBWT;                  // output and input file BWT;
    char *filenameIn = new char[12];
    char *filename = new char[8];
    const char *ext = "";

    //std::cerr << "Compute new posN" << std::endl;

    LetterCount counters;

    dataTypeNChar numchar=0;
    dataTypeNChar toRead = 0;
    //Find the positions of the new symbols
    dataTypeNSeq j = startIndex;

    BwtReaderBase *pReader( NULL );
    sortElement newVectTripleItem;

    counters.clear();

    if ( parallelPile == 6 )
    {
        assert( startIndex==endIndex );
        return;
    }
    if ( parallelPile == 0 && startIndex==endIndex )
    {
        return;
    }
    //    clog << "---------- " << (int)vectTriple[j].pileN << "|" << parallelPile << endl;
    dataTypedimAlpha currentPile = parallelPile; //vectTriple[j].pileN;
    assert ( currentPile<alphabetSize-1 );
    numchar=sprintf ( filename, "%d", currentPile );
    numchar=sprintf ( filenameIn,"%s%s",filename,ext );
    //printf("===Current BWT-partial= %d\n",currentPile);

    //#define DUMP_EACH_CYCLE
#ifdef DUMP_EACH_CYCLE
    if ( bwtParams_->getValue( OPTION_INTERMEDIATE_STORAGE_MEDIUM ) == INTERMEDIATE_STORAGE_MEDIUM_RAM )
    {
        char debugFilename[1000];
        sprintf ( debugFilename, "%d.%s.debug", debugCycle, filenameIn );
        dumpRamFileToFile( filenameIn, debugFilename );
    }
#endif //ifdef DUMP_EACH_CYCLE

    if ( posSymb==( lengthRead-2 ) )
        pReader = instantiateBwtReaderForFirstCycle( filenameIn );
    else
        pReader = instantiateBwtReaderForIntermediateCycle( filenameIn );
    assert( pReader!=NULL );

    if ( j < endIndex )
    {
        //    BwtReader reader(filenameIn);

        dataTypeNSeq k=j;
        //For each pile, we have a different counter of characters
        for ( dataTypedimAlpha i = 0 ; i < alphabetSize; i++ )
            counters.count_[i]=0;
        dataTypeNChar cont = 0;   //number of the read symbols
        uchar foundSymbol;
        dataTypeNChar numberRead=0;
        while ( ( k< endIndex ) && ( vectTriple[k].pileN == currentPile ) )
        {
            if ( verboseEncode == 1 )
                std::cerr << "j-1: Q["<<k<<"]=" << ( int )vectTriple[k].pileN << " P["<<k<<"]=" << ( dataTypeNChar )vectTriple[k].posN << " N["<<k<<"]=" << ( dataTypeNSeq )vectTriple[k].seqN << "\t";

            //std::cerr << "--k= " << k << " pileN[k]= " << pileN[k] << " posN[k]= " << posN[k] <<  " seqN[k]= " << seqN[k] << std::endl;
            //For any character (of differents sequences) in the same pile
            foundSymbol = '\0';

            //cont is the number of symbols already read!
            //      toRead = vectTriple[k].posN - cont -1;



            toRead = vectTriple[k].posN - cont;
            //      cout << "toRead: " << toRead << endl;
            if ( toRead>0 )
            {
                if ( toRead>1 )
                {
                    numberRead = ( *pReader ).readAndCount( counters, toRead-1 );
                    assert ( toRead-1 == numberRead );
                }
                assert( ( *pReader )( ( char * )&foundSymbol, 1 )==1 );
                if ( whichPile[( int )foundSymbol]<alphabetSize-1 ) {}
                else
                {
                    cout << ( int )foundSymbol << " " << foundSymbol << endl;
                    assert( 1==0 );
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
                std::cerr << "\nInit New P["<< k <<"]= " << vectTriple[k].posN <<std::endl; //TODO: update this to newVectTripleItem

            for ( dataTypedimAlpha g = 0 ; g < currentPile; g++ )  //I have to count in each pile g= 0... (currentPile-1)-pile
            {
                //                vectTriple[k].posN = vectTriple[k].posN + tableOcc_[g].count_[whichPile[(int)foundSymbol]];
                newVectTripleItem.posN += tableOcc_[g].count_[whichPile[( int )foundSymbol]];
                //std::cerr << "--New posN[k]=" << (int)posN[k] << " tableOcc[g][whichPile[(int)symbol]] " << tableOcc[g][whichPile[(int)symbol]] <<std::endl;
                if ( verboseEncode == 1 )
                {
                    std::cerr << "g= " << ( int )g << " symbol= " << ( int )foundSymbol << " whichPile[symbol]= "
                              << ( int )whichPile[( int )foundSymbol] <<std::endl;
                    std::cerr << "Add New posN[k]=" << vectTriple[k].posN << " tableOcc[g][whichPile[(int)symbol]] "
                              << tableOcc_[g].count_[whichPile[( int )foundSymbol]] <<std::endl;
                }

            }
            //I have to insert the new symbol in the symbol-pile
            assert( whichPile[( int )foundSymbol]<alphabetSize-1 );
            //            vectTriple[k].pileN=whichPile[(int)foundSymbol];
            newVectTripleItem.pileN=whichPile[( int )foundSymbol];
            //std::cerr << "New posN[k]=" << (int)posN[k] << " New pileN[k]=" << (int)pileN[k] << std::endl;
            if ( verboseEncode == 1 )
                std::cerr << "j  : Q[q]=" << ( int )vectTriple[k].pileN << " P[q]=" << ( dataTypeNChar )vectTriple[k].posN <<  " N[q]=" << ( dataTypeNSeq )vectTriple[k].seqN << std::endl;

            newVectTripleItem.seqN = vectTriple[k].seqN;
            //            vectTriple[k] = newVectTripleItem;
            newVectTriplePerNewPile[ newVectTripleItem.pileN ].push_back( newVectTripleItem );

            k++;
        }

        //    fclose(InFileBWT);
        delete pReader;
        pReader = NULL;
        j=k;

        assert ( j == endIndex ); // while loop removed, as we now only process one pile here
    } // ~while j
}

void BCRexternalBWT::storeBWT( uchar const *newSymb, uchar const *newQual )
{

    //I have found the position where I have to insert the chars in the position t of the each text
    //Now I have to update the BWT in each file.

    // We first calculate at which index each pile starts
    vector<dataTypeNSeq> pileStarts( alphabetSize );
    int lastPile = -1;
    for ( dataTypeNSeq j=0; j < nText; ++j )
    {
        int currentPile = vectTriple[j].pileN;
        if ( currentPile == lastPile )
        {
        }
        else
        {
            Logger::out( LOG_FOR_DEBUGGING ) << "pile " << currentPile << " starts at index " << j << endl;
            pileStarts[currentPile] = j;
            lastPile = currentPile;
        }
    }
    Logger::out( LOG_FOR_DEBUGGING ) << "piles finish at index " << nText << endl;
    pileStarts[alphabetSize-1] = nText;
    for ( unsigned int j = alphabetSize-2; j >= 1; --j )
    {
        if ( pileStarts[j] == 0 )
            pileStarts[j] = pileStarts[j+1];
    }


    int parallelPile;
    //    for (int parallelPile = alphabetSize-2; parallelPile >= 0; --parallelPile)
    #pragma omp parallel for
    for ( parallelPile = 0; parallelPile < alphabetSize-1; ++parallelPile )
    {
        storeBWT_parallelPile( newSymb, newQual, parallelPile, pileStarts[parallelPile], pileStarts[parallelPile+1] );
    }


    //static FILE *tmpFile;
    //tmpFile = fopen("sizeFileBwt.txt", "a");
    static FILE *OutFileBWT;                  // output and input file BWT;
    char *filenameOut = new char[16];
    char *filenameIn = new char[12];
    char *filename = new char[8];
    char *filenameQualOut = new char[21];
    char *filenameQualIn = new char[17];
    const char *ext = "";
    dataTypeNChar numchar=0;

    //Renaming new to old
    for ( dataTypedimAlpha g = 0 ; g < alphabetSize-1; g++ )
    {
        numchar=sprintf ( filename, "%d", g );
        numchar=sprintf ( filenameIn,"%s%s",filename,ext );
        numchar=sprintf ( filenameOut,"new_%s%s",filename,ext );
        //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
        OutFileBWT = fopen( filenameOut, "rb" );

        if ( OutFileBWT!=NULL ) //If it exists
        {
            fclose( OutFileBWT );
            if ( remove( filenameIn )!=0 )
                std::cerr << filenameIn <<": Error deleting file" << std::endl;
            else if ( rename( filenameOut,filenameIn ) )
                std::cerr << filenameOut <<": Error renaming " << std::endl;
        }

        if ( verboseEncode == 1 )
        {
            struct stat results;
            if ( stat( filenameIn, &results ) == 0 )
                // The size of the file in bytes is in results.st_size
                //fprintf(tmpFile,"%s\t%u\n", filenameIn, results.st_size);
                std::cerr << filenameIn <<"\t" << results.st_size << std::endl;
            else
                //fprintf(tmpFile,"An error occurred %s\n", filenameIn);
                std::cerr << "An error occurred" << std::endl;
        }

        if ( newQual )
        {
            sprintf ( filenameQualIn,"%s.qual%s",filename,ext );
            sprintf ( filenameQualOut,"new_%s.qual%s",filename,ext );
            FILE *OutQualFileBWT = fopen( filenameQualOut, "rb" );

            if ( OutQualFileBWT!=NULL ) //If it exists
            {
                fclose( OutQualFileBWT );
                if ( remove( filenameQualIn )!=0 )
                    std::cerr << filenameQualIn <<": Error deleting file" << std::endl;
                else if ( rename( filenameQualOut,filenameQualIn ) )
                    std::cerr << filenameQualOut <<": Error renaming " << std::endl;
            }
        }
    }
    //std::cerr <<  std::endl;
    //fprintf(tmpFile,"\n");
    //fclose(tmpFile);
    //  delete pReader;
    //  delete pWriter;

#ifdef OLD
    delete [] buffer;
#endif
    delete [] filenameIn;
    delete [] filename;
    delete [] filenameOut;
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
    FILE *OutFileBWT;                  // output and input file BWT;
    char *filenameOut = new char[16];
    char *filenameIn = new char[12];
    char *filename = new char[8];
    char *filenameQualOut = new char[21];
    char *filenameQualIn = new char[17];
    const char *ext = "";

    dataTypeNChar numchar=0;
    //  dataTypeNChar numcharWrite=0;
#ifdef OLD
    uchar *buffer = new uchar[SIZEBUFFER];
#endif
    dataTypeNChar toRead = 0;

    dataTypeNSeq j;
    dataTypedimAlpha currentPile;
    uchar symbol='\0';

    BwtReaderBase *pReader( NULL );
    BwtWriterBase *pWriter( NULL );
    BwtReaderBase *pQualReader( NULL );
    BwtWriterBase *pQualWriter( NULL );


    j=startIndex;
    while ( j < endIndex )
    {
        currentPile = vectTriple[j].pileN;
        if ( verboseEncode==1 )
            std::cerr << "index j= " << j << " current BWT segment " << ( int )currentPile << std::endl;

        //std::cerr << "Pile " << (int)currentPile << std::endl;
        numchar=sprintf ( filename, "%d", currentPile );
        numchar=sprintf ( filenameIn,"%s%s",filename,ext );
        numchar=sprintf ( filenameOut,"new_%s%s",filename,ext );

        if ( pReader!=NULL ) delete pReader;
        if ( pWriter!=NULL ) delete pWriter;

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


        if ( newQual )
        {
            if ( pQualReader!=NULL ) delete pQualReader;
            if ( pQualWriter!=NULL ) delete pQualWriter;
            sprintf ( filename, "%d", currentPile );
            sprintf ( filenameQualIn,"%s.qual%s",filename,ext );
            sprintf ( filenameQualOut,"new_%s.qual%s",filename,ext );
            pQualReader = new BwtReaderASCII( filenameQualIn );
            pQualWriter = new BwtWriterASCII( filenameQualOut );
        }


        //For each new symbol in the same pile
        dataTypeNSeq k=j;
        dataTypeNChar cont = 0;
        while ( ( k< nText ) && ( vectTriple[k].pileN == currentPile ) )
        {
            if ( verboseEncode==1 )
                std::cerr << "k= " << k << " Q[k]= " << ( int )vectTriple[k].pileN << " P[k]= " << vectTriple[k].posN << " cont = "<< cont << std::endl;
            //std::cerr << "k= " << k << " pileN[k]= " << pileN[k] << " posN[k]= " << posN[k] << std::endl;
            //So I have to read the k-BWT and I have to count the number of the symbols up to the position posN.
            symbol = '\0';
            //As PosN starts to the position 1 and I have to insert the new symbol in position posN[k]
            // I have to read posN[k]-1 symbols
            //cont is the number of symbols already read!
            toRead = ( vectTriple[k].posN-1 ) - cont;
            if ( verboseEncode == 1 )
                std::cerr << "Start: to Read " << toRead << "\n";

            ( *pReader ).readAndSend( ( *pWriter ), toRead );
            if ( newQual )
                ( *pQualReader ).readAndSend( ( *pQualWriter ), toRead );

            cont+=toRead;

            //      toRead=0;


            //Now I have to insert the new symbol associated with the suffix of the sequence k
            //And I have to update the number of occurrences of each symbol
            //   if (toRead==0) {

            ( *pWriter )( ( char * )&newSymb[vectTriple[k].seqN], 1 );
            if ( newQual )
                ( *pQualWriter )( ( char * )&newQual[vectTriple[k].seqN], 1 );

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
        if ( newQual )
            ( *pQualReader ).readAndSend( ( *pQualWriter ) );

        j=k;
    }

    delete pReader;
    pReader=NULL; // %%%
    delete pWriter;
    pWriter=NULL; // %%%
    if ( newQual )
    {
        delete pQualReader;
        pQualReader=NULL;
        delete pQualWriter;
        pQualWriter=NULL;
    }
}

void BCRexternalBWT::storeEntireBWT( const char *fn )
{

    static FILE *OutFileBWT, *InFileBWT;                  // output and input file BWT;
    char *filenameIn = new char[12];
    char *filename = new char[8];
    const char *ext = "";
    dataTypeNChar numchar=0;
    dataTypeNChar numcharWrite=0;

    uchar *buffer = new uchar[SIZEBUFFER];

    dataTypeNChar *freqOut = new dataTypeNChar [256];
    for ( unsigned i = 0; i < 255; ++i )
        freqOut[i] = 0;

    OutFileBWT = fopen( fn, "wb" );
    if ( OutFileBWT==NULL )
    {
        std::cerr << "storeEntireBWT: Error opening " << std::endl;
        exit ( EXIT_FAILURE );
    }

    if ( verboseEncode==1 )
    {
        std::cerr << "\nThe last BWT-segment:"<< std::endl;
        unsigned int mmm=0;
        while ( mmm < alphabetSize )
        {
            numchar=sprintf ( filename, "%d", mmm );
            numchar=sprintf ( filenameIn,"%s%s",filename,ext );
            //printf("===Current BWT-partial= %d\n",mmm);
            InFileBWT = fopen( filenameIn, "rb" );
            for ( dataTypeNChar g = 0 ; g < SIZEBUFFER; g++ )
                buffer[g] = '\0';
            numchar = fread( buffer,sizeof( uchar ),SIZEBUFFER,InFileBWT );
            std::cerr << "B[" << ( int )mmm << "]:\t";
            if ( numchar==0 )
                std::cerr  << "empty";
            else
                std::cerr  << buffer;
            while ( numchar!=0 )
            {
                for ( dataTypeNChar g = 0 ; g < SIZEBUFFER; g++ )
                    buffer[g] = '\0';
                numchar = fread( buffer,sizeof( uchar ),SIZEBUFFER,InFileBWT );
                if ( numchar!=0 )
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
        numchar=sprintf ( filename, "%d", g );
        numchar=sprintf ( filenameIn,"%s%s",filename,ext );
        InFileBWT = fopen( filenameIn, "rb" );
        if ( InFileBWT==NULL )
        {
            std::cerr << "storeEntireBWT: " <<"BWT file " << ( int )g <<": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        //std::cerr << "BWT file " << (int)g << "= ";
        while ( numchar!=0 )
        {
            numchar = fread( buffer,sizeof( uchar ),SIZEBUFFER,InFileBWT );
            //std::cerr << "number read " << numchar << "\n";
            numcharWrite = fwrite ( buffer, sizeof( uchar ), numchar , OutFileBWT );
            checkIfEqual( numchar,numcharWrite ); // we should always read/write the same number of characters

            for ( unsigned j = 0 ; j < numchar; j++ )
                freqOut[( int )( buffer[j] )]++;
        }

        fclose( InFileBWT );
    }
    fclose( OutFileBWT );

    if ( verboseEncode==1 )
    {
        std::cerr << "\nThe Entire BWT:"<< std::endl;
        OutFileBWT = fopen( fn, "rb" );
        if ( OutFileBWT==NULL )
        {
            std::cerr << "storeEntireBWT: Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }

        for ( dataTypeNSeq g = 0 ; g < SIZEBUFFER; g++ )
            buffer[g] = '\0';
        numchar = fread( buffer,sizeof( uchar ),SIZEBUFFER,OutFileBWT );
        if ( numchar==0 )
            std::cerr  << "empty\n";
        else
            std::cerr  << buffer << "\n";
        fclose( OutFileBWT );
    }
    delete [] buffer;



    // freqOut[256]=0;
    std::cerr << "Distribution in BWT\n";
    for ( dataTypedimAlpha i = 0; i < 255; ++i )
        if ( freqOut[i] > 0 )
            std::cerr << i << " " << freqOut[i] << "\n";
    delete [] freqOut;

    delete [] filenameIn;
    delete [] filename;
}

void BCRexternalBWT::storeSA( dataTypelenSeq posSymb )
{

    //I have found the position where I have to insert the chars in the position t of the each text
    //Now I have to update the SA in each file.
    static FILE *OutFileSA, *InFileSA;                  // output and input file SA;
    char *filenameOut = new char[16];
    char *filenameIn = new char[12];
    char *filename = new char[8];
    const char *ext = "";

    dataTypeNChar numchar=0;
    dataTypeNChar numcharWrite=0;
    ElementType *buffer = new ElementType[SIZEBUFFER];
    dataTypeNChar toRead = 0;

    dataTypeNSeq j;
    dataTypedimAlpha currentPile;
    // uchar symbol='\0';
    j=0;
    while ( j < nText )
    {
        currentPile = vectTriple[j].pileN;
        //if (verboseEncode==1)
        // std::cerr << "\nNew Segment; index text j= " << j << " current SA segment is " << (int)currentPile << std::endl;
        //std::cerr << "Pile " << (int)currentPile << std::endl;
        numchar=sprintf ( filename, "sa_%d", currentPile );
        numchar=sprintf ( filenameIn,"%s%s",filename,ext );
        InFileSA = fopen( filenameIn, "rb" );
        if ( InFileSA==NULL )
        {
            std::cerr << "In SA file " << ( int )j <<": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }

        numchar=sprintf ( filenameOut,"new_%s%s",filename,ext );
        OutFileSA = fopen( filenameOut, "wb" );
        if ( OutFileSA==NULL )
        {
            std::cerr << "Out SA file " << ( int )j <<": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        //std::cerr << "In File " << filenameIn << std::endl;
        //std::cerr << "Out File " << filenameOut << std::endl;

        //For each new symbol in the same pile
        dataTypeNSeq k=j;
        dataTypeNChar cont = 0;
        while ( ( k< nText ) && ( vectTriple[k].pileN == currentPile ) )
        {

            //if (verboseEncode==1)
            // std::cerr << "k= " << k << " Q[k]= " << (int)vectTriple[k].pileN << " P[k]= " << vectTriple[k].posN << " cont = "<< cont << std::endl;
            //std::cerr << "k= " << k << " pileN[k]= " << pileN[k] << " posN[k]= " << posN[k] << std::endl;
            //So I have to read the k-SA and I have to count the number of the symbols up to the position posN.
            //As PosN starts to the position 1 and I have to insert the new symbol in position posN[k]
            // I have to read posN[k]-1 symbols
            //cont is the number of symbols already read!
            toRead = ( vectTriple[k].posN-1 ) - cont;
            /*
            if (verboseEncode == 1)
                std::cerr << "Start: to Read " << toRead << "\n";
            */
            while ( toRead > 0 )            //((numchar!=0) && (toRead > 0)) {
            {
                if ( toRead < SIZEBUFFER ) //The last reading for this sequence
                {
                    numchar = fread( buffer,sizeof( ElementType ),toRead,InFileSA );
                    /*
                    if (verboseEncode == 1)
                        std::cerr << "number read " << numchar << " to Read " << toRead << "\n";
                    */
                    checkIfEqual( numchar,toRead ); // we should always read/write the same number of characters

                    numcharWrite = fwrite ( buffer, sizeof( ElementType ), numchar , OutFileSA );
                    checkIfEqual( numchar,numcharWrite ); // we should always read/write the same number of characters
                    //std::cerr << "toread number write " << numcharWrite << "\n";
                }
                else
                {
                    numchar = fread( buffer,sizeof( ElementType ),SIZEBUFFER,InFileSA );
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
            if ( toRead==0 )
            {
                ElementType newEle;
                newEle.sa=( posSymb + 1 ) % ( lengthRead+1 );
                newEle.numSeq=vectTriple[k].seqN;

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
        while ( numchar!=0 )
        {
            numchar = fread( buffer,sizeof( ElementType ),SIZEBUFFER,InFileSA );
            //std::cerr << "After insert: " << numchar << "\n";
            numcharWrite = fwrite ( buffer, sizeof( ElementType ), numchar , OutFileSA );
            checkIfEqual( numchar, numcharWrite ); // we should always read/write the same number of characters
        }

        fclose( InFileSA );
        fclose( OutFileSA );
        j=k;
    }

    //Renaming new to old
    for ( dataTypedimAlpha g = 0 ; g < alphabetSize; g++ )
    {
        numchar=sprintf ( filename, "sa_%d", g );
        numchar=sprintf ( filenameIn,"%s%s",filename,ext );
        numchar=sprintf ( filenameOut,"new_%s%s",filename,ext );
        //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
        OutFileSA = fopen( filenameOut, "rb" );

        if ( OutFileSA!=NULL ) //If it exists
        {
            fclose( OutFileSA );
            if ( remove( filenameIn )!=0 )
                std::cerr << filenameIn <<": Error deleting file" << std::endl;
            else if ( rename( filenameOut,filenameIn ) )
                std::cerr << filenameOut <<": Error renaming " << std::endl;
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
    delete [] filenameIn;
    delete [] filename;
    delete [] filenameOut;
}

void BCRexternalBWT::storeEntirePairSA( const char *fn )
{

    std::cerr << "\nEntire Pairs SA file (position, number of sequence)" << std::endl;

    static FILE *OutFileSA, *InFileSA;                  // output and input file SA;
    char *filenameIn = new char[12];
    char *filename = new char[8];
    const char *ext = "";
    dataTypeNChar numcharWrite, numcharRead;
    ElementType *buffer = new ElementType[SIZEBUFFER];

    int lung = strlen( fn );
    char *fnSA = new char[lung+7];
    numcharRead=sprintf ( fnSA,"%s%s",fn,".pairSA" );

    OutFileSA = fopen( fnSA, "wb" );
    if ( OutFileSA==NULL )
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
        numcharRead=sprintf ( filename, "sa_%d", g );
        numcharRead=sprintf ( filenameIn,"%s%s",filename,ext );
        InFileSA = fopen( filenameIn, "rb" );
        if ( InFileSA==NULL )
        {
            std::cerr << "SA file " << ( int )g <<": Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }
        //std::cerr << "SA file " << (int)g << "= ";

        numcharRead = fread( buffer,sizeof( ElementType ),SIZEBUFFER,InFileSA );

        /* //it will be useful for varying length reads
        //Correction of the length of the sequences.
        for (dataTypeNSeq num = 0; num < numcharRead; num++) {
            buffer[num].sa = buffer[num].sa - (lengthRead - vectLen[buffer[num].numSeq]);
        }
        */

        numcharWrite = fwrite ( buffer, sizeof( ElementType ), numcharRead , OutFileSA );
        checkIfEqual ( numcharRead , numcharWrite );

        while ( numcharRead!=0 )
        {
            numcharRead = fread( buffer,sizeof( ElementType ),SIZEBUFFER,InFileSA );
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
        if ( remove( filenameIn )!=0 )
            std::cerr << filenameIn <<": Error deleting file" << std::endl;
    }

    fclose( OutFileSA );

    if ( verboseEncode==1 )
    {
        OutFileSA = fopen( fnSA, "rb" );
        if ( OutFileSA==NULL )
        {
            std::cerr << "Entire SA file: Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }

        numcharRead = fread( buffer,sizeof( ElementType ),SIZEBUFFER,OutFileSA );
        if ( numcharRead==0 )
            std::cerr  << "empty\n";
        else
            for ( dataTypeNChar g = 0 ; g < numcharRead; g++ )
            {
                std::cerr  << "(" << ( int )buffer[g].sa << "," << buffer[g].numSeq << ") ";
            }
        while ( numcharRead!=0 )
        {
            numcharRead = fread( buffer,sizeof( ElementType ),SIZEBUFFER,OutFileSA );
            for ( dataTypeNChar g = 0 ; g < numcharRead; g++ )
            {
                std::cerr  << "(" << buffer[g].sa << "," << buffer[g].numSeq << ") ";
            }
        }
        std::cerr << std::endl;

        fclose( OutFileSA );
    }

    delete [] buffer;
    delete [] filenameIn;
    delete [] filename;
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
    int lung = strlen( fn );
    char *fnSA = new char[lung+3];
    char *fnPairSA = new char[lung+7];
    numchar=sprintf ( fnSA,"%s%s",fn,".sa" );
    numchar=sprintf ( fnPairSA,"%s%s",fn,".pairSA" );

    InFilePairSA = fopen( fnPairSA, "rb" );
    if ( InFilePairSA==NULL )
    {
        std::cerr << "Entire Pairs SA file: Error opening " << fnPairSA << std::endl;
        exit ( EXIT_FAILURE );
    }

    OutFileSA = fopen( fnSA, "wb" );
    if ( OutFileSA==NULL )
    {
        std::cerr << "Entire SA file: Error opening " << fnSA << std::endl;
        exit ( EXIT_FAILURE );
    }

    ElementType *buffer = new ElementType[SIZEBUFFER];
    dataTypeNChar *bufferNChar = new dataTypeNChar[SIZEBUFFER];

    while ( !feof( InFilePairSA ) )
    {
        numchar = fread( buffer,sizeof( ElementType ),SIZEBUFFER,InFilePairSA );
        //std::cerr << "number read " << numchar << "\n";
        if ( numchar > 0 )
        {
            for ( dataTypeNChar i=0; i < numchar; i++ )
            {
                bufferNChar[i] = ( dataTypeNChar )( buffer[i].numSeq * ( lengthRead+1 ) + buffer[i].sa );
                //std::cerr << buffer[i].numSeq << " " << (int)lengthRead << " " << (int)buffer[i].sa << " --> " << (int)bufferNChar[i] << "\n";
                //bufferNChar[i] = (dataTypeNChar)(vectSumCumLen[buffer[i].numSeq] + buffer[i].sa);       //it will be useful for varying length reads
                //std::cerr << "vectSumCumLen["<< buffer[i].numSeq<< "]= " << (int)vectSumCumLen[buffer[i].numSeq] << " + " << (int)buffer[i].sa << " --> " << (int)bufferNChar[i] << "\n";
            }
            numcharWrite = fwrite ( bufferNChar, sizeof( dataTypeNChar ), numchar, OutFileSA );
            //std::cerr << "number write " << numcharWrite << "\n";
        }
    }
    fclose( InFilePairSA );
    fclose( OutFileSA );

    if ( verboseEncode==1 )
    {
        std::cerr << "\nThe Entire SA. The file is "<< fnSA << std::endl;
        OutFileSA = fopen( fnSA, "rb" );
        if ( OutFileSA==NULL )
        {
            std::cerr << "Entire SA file: Error opening " << std::endl;
            exit ( EXIT_FAILURE );
        }

        numchar = fread( bufferNChar,sizeof( dataTypeNChar ),SIZEBUFFER,OutFileSA );
        if ( numchar==0 )
            std::cerr  << "empty\n";
        else
            for ( dataTypeNChar g = 0 ; g < numchar; g++ )
            {
                std::cerr  << bufferNChar[g] << " ";
            }
        while ( numchar!=0 )
        {
            numchar = fread( buffer,sizeof( ElementType ),SIZEBUFFER,OutFileSA );
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
    delete [] fnSA;
    delete [] fnPairSA;
}
