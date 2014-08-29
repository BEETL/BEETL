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
#include "BwtReader.hh"
#include "BwtWriter.hh"
#include "Filename.hh"
#include "LetterCount.hh"
#include "PredictiveEncoding.hh"
#include "SeqReader.hh"
#include "Timer.hh"
#include "Tools.hh"
#include "TransposeFasta.hh"
#include "libzoo/util/Logger.hh"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>

#ifdef _OPENMP
# include <omp.h>
#endif //ifdef _OPENMP

using namespace std;
using namespace BeetlBwtParameters;


bool SAPstopped = false;
vector<SequenceNumber> sapCount;
const int sizeAlphaM1 = 4;
vector<char> whichPileSAP;


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
        ofstream os( filename.str() );
        for ( unsigned int j = 0; j < ramFiles[i].size(); j += 2 )
        {
            os << ( ( ( unsigned int )ramFiles[i][j] ) & 0xFF ) << ", " << ( ( unsigned int )( ramFiles[i][j + 1] ) & 0xFF ) << endl;
        }
    }
}

void debugRamFile( const string &filenameIn, size_t n, const string &filenameOut = "tmp.debug" )
{
    BwtReaderIncrementalRunLength pReader( filenameIn );

    Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "Writing " << filenameOut << endl;
    BwtWriterASCII pWriter( filenameOut );
    for ( size_t i = 0; i < n; ++i )
    {
        pReader.readAndSend( pWriter, 1 );
    }
    fflush( 0 );

    while ( 1 )
    {
        pReader.readAndSend( pWriter, 1 );
    }
}

void BCRexternalBWT::convertFileFromIntermediateToFinalFormat( const char *filenameIn, const char *filenameOut )
{
    BwtReaderBase *pReader = instantiateBwtReaderForIntermediateCycle( filenameIn );
    assert( pReader != NULL );

    BwtWriterBase *pWriter = instantiateBwtWriterForLastCycle( filenameOut );
    assert( pWriter != NULL );

    Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "Writing " << filenameOut << endl;
    pReader->readAndSend( *pWriter );
    Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << " ... done" << endl;

    delete pReader;
    delete pWriter;
}

void BCRexternalBWT::ReadFilesForCycle( const char *prefix, const SequenceLength cycle1, const SequenceLength readLength, const SequenceNumber nText, uchar *newSymb, const bool processQualities, uchar *newQual )
{
    SequenceNumber count( nText );
    if ( bwtParams_->getValue( PARAMETER_ADD_REV_COMP ) == 1 )
    {
        count /= 2;
    }

    for ( int revComp = 0; revComp < 2; ++revComp )
    {
        SequenceLength cycle( cycle1 );
        if ( revComp == 1 )
        {
            if ( bwtParams_->getValue( PARAMETER_ADD_REV_COMP ) == 0 )
                break;
            // This takes care of adding reversed reads
            cycle = readLength - 1 - cycle1;
        }

        Filename filename( prefix, cycle, "" );
        Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Opening file " << filename << " for reading" << endl;
        FILE *InFileInputText = fopen( filename, "rb" );
        if ( InFileInputText == NULL )
        {
            cerr << filename << " : Error opening " << endl;
            exit ( EXIT_FAILURE );
        }
        //#define ALL_A_FOR_DEBUGGING
#ifdef ALL_A_FOR_DEBUGGING
        for ( SequenceNumber i = 0; i < count; ++i )
            newSymb[i + ( revComp ? count : 0 )] = 'A';
#else
        SequenceNumber num = fread( newSymb + ( revComp ? count : 0 ), sizeof( uchar ), count, InFileInputText );
        checkIfEqual( num, count ); // we should always read the same number of characters
#endif
        fclose( InFileInputText );

        if ( revComp == 1 )
        {
            // This takes care of changing bases to their complement
            uchar *startBase = newSymb + count;
            uchar *endBase = startBase + count;
            for ( uchar *ptr = startBase; ptr < endBase; ++ptr )
            {
                // TODO: Use EAGLE's special class to do this comp operation, after moving it to libzoo
                switch ( *ptr )
                {
                    case 'A':
                    case 'a':
                        *ptr = 'T';
                        break;
                    case 'C':
                    case 'c':
                        *ptr = 'G';
                        break;
                    case 'G':
                    case 'g':
                        *ptr = 'C';
                        break;
                    case 'T':
                    case 't':
                        *ptr = 'A';
                        break;
                }
            }
        }

        if ( processQualities )
        {
            Filename qualFilename( prefix, "qual.", cycle, "" );
            Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Opening file " << qualFilename << " for reading" << endl;
            FILE *InQualFileInputText = fopen( qualFilename, "rb" );
            if ( InQualFileInputText == NULL )
            {
                cerr << "buildBCR: " << qualFilename << " : Error opening " << endl;
                exit ( EXIT_FAILURE );
            }
            size_t numQual = fread( newQual + ( revComp ? count : 0 ), sizeof( uchar ), count, InQualFileInputText );
            checkIfEqual( numQual, count );
            fclose( InQualFileInputText );
        }
    }
}

BwtReaderBase *BCRexternalBWT::instantiateBwtReaderForIntermediateCycle( const char *filenameIn, bool allowDefrag )
{
    BwtReaderBase *pReader = NULL;
    int intermediateFormat = bwtParams_->getValue( PARAMETER_INTERMEDIATE_FORMAT );
    if ( intermediateFormat == INTERMEDIATE_FORMAT_MULTIRLE && filenameIn[strlen( filenameIn ) - 1] == '0' )
        intermediateFormat = INTERMEDIATE_FORMAT_RLE;
    switch ( intermediateFormat )
    {
        case INTERMEDIATE_FORMAT_ASCII:
            pReader = new BwtReaderASCII( filenameIn );
            break;
        case INTERMEDIATE_FORMAT_RLE:
            pReader = new BwtReaderRunLengthV3( filenameIn );
            break;
        case INTERMEDIATE_FORMAT_MULTIRLE:
            pReader = new BwtReaderIncrementalRunLength( filenameIn );
            break;
#ifdef ACTIVATE_HUFFMAN
        case INTERMEDIATE_FORMAT_HUFFMAN:
            pReader = new BwtReaderHuffman( filenameIn );
            break;
#endif
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

            Logger_if( LOG_SHOW_IF_VERBOSE )
            {
                #pragma omp critical
                {
                    Logger::out() << "After defrag, time now: " << timer.timeNow();
                }
            }

            delete( pReader );
            pReader = new BwtReaderIncrementalRunLength( filenameIn );


#ifdef DUMP_EACH_CYCLE
            {
                Filename debugFilename( "", debugCycle, string( filenameIn ) + ".afterDefrag.debug" );
                convertFileFromIntermediateToFinalFormat( filenameIn, debugFilename );
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
    int intermediateFormat = bwtParams_->getValue( PARAMETER_INTERMEDIATE_FORMAT );
    if ( intermediateFormat == INTERMEDIATE_FORMAT_MULTIRLE && filenameOut[strlen( filenameOut ) - 1] == '0' )
        intermediateFormat = INTERMEDIATE_FORMAT_RLE;
    switch ( intermediateFormat )
    {
        case INTERMEDIATE_FORMAT_ASCII:
            pWriter = new BwtWriterASCII( filenameOut );
            break;
        case INTERMEDIATE_FORMAT_RLE:
            pWriter = new BwtWriterRunLengthV3( filenameOut );
            break;
        case INTERMEDIATE_FORMAT_MULTIRLE:
            pWriter = new BwtWriterIncrementalRunLength( filenameOut );
            break;
#ifdef ACTIVATE_HUFFMAN
        case INTERMEDIATE_FORMAT_HUFFMAN:
            pWriter = new BwtWriterHuffman( filenameOut );
            break;
#endif
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
    int outputFormat = bwtParams_->getValue( PARAMETER_OUTPUT_FORMAT );
    switch ( outputFormat )
    {
        case OUTPUT_FORMAT_ASCII:
            pWriter = new BwtWriterASCII( filenameOut );
            break;
        case OUTPUT_FORMAT_RLE:
            pWriter = new BwtWriterRunLengthV3( filenameOut );
            break;
#ifdef ACTIVATE_HUFFMAN
        case OUTPUT_FORMAT_HUFFMAN:
            pWriter = new BwtWriterHuffman( filenameOut );
            break;
#endif
        default:
            cerr << "Error in BCRexternalBWT::instantiateBwtWriterForLastCycle: unknown output format: " << outputFormat << endl;
            exit ( EXIT_FAILURE );
    }
    assert( pWriter );
    return pWriter;
}

BwtReaderBase *BCRexternalBWT::instantiateBwtReaderForLastCycle( const char *filenameOut )
{
    BwtReaderBase *pReader = NULL;
    int outputFormat = bwtParams_->getValue( PARAMETER_OUTPUT_FORMAT );
    switch ( outputFormat )
    {
        case OUTPUT_FORMAT_ASCII:
            pReader = new BwtReaderASCII( filenameOut );
            break;
        case OUTPUT_FORMAT_RLE:
            pReader = new BwtReaderRunLengthV3( filenameOut );
            break;
#ifdef ACTIVATE_HUFFMAN
        case OUTPUT_FORMAT_HUFFMAN:
            pReader = new BwtReaderHuffman( filenameOut );
            break;
#endif
        default:
            cerr << "Error in BCRexternalBWT::instantiateBwtReaderForLastCycle: unknown output format: " << outputFormat << endl;
            exit ( EXIT_FAILURE );
    }
    assert( pReader );
    return pReader;
}


int BCRexternalBWT::buildBCR( const string &file1, const string &fileOut, const BwtParameters *bwtParams )
{
#ifdef _OPENMP
    //    if ( bwtParams->getValue( PARAMETER_PARALLEL_PROCESSING ) != PARALLEL_PROCESSING_OFF )
    {
        // Use nested openmp parallelisation
        omp_set_nested( 1 );
    }
#endif //ifdef _OPENMP

    const bool permuteQualities = ( bwtParams_->getValue( PARAMETER_PROCESS_QUALITIES ) == PROCESS_QUALITIES_PERMUTE );
    const bool generateCycleQualities = ( bwtParams_->getValue( PARAMETER_GENERATE_CYCLE_QUAL ) != GENERATE_CYCLE_QUAL_OFF );
    const bool readQualities = permuteQualities || generateCycleQualities;

    string cycFilesPrefix;
    TransposeFasta transp;
    if ( bwtParams->getValue( PARAMETER_INPUT_FORMAT ) == INPUT_FORMAT_CYC )
    {
        cycFilesPrefix = string( file1 );
        transp.inputCycFile( cycFilesPrefix );
    }
    else
    {
        TmpFilename cycFilesPrefix2( fileOut ); // "move" filename to temp directory
        cycFilesPrefix = cycFilesPrefix2.str();
        FILE *f;
        if ( file1 == "-" )
            f = stdin;
        else
            f = fopen( file1.c_str(), "rb" );
        SeqReaderFile *pReader( SeqReaderFile::getReader( f ) );
        transp.init( pReader, readQualities );
        transp.convert( cycFilesPrefix );
        delete pReader;
        if ( f != stdin )
            fclose( f );
    }

    nText = transp.nSeq;
    lengthRead = transp.lengthRead;
    lengthTot = transp.lengthTexts;
    bool processQualities = transp.hasProcessedQualities();

    SequenceLength currentIteration = 0; // iteration counter of symbols insertion
    SequenceLength currentCycleFileNum; // file num processed in the current iteration (may not be equal to file num being read, in case of prefetch)
    int cycleFileNumIncrement;
    if ( bwtParams->getValue( PARAMETER_REVERSE ) == true )
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
    for ( AlphabetSymbol i = 0; i < 255; ++i )
        if ( transp.freq[i] > 0 )
        {
            alpha[i] = sizeAlpha;
            sizeAlpha++;
        }
    lengthTot_plus_eof = lengthTot + nText;

    Logger_if( LOG_FOR_DEBUGGING )
    {
        Logger::out() << "We supposed that the symbols in the input file are:\n";
        for ( AlphabetSymbol i = 0; i < 255; ++i )
            if ( transp.freq[i] > 0 )
                Logger::out() << i << " " << transp.freq[i] << " " << ( int )alpha[i] << "\n";
        //#endif

        Logger::out() << "sizeof(type size of alpha): " << sizeof( AlphabetSymbol ) << "\n";
        Logger::out() << "sizeof(type of #sequences): " << sizeof( SequenceNumber ) << "\n";
        Logger::out() << "sizeof(type of #characters): " << sizeof( LetterNumber ) << "\n";

        Logger::out() << "\nalphabetSize: " << ( int )alphabetSize << "\n";
        Logger::out() << "Number of sequences: " << nText << "\n";
        Logger::out() << "Length of each sequence: " << lengthRead << "\n\n";
        Logger::out() << "Total length (without $): " << lengthTot << "\n";
        Logger::out() << "Total length (with $): " << lengthTot_plus_eof << "\n";
    }

    Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Partial File name for input: " << cycFilesPrefix << " \n\n";

    // Make sure that our temporary filenames are not already used
    for ( AlphabetSymbol i = 0 ; i < sizeAlpha; ++i )
    {
        for ( int pass = 0; pass < 2; ++pass )
        {
            TmpFilename filename( pass ? "" : "new_", i, "" );
            ofstream os( filename.str() );
            os.close();
            if ( unlink( filename ) != 0 )
            {
                cerr << "Error preparing temporary file " << filename << endl;
                exit ( -1 );
            }
        }
    }

    // Reverse-complemented reads
    if ( bwtParams_->getValue( PARAMETER_ADD_REV_COMP ) == 1 )
    {
        SequenceNumber newNText = nText * 2;
        if (newNText < nText)
        {
            Logger::error() << "Error: Too many sequences. This version of BEETL was compiled for a maximum of " << maxSequenceNumber << " sequences, but this input has " << (static_cast<uint64_t>(nText)*2) << " sequences. You can increase this limit by changing the type definition of 'SequenceNumber' in Types.hh and recompiling BEETL." << endl;
            exit( -1 );
        }
        nText = newNText;
    }

    // Prepare reset point between multiple reads
    int nextIterationReset = -1;
    if ( ( *bwtParams_ )[ PARAMETER_SUB_SEQUENCE_LENGTH ].isSet() )
        nextIterationReset = ( *bwtParams_ )[ PARAMETER_SUB_SEQUENCE_LENGTH ];

    // Real start
    uchar *newSymb = new uchar[nText];
    uchar *newQual = processQualities ? ( new uchar[nText] ) : NULL;
    uchar *nextSymb = new uchar[nText];
    uchar *nextQual = processQualities ? ( new uchar[nText] ) : NULL;
    vectTriple.resize( nText );

#ifdef REPLACE_TABLEOCC
    tableOcc = new LetterNumber*[sizeAlpha];
    for ( AlphabetSymbol j = 0 ; j < sizeAlpha; j++ )   //Counting for each pile: $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
    {
        tableOcc[j] = new LetterNumber[sizeAlpha];
    }
    for ( AlphabetSymbol j = 0 ; j < sizeAlpha; j++ )
        for ( AlphabetSymbol h = 0 ; h < sizeAlpha; h++ )
            tableOcc[j][h] = 0;
#endif
    tableOcc_.clear();

    Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "\nFirst symbols: " << "Iteration " << 0 << " - symbols in (zero-based) position " << currentCycleFileNum << "\n";
    Logger::out() << "Starting iteration " << currentIteration << ", time now: " << timer.timeNow();
    Logger::out() << "Starting iteration " << currentIteration << ", usage: " << timer << endl;

    ReadFilesForCycle( cycFilesPrefix.c_str(), currentCycleFileNum, lengthRead, nText, newSymb, processQualities, newQual );
    InitialiseTmpFiles();


    if ( ( *bwtParams_ )[PARAMETER_SAP_ORDERING] == true )
    {
        uchar *newSymb3 = new uchar[nText];

        whichPileSAP.resize( 256, -1 );
        whichPileSAP['A'] = 0;
        whichPileSAP['C'] = 1;
        whichPileSAP['G'] = 2;
        whichPileSAP['T'] = 3;

        sapCount.clear();
        sapCount.resize( sizeAlphaM1, 0 );

        Logger_if( LOG_FOR_DEBUGGING )
        {
            cerr << "Before1:" << endl;
            for ( SequenceNumber i = 0; i < nText; i++ )
                cerr << "Triple[" << i << "]: " << vectTriple[i].seqN << " " << vectTriple[i].posN << " " << ( int )vectTriple[i].pileN << ", newSymb=" << newSymb[i] << endl;
        }

        //        vector <sortElement> vectTriple2( nText );

        SequenceNumber pos = 0;
        for ( AlphabetSymbol j = 0 ; j < sizeAlphaM1; j++ )
        {
            for ( SequenceNumber i = 0; i < nText; ++i )
            {
                if ( whichPileSAP[( int )newSymb[i]] == -1 )
                {
                    cerr << "Error SAP with char " << newSymb[i] << " at position " << i << endl;
                    assert( false );
                }
                if ( whichPileSAP[( int )newSymb[i]] == j )
                {
                    newSymb3[pos] = newSymb[i];
                    ++pos;
                    ++sapCount[j];
                }
            }
        }
        assert( pos == nText );

        Logger_if( LOG_FOR_DEBUGGING )
        {
            cerr << "After1:" << endl;
            for ( SequenceNumber i = 0; i < nText; i++ )
                cerr << "Triple[" << i << "]: " << vectTriple[i].seqN << " " << vectTriple[i].posN << " " << ( int )vectTriple[i].pileN << ", newSymb=" << newSymb[i] << ", newSymb3=" << newSymb3[i] << endl;
            cerr << "sapCount=";
            for ( unsigned int i = 0; i < sapCount.size(); ++i )
                cerr << sapCount[i] << ",";
            cerr << endl;
        }
        InsertFirstsymbols( newSymb3, newQual );
    }
    else
        InsertFirstsymbols( newSymb, newQual );



    // Update iteration counters
    ++currentIteration;
    currentCycleFileNum += cycleFileNumIncrement;
    debugCycle = currentIteration;

    if ( lengthRead >= 2 )
    {
        Logger_if( LOG_SHOW_IF_VERBOSE )
        {
            Logger::out() << "Reading next cycle files, time now: " << timer.timeNow();
            Logger::out() << "Reading next cycle files, usage: " << timer << endl;
        }
        ReadFilesForCycle( cycFilesPrefix.c_str(), currentCycleFileNum, lengthRead, nText, nextSymb, processQualities, nextQual );


        if ( ( *bwtParams_ )[PARAMETER_SAP_ORDERING] == true )
        {
            Logger_if( LOG_FOR_DEBUGGING )
            {
                cerr << "Before:" << endl;
                for ( SequenceNumber i = 0; i < nText; i++ )
                    cerr << "Triple[" << i << "]: " << vectTriple[i].seqN << " " << vectTriple[i].posN << " " << ( int )vectTriple[i].pileN << ", newSymb=" << newSymb[i] << ", nextSymb=" << nextSymb[i] << endl;
            }

            vector <sortElement> vectTriple2( nText );
            uchar *nextSymb2 = new uchar[nText];

            SequenceNumber pos = 0;
            for ( AlphabetSymbol j = 0 ; j < sizeAlpha; j++ )
            {
                for ( SequenceNumber i = 0; i < nText; ++i )
                {
                    if ( whichPile[( int )newSymb[i]] == j )
                    {
                        nextSymb2[pos] = nextSymb[i];
                        vectTriple2[pos] = vectTriple[i];
                        vectTriple2[pos].posN = pos + 1;
                        ++pos;
                    }
                }
            }
            assert( pos == nText );
            /*
                        {
                            uchar *tmp;
                            tmp = nextSymb;
                            nextSymb = nextSymb2;
                            nextSymb2 = tmp;
                        }
            */
            vectTriple.swap( vectTriple2 );
            delete [] nextSymb2;

            Logger_if( LOG_FOR_DEBUGGING )
            {
                cerr << "After:" << endl;
                for ( SequenceNumber i = 0; i < nText; i++ )
                    cerr << "Triple[" << i << "]: " << vectTriple[i].seqN << " " << vectTriple[i].posN << " " << ( int )vectTriple[i].pileN << ", newSymb=" << newSymb[i] << ", nextSymb=" << nextSymb[i] << endl;
            }
        }

        {
            uchar *tmp;
            tmp = newSymb;
            newSymb = nextSymb;
            nextSymb = tmp;

            tmp = newQual;
            newQual = nextQual;
            nextQual = tmp;
        }

    }

    while ( currentIteration <= lengthRead - 2 )
        //    for ( SequenceLength t = lengthRead - 2 ; t > 0; t-- ) //SequenceLength is unsigned
    {
        pauseBetweenCyclesIfNeeded ();

        Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Iteration " << ( int ) currentIteration << " - symbols in position " << ( int ) currentCycleFileNum << endl;
        Logger::out() << "Starting iteration " << currentIteration << ", time now: " << timer.timeNow();
        Logger::out() << "Starting iteration " << currentIteration << ", usage: " << timer << endl;







        //To insert the symbol from position m-3 to position 1
        //The last inserted symbol is in position i+1 (or it is newSymb[j]),
        //the next symbol (to insert) is in position i

        #pragma omp parallel sections
        {
            #pragma omp section
            {
                if ( ( int )currentIteration == nextIterationReset )
                {
                    Logger::out() << "Resetting counters" << endl;
                    nextIterationReset += ( *bwtParams_ )[ PARAMETER_SUB_SEQUENCE_LENGTH ];

                    Logger_if( LOG_SHOW_IF_VERBOSE )
                    {
                        Logger::out() << "Inserting new '$' symbols" << endl;
                    }
                    // Standalone block
                    {
                        vector<uchar> newSymb2( nText, '$' );
                        assert( newSymb2.size() == nText );
                        assert( newSymb2[nText - 1] == '$' );
                        InsertNsymbols( newSymb2.data(), currentIteration, NULL );
                    }

                    if ( bwtParams->getValue( PARAMETER_GENERATE_ENDPOSFILE ) || BUILD_SA )
                        writeEndPosFile( 1, false );

                    Logger_if( LOG_SHOW_IF_VERBOSE )
                    {
                        Logger::out() << "Continuing iteration " << currentIteration << ", time now: " << timer.timeNow();
                        Logger::out() << "Continuing iteration " << currentIteration << ", usage: " << timer << endl;
                    }
                    int subSequenceCount = currentIteration / ( *bwtParams_ )[ PARAMETER_SUB_SEQUENCE_LENGTH ];
                    InsertFirstsymbols( newSymb, newQual, subSequenceCount );
                }
                else
                {
                    InsertNsymbols( newSymb, currentIteration, newQual );
                }
            }
            #pragma omp section
            {
                ReadFilesForCycle( cycFilesPrefix.c_str(), currentCycleFileNum + cycleFileNumIncrement, lengthRead, nText, nextSymb, processQualities, nextQual );
                Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Finished reading input file, time now: " << timer.timeNow();
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
    Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Iteration " << ( int ) currentIteration << " - symbols in position " << ( int ) currentCycleFileNum << endl;
    Logger::out() << "Starting iteration " << currentIteration << ", time now: " << timer.timeNow();
    Logger::out() << "Starting iteration " << currentIteration << ", usage: " << timer << endl;
    assert( currentIteration == lengthRead - 1 );
    assert( currentCycleFileNum == 0 || currentCycleFileNum == lengthRead - 1 ); // depending on the --reverse flag
    InsertNsymbols( newSymb, currentIteration, newQual );
    // Update iteration counters
    ++currentIteration;
    debugCycle = currentIteration;

    pauseBetweenCyclesIfNeeded ();

    //The last inserted symbol was in position 0 (or it is newSymb[j]),
    //the next symbol (to insert) is in position m-1, that is, I have to inserted the symbols $
    Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Final iteration " << ( int ) currentIteration << " - Inserting $=" << ( int )terminatorChar << "=" << terminatorChar << " symbols" << endl;
    Logger::out() << "Starting iteration " << currentIteration << ", time now: " << timer.timeNow();
    Logger::out() << "Starting iteration " << currentIteration << ", usage: " << timer << endl;
    assert( currentIteration == lengthRead );
    for ( SequenceNumber j = 0 ; j < nText; j++ )
    {
        newSymb[j] = '$';
        if ( processQualities )
            newQual[j] = 0;
    }
    InsertNsymbols( newSymb, currentIteration, newQual );

    Logger::out() << "Final iteration complete, time now: " << timer.timeNow();
    Logger::out() << "Final iteration complete, usage: " << timer << endl;


    // We don't need this one anymore
    delete pWriterBwt0_;
    pWriterBwt0_ = 0;

    if ( bwtParams->getValue( PARAMETER_GENERATE_ENDPOSFILE ) || BUILD_SA )
        writeEndPosFile( 0, true );

    /* We shouldn't need this anymore as long as we never send the last cycle to ram
      if (bwtParams->getValue( PARAMETER_INTERMEDIATE_STORAGE_MEDIUM ) == INTERMEDIATE_STORAGE_MEDIUM_RAM)
      {
    #pragma omp parallel for
          for (AlphabetSymbol i = 1; i < alphabetSize; i++)
          {
              char inputFilename[1000];
              char outputFilename[1000];
              TmpFilename inputFilename( i );
              Filename outputFilename( "final.", i );
              convertFileFromIntermediateToFinalFormat( outputCompression_, inputFilename, compressionRunLength, outputFilename);
          }
          cout << "Disk output complete, usage: " << timer << endl;
      }
    */

    // to delete those:
    delete [] newSymb;
    delete [] nextSymb;
    // vectTriple.~vector<sortElement>();
    /*
      delete [] seqN;
      delete [] pileN;
      delete [] posN;
    */
    //  cerr << endl;
    //cerr << "The input text is long " << lengthTot << endl;

    LetterNumber numCharInTable = 0;
    for ( AlphabetSymbol r = 0; r < alphabetSize; r++ )
    {
        for ( AlphabetSymbol t = 0; t < alphabetSize; t++ )
        {
            numCharInTable += tableOcc_[r].count_[t];
        }
    }
    Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "In tableOcc_, there are " << numCharInTable << " letters" << endl;

#ifdef REPLACING_TABLEOCC
    for ( AlphabetSymbol j = 0 ; j < sizeAlpha; j++ )
    {
        delete [] tableOcc[j];
        tableOcc[j] = NULL;
    }
    delete [] tableOcc;
#endif

    if ( bwtParams_->getValue( PARAMETER_GENERATE_LCP ) == true )
    {
        if ( bwtParams_->getValue( PARAMETER_OUTPUT_FORMAT ) != OUTPUT_FORMAT_ASCII )
        {
            // LCP source code only generate ASCII BWT, so we convert it here if needed
            Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Converting ASCII BWT to " << bwtParams_->getStringValue( PARAMETER_OUTPUT_FORMAT ) << endl;
            for ( AlphabetSymbol g = 0 ; g < alphabetSize; g++ )
            {
                TmpFilename filename1( "", g, "" );
                TmpFilename filename2( "new_", g, "" );

                BwtReaderBase *pReader = new BwtReaderASCII( filename1 );
                BwtWriterBase *pWriter = instantiateBwtWriterForLastCycle( filename2 );

                pReader->readAndSend( *pWriter );

                delete pReader;
                delete pWriter;

                if ( remove( filename1 ) != 0 )
                    cerr << "Error deleting file " << filename1 << endl;
                else if ( safeRename( filename2, filename1 ) )
                    cerr << "Error renaming " << filename2 << " to " << filename1 << endl;
            }
        }
    }
    else
    {
        if ( bwtParams_->getStringValue( PARAMETER_OUTPUT_FORMAT ) != bwtParams_->getStringValue( PARAMETER_INTERMEDIATE_FORMAT ) )
        {
            // We only convert BWT pile 0 here, as other piles are converted during the final iteration
            TmpFilename filename1( "", 0, "" );
            TmpFilename filename2( "new_", 0, "" );
            convertFileFromIntermediateToFinalFormat( filename1, filename2 );

            if ( remove( filename1 ) != 0 )
                cerr << "Error deleting file " << filename1 << endl;
            else if ( safeRename( filename2, filename1 ) )
                cerr << "Error renaming " << filename2 << " to " << filename1 << endl;
        }
    }

    return permuteQualities ? 2 : 1;
} // ~buildBCR

void BCRexternalBWT::InitialiseTmpFiles()
{
    //Creates empty files for each letter in the alphabet
    for ( AlphabetSymbol i = 0; i < alphabetSize; ++i )
    {
        TmpFilename filenameOut( i );
        if ( i == 0 ) // BWT0 file is created separately
        {
            pWriterBwt0_ = instantiateBwtWriterForIntermediateCycle( filenameOut );
        }
        else
        {
            unique_ptr<BwtWriterBase> emptyPileFile( instantiateBwtWriterForIntermediateCycle( filenameOut ) );
        }

        const bool permuteQualities = ( bwtParams_->getValue( PARAMETER_PROCESS_QUALITIES ) == PROCESS_QUALITIES_PERMUTE );
        if ( permuteQualities )
        {
            TmpFilename filenameQualOut( "", i, ".qual" );
            FILE *OutFileBWTQual = fopen( filenameQualOut, "wb" );
            if ( OutFileBWTQual == NULL )
            {
                cerr << "BWT file $: Error opening " << filenameQualOut << endl;
                exit ( EXIT_FAILURE );
            }
            fclose( OutFileBWTQual );
        }

        if ( bwtParams_->getValue( PARAMETER_GENERATE_LCP ) == true )
        {
            TmpFilename filenameOut( "", i, ".lcp" );
            FILE *OutFileLCP = fopen( filenameOut, "wb" );
            if ( OutFileLCP == NULL )
            {
                cerr << "LCP file: " << filenameOut << " : Error opening " << endl;
                exit ( EXIT_FAILURE );
            }
            fclose( OutFileLCP );
        }

        //Do we want compute the extended suffix array (position and number of sequence)?
        if ( BUILD_SA == 1 )  //To store the SA
        {
            TmpFilename filenameOut( "sa_", i );
            FILE *OutFileSA = fopen( filenameOut, "wb" );
            if ( OutFileSA == NULL )
            {
                cerr << "SA file " << ( int )i << " : Error opening " << endl;
                exit ( EXIT_FAILURE );
            }
            fclose( OutFileSA );
        }
    }

}

void BCRexternalBWT::InsertFirstsymbols( uchar const *newSymb, uchar const *newSymbQual, const int subSequenceNum )
{
    for ( SequenceNumber j = 0 ; j < nText; j++ )
    {
        vectTriple[j] = sortElement( 0, nText * subSequenceNum + j + 1, j );
    }
    if ( verboseEncode == 1 )
    {
        cerr << "First step" << endl;
        cerr << "Q  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            cerr << ( int )vectTriple[g].pileN << " ";
        }
        cerr << endl;
        cerr << "P  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            cerr << vectTriple[g].posN  << " ";
        }
        cerr << endl;
        cerr << "N  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            cerr << vectTriple[g].seqN  << " ";
        }
        cerr << endl;

        if ( bwtParams_->getValue( PARAMETER_GENERATE_LCP ) == true )
        {
            cerr << "C  ";             //LCP current
            for ( SequenceNumber g = 0 ; g < nText; g++ )
            {
                cerr << ( int )vectTriple[g].getLcpCurN()  << " ";
            }
            cerr << endl;
            cerr << "S  ";                     //LCP successive
            for ( SequenceNumber g = 0 ; g < nText; g++ )
            {
                cerr << ( int )vectTriple[g].getLcpSucN()  << " ";
            }
            cerr << endl;
        }
    }

    for ( SequenceNumber j = 0 ; j < nText; j++ )
    {
        tableOcc_[0].count_[whichPile[( int )newSymb[j]]]++;     //counting the number of occurrences in BWT of the $-pile
    }
    //Store newSymb into $-pile BWT
    ( *pWriterBwt0_ )( ( const char * )newSymb, nText );
    pWriterBwt0_->flush();

    const bool permuteQualities = ( bwtParams_->getValue( PARAMETER_PROCESS_QUALITIES ) == PROCESS_QUALITIES_PERMUTE );
    if ( permuteQualities )
    {
        assert ( newSymbQual );
        TmpFilename filenameQualOut( "0.qual" );
        FILE *OutFileBWTQual = fopen( filenameQualOut, "ab" );
        if ( OutFileBWTQual == NULL )
        {
            cerr << "BWT file $: Error opening " << filenameQualOut << endl;
            exit ( EXIT_FAILURE );
        }
        LetterNumber numQual = fwrite ( newSymbQual, sizeof( uchar ), nText, OutFileBWTQual );
        checkIfEqual( numQual, nText );
        fclose( OutFileBWTQual );
    }


    if ( bwtParams_->getValue( PARAMETER_GENERATE_LCP ) == true )
    {
        SequenceLength *vectLCP = new SequenceLength[nText];
        for ( SequenceNumber j = 0 ; j < nText; j++ )
        {
            vectLCP[j] = 0;
        }
        FILE *OutFileLCP;                  // output and input file LCP;
        TmpFilename filenameOut( "0.lcp" );
        OutFileLCP = fopen( filenameOut, "ab" );
        if ( OutFileLCP == NULL )
        {
            cerr << "LCP file $: " << filenameOut << " Error opening " << endl;
            exit ( EXIT_FAILURE );
        }

        LetterNumber num = fwrite ( vectLCP, sizeof( SequenceLength ), nText , OutFileLCP );
        checkIfEqual( num , nText ); // we should always read the same number of integers
        //vectLCP.clear();
        fclose( OutFileLCP );
        delete [] vectLCP;
    }

    //Do we want compute the extended suffix array (position and number of sequence)?
    if ( BUILD_SA == 1 )  //To store the SA
    {
        TmpFilename filenameOut( "sa_0" );
        FILE *OutFileSA = fopen( filenameOut, "ab" );    // output file SA;
        if ( OutFileSA == NULL )
        {
            cerr << "SA file: Error opening: " << filenameOut << " (SA file $)" << endl;
            exit ( EXIT_FAILURE );
        }

        ElementType *newEle = new ElementType[nText];
        for ( SequenceNumber j = 0 ; j < nText; j++ )
        {
            //newEle[j].sa=(iterationNum + 1) % (lengthRead + 1);
            newEle[j].sa = lengthRead;
            newEle[j].numSeq = j;
            //cerr << "(" << (int)newEle[j].sa << ", " << newEle[j].numSeq << ")\n";
        }
        //Store into $-pile SA
        LetterNumber num = fwrite ( newEle, sizeof( ElementType ), nText , OutFileSA );
        if ( num != nText )
            cerr << "Error: The written characters is not equal to number of the texts in SA" << num << " and " << nText << "\n";
        assert( num == nText );

        fclose( OutFileSA );
        delete [] newEle;
    }
}

void BCRexternalBWT::InsertNsymbols( uchar const *newSymb, SequenceLength iterationNum, uchar const *newQual )
{
    LetterNumber numchar = 0;

    // We first calculate at which index each pile starts
    vector<SequenceNumber> pileStarts( alphabetSize + 1 );
    pileStarts[0] = 0;
    SequenceNumber index = 0;
    for ( int pile = 1; pile < alphabetSize + 1; ++pile )
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
    vector< vector< FragmentedVector< sortElement > > > parallelVectTriplePerNewPile( alphabetSize, vector< FragmentedVector< sortElement > >( alphabetSize ) );

    int parallelPile;
    //    for (int parallelPile = alphabetSize-2; parallelPile >= 0; --parallelPile)
    #pragma omp parallel for
    for ( parallelPile = 0; parallelPile < alphabetSize; ++parallelPile )
    {
        InsertNsymbols_parallelPile( newSymb, iterationNum, newQual, parallelPile, pileStarts[parallelPile], pileStarts[parallelPile + 1], parallelVectTriplePerNewPile[parallelPile] );
    }

    Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Finished inserting symbols in RAM, time now: " << timer.timeNow();

    if ( verboseEncode == 1 )
    {
        cerr << "The segments before inserting are:\n";
        uchar *buffer = new uchar[SIZEBUFFER];
        AlphabetSymbol mmm = 0;
        while ( mmm < alphabetSize )
        {
            TmpFilename filenameIn( mmm );
            //printf("===currentPile= %d\n",mmm);
            FILE *InFileBWT = fopen( filenameIn, "r" );
            for ( SequenceNumber g = 0 ; g < SIZEBUFFER; g++ )
                buffer[g] = '\0';
            numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, InFileBWT );
            cerr << "B[" << mmm << "]:\t";
            if ( numchar == 0 )
                cerr  << "empty\n";
            else
                cerr  << buffer << "\n";
            fclose( InFileBWT );
            mmm++;
        }
        delete [] buffer;
    } // ~if verboseEncode

    if ( verboseEncode == 1 )
    {
        if ( bwtParams_->getValue( PARAMETER_GENERATE_LCP ) == true )
        {
            SequenceLength *bufferLCP = new SequenceLength[SIZEBUFFER];
            AlphabetSymbol mmm = 0;
            while ( mmm < alphabetSize )
            {
                TmpFilename filenameInLCP( "", mmm, ".lcp" );
                FILE *InFileLCP = fopen( filenameInLCP, "rb" );
                for ( LetterNumber g = 0 ; g < SIZEBUFFER; g++ )
                    bufferLCP[g] = 0;
                numchar = fread( bufferLCP, sizeof( SequenceLength ), SIZEBUFFER, InFileLCP );
                cerr << "L[" << ( int )mmm << "]:\t";
                if ( numchar == 0 )
                    cerr  << "empty";
                else
                    for ( SequenceNumber g = 0 ; g < numchar; g++ )
                        cerr  << ( int )bufferLCP[g] << " ";
                cerr  << "\n";
                fclose( InFileLCP );
                mmm++;
            }
            delete [] bufferLCP;
        }

        cerr << "NewSymbols " ;
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            cerr << newSymb[g] << " ";
        }
        cerr << endl;
        cerr << "Before Sorting" << endl;
        cerr << "Q  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            cerr << ( int )vectTriple[g].pileN << " ";
        }
        cerr << endl;
        cerr << "P  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            cerr << vectTriple[g].posN  << " ";
        }
        cerr << endl;
        cerr << "N  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            cerr << vectTriple[g].seqN  << " ";
        }
        cerr << endl;
    } // ~if verboseEncode

#ifdef XXX
    delete [] counters;
#endif

    //    quickSort(vectTriple);
    vectTriple.clear();
    for ( uint newPile = 0; newPile < alphabetSize; ++newPile )
    {
        for ( uint prevPile = 0; prevPile < alphabetSize; ++prevPile )
        {
            if ( !parallelVectTriplePerNewPile[prevPile][newPile].empty() )
            {
                //                vectTriple.insert( vectTriple.end(), parallelVectTriplePerNewPile[prevPile][newPile].begin(), parallelVectTriplePerNewPile[prevPile][newPile].end() );
                parallelVectTriplePerNewPile[prevPile][newPile].appendTo( vectTriple );
            }
        }
    }

    Logger_if( LOG_FOR_DEBUGGING )
    {
        Logger::out() << "After Sorting" << endl;
        Logger::out() << "U  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            Logger::out() << newSymb[g] << " ";
        }
        Logger::out() << endl;
        Logger::out() << "Q  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            Logger::out() << ( int )vectTriple[g].pileN << " ";
        }
        Logger::out() << endl;
        Logger::out() << "P  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            Logger::out() << vectTriple[g].posN  << " ";
        }
        Logger::out() << endl;
        Logger::out() << "N  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            Logger::out() << vectTriple[g].seqN  << " ";
        }
        Logger::out() << endl;

        if ( bwtParams_->getValue( PARAMETER_GENERATE_LCP ) == true )
        {
            Logger::out() << "C  ";
            for ( SequenceNumber g = 0 ; g < nText; g++ )
            {
                Logger::out() << ( int )vectTriple[g].getLcpCurN()  << " ";
            }
            Logger::out() << endl;
            Logger::out() << "S  ";
            for ( SequenceNumber g = 0 ; g < nText; g++ )
            {
                Logger::out() << ( int )vectTriple[g].getLcpSucN()  << " ";
            }
            Logger::out() << endl;
        }
    } // ~if verboseEncode



    if ( ( *bwtParams_ )[PARAMETER_SAP_ORDERING] == true && !SAPstopped )
    {
        Logger_if( LOG_FOR_DEBUGGING )
        {
            Logger::out() << "Before2:" << endl;
            for ( SequenceNumber i = 0; i < nText; i++ )
                Logger::out() << "Triple[" << i << "]: " << vectTriple[i].seqN << " " << vectTriple[i].posN << " " << ( int )vectTriple[i].pileN << ", newSymb=" << newSymb[i] << endl;
        }

        vector <sortElement> vectTriple2( nText );
        //            uchar *nextSymb2 = new uchar[nText];
        //            uchar *newSymb3 = new uchar[nText];

        vector<SequenceNumber> sapAccumulatedCount( sapCount.size(), 0 );
        vector<SequenceNumber> sapCount2;
        try
        {
            sapCount2.resize( sapCount.size() * sizeAlphaM1, 0 );
        }
        catch ( const std::exception &e )
        {
            cerr << "not enough RAM. Stopping SAP" << endl;
            vector<SequenceNumber> emptySapCount;
            sapCount.swap( emptySapCount );
            SAPstopped = true;
        }

        if ( !SAPstopped )
        {
            for ( unsigned int sapSet = 0; sapSet < sapCount.size(); ++sapSet )
                sapAccumulatedCount[sapSet] = ( sapSet > 0 ? sapAccumulatedCount[sapSet - 1] : 0 ) + sapCount[sapSet];

            SequenceNumber sapActiveSetCount = 0;
            SequenceNumber sapActiveBaseCount = 0;
            SequenceNumber sapVeryActiveSetCount = 0;
            SequenceNumber sapVeryActiveBaseCount = 0;
            #pragma omp parallel for
            for ( unsigned int sapSet = 0; sapSet < sapCount.size(); ++sapSet )
            {
                if ( sapCount[sapSet] == 0 )
                    continue;

                if ( sapCount[sapSet] > 1 )
                {
                    #pragma omp critical
                    {
                        ++sapActiveSetCount;
                        sapActiveBaseCount += sapCount[sapSet];
                    }
                }

                Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "sapCount[" << sapSet << "]=" << sapCount[sapSet] << endl;
                SequenceNumber startSeqNum = ( sapSet > 0 ? sapAccumulatedCount[sapSet - 1] : 0 ); //currentSeq;
                SequenceNumber endSeqNum = sapAccumulatedCount[sapSet]; //currentSeq + sapCount[sapSet];
                SequenceNumber pos = startSeqNum;

                for ( AlphabetSymbol j = 0 ; j < sizeAlphaM1; j++ )
                {
                    for ( SequenceNumber i = startSeqNum; i < endSeqNum; ++i )
                    {
                        int s = vectTriple[i].seqN;
                        if ( whichPileSAP[( int )newSymb[s]] == -1 )
                        {
                            cerr << "Error SAP with char " << newSymb[s] << " at position " << s << endl;
                            assert( false );
                        }
                        if ( whichPileSAP[( int )newSymb[s]] == j )
                        {
                            vectTriple2[pos] = vectTriple[i];
                            vectTriple2[pos].posN = vectTriple[pos].posN;
                            vectTriple2[pos].pileN = vectTriple[pos].pileN;
                            ++pos;
                            ++sapCount2[sapSet + j * sapCount.size()];
                        }
                    }
                }

                if ( newSymb[vectTriple2[startSeqNum].seqN] != newSymb[vectTriple2[endSeqNum - 1].seqN] )
                {
                    #pragma omp critical
                    {
                        ++sapVeryActiveSetCount;
                        sapVeryActiveBaseCount += sapCount[sapSet];
                    }
                }

                //            currentSeq += sapCount[sapSet];
                assert( pos == endSeqNum );
            }

            Logger::out() << "SAP active sets=" << sapActiveSetCount << ", active bases=" << sapActiveBaseCount << ", veryActive sets=" << sapVeryActiveSetCount << ", veryActive bases=" << sapVeryActiveBaseCount << endl;
            if ( sapActiveSetCount == 0 )
            {
                SAPstopped = true;
            }
            vectTriple.swap( vectTriple2 );
            sapCount.swap( sapCount2 );

            Logger_if( LOG_FOR_DEBUGGING )
            {
                Logger::out() << "After2:" << endl;
                for ( SequenceNumber i = 0; i < nText; i++ )
                    Logger::out() << "Triple[" << i << "]: " << vectTriple[i].seqN << " " << vectTriple[i].posN << " " << ( int )vectTriple[i].pileN << ", newSymb=" << newSymb[i] << endl;
                Logger::out() << "sapCount=";
                for ( unsigned int i = 0; i < sapCount.size(); ++i )
                {
                    if ( i % sizeAlphaM1 == 0 )
                        Logger::out() << "   ";
                    if ( i && i % ( sizeAlphaM1 * sizeAlphaM1 ) == 0 )
                        Logger::out() << endl << "  ";
                    Logger::out() << sapCount[i] << ",";
                }
                Logger::out() << endl;
            }
        }
    }

    Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Writing intermediate BWT files to disk, time now: " << timer.timeNow();

    if ( bwtParams_->getValue( PARAMETER_GENERATE_LCP ) == false )
    {
        storeBWT( newSymb, newQual );
    }
    else
    {
        storeBWTandLCP( newSymb );
    }

    //cerr << "End storing BWT" << endl;

    //Do we want to compute the generalized suffix array (position and number of sequence)?
    if ( BUILD_SA == 1 )
    {
        storeSA( iterationNum );
    }

    //  delete pReader;
}

void BCRexternalBWT::InsertNsymbols_parallelPile( uchar const *newSymb, SequenceLength iterationNum, uchar const *newQual, unsigned int parallelPile, SequenceNumber startIndex, SequenceNumber endIndex, vector< FragmentedVector< sortElement > > &newVectTriplePerNewPile )
{
    Logger_if( LOG_FOR_DEBUGGING )
    {
        #pragma omp critical
        {
            Logger::out() << "InsertNsymbols_parallelPile: pile=" << parallelPile << " from " << startIndex << " to " << endIndex << endl;
        }
    }

    //<<<<<<< BuildBCR.cpp
    //they are not first symbols
    //cerr << "Compute new posN" << endl;

    LetterCount counters;

    //Find the positions of the new symbols
    SequenceNumber j = startIndex;

    BwtReaderBase *pReader( NULL );

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
    AlphabetSymbol currentPile = parallelPile; //vectTriple[j].pileN;
    assert ( currentPile < alphabetSize );
    TmpFilename filename( "", currentPile, "" );
    //printf("===Current BWT-partial= %d\n",currentPile);

    //#define DUMP_EACH_CYCLE
#ifdef DUMP_EACH_CYCLE
    if ( bwtParams_->getValue( PARAMETER_INTERMEDIATE_STORAGE_MEDIUM ) == INTERMEDIATE_STORAGE_MEDIUM_RAM )
    {
        Filename debugFilename( "", debugCycle, "." + filename.str() + ".debug" );
        convertFileFromIntermediateToFinalFormat( filename, debugFilename );
    }
#endif //ifdef DUMP_EACH_CYCLE

    pReader = instantiateBwtReaderForIntermediateCycle( filename );
    assert( pReader != NULL );

    if ( j < endIndex )
    {
        sortElement newVectTripleItem;
        //    BwtReader reader(filename.str().c_str());

        SequenceNumber k = j;
        //For each pile, we have a different counter of characters
        for ( AlphabetSymbol i = 0 ; i < alphabetSize; i++ )
            counters.count_[i] = 0;
        LetterNumber cont = 0;   //number of the read symbols
        uchar foundSymbol;
        LetterNumber numberRead = 0;
        while ( ( k < endIndex ) && ( vectTriple[k].pileN == currentPile ) )
        {
            //For any character (of differents sequences) in the same pile
            if ( verboseEncode == 1 )
            {
                cerr << "j-1: Q[" << k << "]=" << ( int )vectTriple[k].pileN << " P[" << k << "]=" << ( LetterNumber )vectTriple[k].posN << " N[" << k << "]=" << ( SequenceNumber )vectTriple[k].seqN << "\t";
                //cerr << "--k= " << k << " pileN[k]= " << pileN[k] << " posN[k]= " << posN[k] <<  " seqN[k]= " << seqN[k] << endl;
            }
            foundSymbol = '\0';

            //cont is the number of symbols already read!
            LetterNumber toRead = vectTriple[k].posN - cont;
            if ( toRead > 0 )
            {
                Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "toRead: " << toRead << endl;
                if ( toRead > 1 )
                {
                    numberRead = ( *pReader ).readAndCount( counters, toRead - 1 );
                    if ( toRead - 1 != numberRead )
                    {
                        cerr << "ERROR: toRead - 1 != numberRead" << endl;
                        cerr << "  toRead-1 = " << ( toRead - 1 ) << endl;
                        cerr << "  numberRead=" << numberRead << endl;
                        cerr << "  cont=" << cont << endl;
                        cerr << "  k=" << k << endl;
                        cerr << "  vectTriple[k].posN=" << vectTriple[k].posN << endl;
                        cerr << "  newSymb= " << newSymb[0] << newSymb[1] << newSymb[2] << "..." << endl;
                        cerr << "  iterationNum=" << iterationNum << endl;
                        if ( newQual )
                            cerr << "  newQual=" << newQual[0] << newQual[1] << newQual[2] << "..." << endl;
                        cerr << "  parallelPile=" << parallelPile << endl;
                        cerr << "  startIndex=" << startIndex << endl;
                        cerr << "  endIndex=" << endIndex << endl;
                        cerr << "  j=" << j << endl;
                        assert( false );
                    }
                }
                assert( ( *pReader )( ( char * )&foundSymbol, 1 ) == 1 );
                if ( whichPile[( int )foundSymbol] < alphabetSize ) {}
                else
                {
                    cout << ( int )foundSymbol << " " << foundSymbol << endl;
                    assert( 1 == 0 );
                }

                counters.count_[whichPile[( int )foundSymbol]]++;
                cont += toRead;
            }

            Logger_if( LOG_FOR_DEBUGGING )
            {
                Logger::out() << "toRead=" << toRead << ", foundSymbol=" << foundSymbol  << ", counters=" << counters << endl;
            }

            //cerr << "toRead " << toRead << "Found Symbol is " << foundSymbol << "\n";


            //I have to update the value in vectTriple[k].posN, it must contain the position of the new symbol
            //#ifdef XXX
            //vectTriple[k].posN = counters.count_[whichPile[(int)foundSymbol]];
            newVectTripleItem.posN = counters.count_[whichPile[( int )foundSymbol]];
            //#endif
            //cerr << "--New posN[k]=" << (int)posN[k] <<endl;
            if ( verboseEncode == 1 )
                cerr << "\nInit New P[" << k << "]= " << vectTriple[k].posN << endl; //TODO: update this to newVectTripleItem

            for ( AlphabetSymbol g = 0 ; g < currentPile; g++ )  //I have to count in each pile g= 0... (currentPile-1)-pile
            {
                //                vectTriple[k].posN = vectTriple[k].posN + tableOcc_[g].count_[whichPile[(int)foundSymbol]];
                newVectTripleItem.posN += tableOcc_[g].count_[whichPile[( int )foundSymbol]];
                //cerr << "--New posN[k]=" << (int)posN[k] << " tableOcc[g][whichPile[(int)symbol]] " << tableOcc[g][whichPile[(int)symbol]] <<endl;
                if ( verboseEncode == 1 )
                {
                    cerr << "g= " << ( int )g << " symbol= " << ( int )foundSymbol << " whichPile[symbol]= "
                         << ( int )whichPile[( int )foundSymbol] << endl;
                    cerr << "Add New posN[k]=" << vectTriple[k].posN << " tableOcc[g][whichPile[(int)symbol]] "
                         << tableOcc_[g].count_[whichPile[( int )foundSymbol]] << endl;
                }

            }
            //I have to insert the new symbol in the symbol-pile
            assert( whichPile[( int )foundSymbol] < alphabetSize );
            //            vectTriple[k].pileN=whichPile[(int)foundSymbol];
            newVectTripleItem.pileN = whichPile[( int )foundSymbol];
            //cerr << "New posN[k]=" << (int)posN[k] << " New pileN[k]=" << (int)pileN[k] << endl;
            if ( verboseEncode == 1 )
                cerr << "j  : Q[q]=" << ( int )vectTriple[k].pileN << " P[q]=" << ( LetterNumber )vectTriple[k].posN <<  " N[q]=" << ( SequenceNumber )vectTriple[k].seqN << endl;

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
    vector<SequenceNumber> pileStarts( alphabetSize + 1 );
    pileStarts[0] = 0;
    SequenceNumber index = 0;
    for ( int pile = 1; pile < alphabetSize + 1; ++pile )
    {
        while ( index < nText && vectTriple[index].pileN < pile )
            ++index;
        pileStarts[pile] = index;
    }


    int parallelPile;
    //    for (int parallelPile = alphabetSize-2; parallelPile >= 0; --parallelPile)
    #pragma omp parallel for
    for ( parallelPile = 0; parallelPile < alphabetSize; ++parallelPile )
    {
        storeBWT_parallelPile( newSymb, newQual, parallelPile, pileStarts[parallelPile], pileStarts[parallelPile + 1] );
    }


    //Renaming new to old
    for ( AlphabetSymbol g = 0 ; g < alphabetSize; g++ )
    {
        TmpFilename filenameIn( g );
        TmpFilename filenameOut( "new_", g );
        //cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << endl;
        FILE *OutFileBWT = fopen( filenameOut, "rb" );

        if ( OutFileBWT != NULL ) //If it exists
        {
            fclose( OutFileBWT );
            if ( remove( filenameIn ) != 0 )
                cerr << filenameIn << ": Error deleting file" << endl;
            else if ( safeRename( filenameOut, filenameIn ) )
                cerr << filenameOut << ": Error renaming " << endl;
        }

        if ( verboseEncode == 1 )
        {
            struct stat results;
            if ( stat( filenameIn, &results ) == 0 )
                // The size of the file in bytes is in results.st_size
                //fprintf(tmpFile,"%s\t%u\n", filenameIn, results.st_size);
                cerr << filenameIn << "\t" << results.st_size << endl;
            else
                //fprintf(tmpFile,"An error occurred %s\n", filenameIn);
                cerr << "An error occurred" << endl;
        }

        const bool permuteQualities = ( bwtParams_->getValue( PARAMETER_PROCESS_QUALITIES ) == PROCESS_QUALITIES_PERMUTE );
        if ( permuteQualities )
        {
            TmpFilename filenameQualIn( "", g, ".qual" );
            TmpFilename filenameQualOut( "new_", g, ".qual" );
            FILE *OutQualFileBWT = fopen( filenameQualOut, "rb" );

            if ( OutQualFileBWT != NULL ) //If it exists
            {
                fclose( OutQualFileBWT );
                if ( remove( filenameQualIn ) != 0 )
                    cerr << filenameQualIn << ": Error deleting file" << endl;
                else if ( safeRename( filenameQualOut, filenameQualIn ) )
                    cerr << filenameQualOut << ": Error renaming " << endl;
            }
        }
    }
    //cerr <<  endl;
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
        case 'N':
        case '$':
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

void BCRexternalBWT::storeBWT_parallelPile( uchar const *newSymb, uchar const *newQual, unsigned int parallelPile, SequenceNumber startIndex, SequenceNumber endIndex )
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
    LetterNumber toRead = 0;
    SequenceNumber j;
    AlphabetSymbol currentPile = parallelPile;
    const bool permuteQualities = ( bwtParams_->getValue( PARAMETER_PROCESS_QUALITIES ) == PROCESS_QUALITIES_PERMUTE );
    const bool generateCycleBwt = ( bwtParams_->getValue( PARAMETER_GENERATE_CYCLE_BWT ) != GENERATE_CYCLE_BWT_OFF );
    const bool generateCycleQualities = ( bwtParams_->getValue( PARAMETER_GENERATE_CYCLE_QUAL ) != GENERATE_CYCLE_QUAL_OFF );
    const bool isCycleBwtPBE = ( bwtParams_->getValue( PARAMETER_GENERATE_CYCLE_BWT ) == GENERATE_CYCLE_BWT_PBE );

    PredictionStatistics predictionStatistics;

    // If there are no letters to add at the last cycle, still process it
    if ( startIndex >= endIndex && debugCycle >= lengthRead && currentPile > 0 )
    {
        TmpFilename filenameIn( "", currentPile, "" );
        TmpFilename filenameOut( "new_", currentPile, "" );
        convertFileFromIntermediateToFinalFormat( filenameIn, filenameOut );
    }

    j = startIndex;
    while ( j < endIndex )
    {
        assert( currentPile == vectTriple[j].pileN );
        if ( verboseEncode == 1 )
            cerr << "index j= " << j << " current BWT segment " << ( int )currentPile << endl;

        //cerr << "Pile " << (int)currentPile << endl;
        TmpFilename filenameIn( "", currentPile, "" );
        TmpFilename filenameOut( "new_", currentPile, "" );

        assert( currentPile > 0 ); // I removed the special case for pile 0
        unique_ptr<BwtReaderBase> pReader( instantiateBwtReaderForIntermediateCycle( filenameIn, true ) );
        unique_ptr<BwtWriterBase> pWriter;
        if ( debugCycle < lengthRead )
        {
            pWriter.reset( instantiateBwtWriterForIntermediateCycle( filenameOut ) );
        }
        else
        {
            pWriter.reset( instantiateBwtWriterForLastCycle( filenameOut ) );
        }
        assert( pReader != NULL );
        assert( pWriter != NULL );


        unique_ptr<BwtReaderBase> pQualReader;
        unique_ptr<BwtWriterBase> pQualWriter;
        if ( permuteQualities )
        {
            TmpFilename filenameQualIn( "", currentPile, ".qual" );
            TmpFilename filenameQualOut( "new_", currentPile, ".qual" );
            pQualReader.reset( new BwtReaderASCII( filenameQualIn ) );
            pQualWriter.reset( new BwtWriterASCII( filenameQualOut ) );
        }

        unique_ptr<BwtWriterBase> pWriterPredictionBasedEncoding;
        unique_ptr<BwtWriterBase> pQualWriterPredictionBasedEncoding;
        if ( generateCycleBwt )
        {
            Filename pbeFilenameOut( "cycle", debugCycle, "_pile", currentPile, isCycleBwtPBE ? ".pbe" : ".ascii" );
            pWriterPredictionBasedEncoding.reset( new BwtWriterASCII( pbeFilenameOut ) );
        }

        if ( generateCycleQualities )
        {
            Filename pbeFilenameQualOut( "cycle", debugCycle, "_pile", currentPile, ".pbe.qual" );
            pQualWriterPredictionBasedEncoding.reset( new BwtWriterASCII( pbeFilenameQualOut ) );
        }

        //For each new symbol in the same pile
        SequenceNumber k = j;
        LetterNumber cont = 0;
        while ( ( k < nText ) && ( vectTriple[k].pileN == currentPile ) )
        {
            if ( verboseEncode == 1 )
                cerr << "k= " << k << " Q[k]= " << ( int )vectTriple[k].pileN << " P[k]= " << vectTriple[k].posN << " cont = " << cont << endl;
            //cerr << "k= " << k << " pileN[k]= " << pileN[k] << " posN[k]= " << posN[k] << endl;
            //So I have to read the k-BWT and I have to count the number of the symbols up to the position posN.
            //As PosN starts to the position 1 and I have to insert the new symbol in position posN[k]
            // I have to read posN[k]-1 symbols
            //cont is the number of symbols already read!
            toRead = ( vectTriple[k].posN - 1 ) - cont;
            if ( toRead )
            {
                Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "Start: to Read " << toRead << "\n";

                ( *pReader ).readAndSend( ( *pWriter ), toRead );
                if ( permuteQualities )
                {
                    ( *pQualReader ).readAndSend( ( *pQualWriter ), toRead );
                }

                cont += toRead;
            }


            //Now I have to insert the new symbol associated with the suffix of the sequence k
            //And I have to update the number of occurrences of each symbol
            //   if (toRead==0) {

            const char charToWrite = newSymb[vectTriple[k].seqN];
            if ( generateCycleBwt || generateCycleQualities )
            {
                char predictedChar = pWriter->getLastChar();
                /*
                                static int correctPredictions = 0;
                                static int incorrectPredictions = 0;
                                if ( charToWrite == predictedChar )
                                {
                                    ++correctPredictions;
                                }
                                else
                                {
                                    ++incorrectPredictions;
                                }
                                clog << "correctPrediction=" << (charToWrite==predictedChar) << " prediction=" << predictedChar << " base=" << charToWrite << endl;
                */
                bool isCorrectlyPredicted;
                const char encodedChar = getPredictionBasedEncoding( charToWrite, predictedChar, isCorrectlyPredicted );
                if ( generateCycleBwt )
                    ( *pWriterPredictionBasedEncoding )( isCycleBwtPBE ? ( char * )&encodedChar : ( char * )&charToWrite, 1 );

                // Prediction-based encoded Qualities
                if ( generateCycleQualities && newQual )
                {
                    const char qual = newQual[vectTriple[k].seqN];
                    bool removeQuality;
                    if ( qual <= 33 + 2 )
                    {
                        // Don't remove QScores 0-2
                        removeQuality = false;
                    }
                    else
                    {
                        removeQuality = isCorrectlyPredicted;
                    }

                    const char qualToWrite = removeQuality ? 255 : qual;
                    ( *pQualWriterPredictionBasedEncoding )( ( char * )&qualToWrite, 1 );
                    if ( qual ) // qual==0 for '$' signs
                        predictionStatistics.add( removeQuality, qual );
                }
            }
            ( *pWriter )( ( char * )&charToWrite, 1 );
            if ( permuteQualities && newQual )
            {
                ( *pQualWriter )( ( char * )&newQual[vectTriple[k].seqN], 1 );
            }

            tableOcc_[currentPile].count_[whichPile[( int )newSymb[vectTriple[k].seqN]]]++;     //update the number of occurrences in BWT of the pileN[k]
            //cerr << "new number write " << numchar << "\n";
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

    if ( generateCycleQualities )
    {
        Filename pbeStatsFilename( "cycle", debugCycle, "_pile", currentPile, ".pbe.qual.stats" );
        predictionStatistics.outputToFile( pbeStatsFilename );
    }
}

void BCRexternalBWT::storeEntireBWT( const string &fn )
{

    LetterNumber numchar = 0;
    LetterNumber numcharWrite = 0;

    uchar *buffer = new uchar[SIZEBUFFER];

    LetterNumber *freqOut = new LetterNumber [256];
    for ( unsigned i = 0; i < 255; ++i )
        freqOut[i] = 0;

    FILE *OutFileBWT = fopen( fn.c_str(), "wb" );
    if ( OutFileBWT == NULL )
    {
        cerr << "storeEntireBWT: Error opening " << endl;
        exit ( EXIT_FAILURE );
    }

    if ( verboseEncode == 1 )
    {
        cerr << "\nThe last BWT-segment:" << endl;
        unsigned int mmm = 0;
        while ( mmm < alphabetSize )
        {
            TmpFilename filenameIn( mmm );
            //printf("===Current BWT-partial= %d\n",mmm);
            FILE *InFileBWT = fopen( filenameIn, "rb" );
            for ( LetterNumber g = 0 ; g < SIZEBUFFER; g++ )
                buffer[g] = '\0';
            numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, InFileBWT );
            cerr << "B[" << ( int )mmm << "]:\t";
            if ( numchar == 0 )
                cerr  << "empty";
            else
                cerr  << buffer;
            while ( numchar != 0 )
            {
                for ( LetterNumber g = 0 ; g < SIZEBUFFER; g++ )
                    buffer[g] = '\0';
                numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, InFileBWT );
                if ( numchar != 0 )
                    cerr  << buffer;
            }
            cerr << endl;

            fclose( InFileBWT );
            mmm++;
        }
    }

    cerr << "Entire BWT file" << endl;
    cerr << "Concatenation of " << ( int )alphabetSize << " segments \n";

    cerr << "Compute the distribution of chars \n";

    for ( AlphabetSymbol g = 0 ; g < alphabetSize; g++ )
    {
        TmpFilename filenameIn( g );
        FILE *InFileBWT = fopen( filenameIn, "rb" );
        if ( InFileBWT == NULL )
        {
            cerr << "storeEntireBWT: " << "BWT file " << ( int )g << ": Error opening " << endl;
            exit ( EXIT_FAILURE );
        }
        //cerr << "BWT file " << (int)g << "= ";
        while ( numchar != 0 )
        {
            numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, InFileBWT );
            //cerr << "number read " << numchar << "\n";
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
        cerr << "\nThe Entire BWT:" << endl;
        OutFileBWT = fopen( fn.c_str(), "rb" );
        if ( OutFileBWT == NULL )
        {
            cerr << "storeEntireBWT: Error opening " << endl;
            exit ( EXIT_FAILURE );
        }

        for ( SequenceNumber g = 0 ; g < SIZEBUFFER; g++ )
            buffer[g] = '\0';
        numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, OutFileBWT );
        if ( numchar == 0 )
            cerr  << "empty\n";
        else
            cerr  << buffer << "\n";
        fclose( OutFileBWT );
    }
    delete [] buffer;


    cerr << "Distribution in BWT\n";
    for ( AlphabetSymbol i = 0; i < 255; ++i )
        if ( freqOut[i] > 0 )
            cerr << i << " " << freqOut[i] << "\n";
    delete [] freqOut;
}

void BCRexternalBWT::storeSA( SequenceLength iterationNum )
{

    //I have found the position where I have to insert the chars in the position t of the each text
    //Now I have to update the SA in each file.

    LetterNumber numchar = 0;
    LetterNumber numcharWrite = 0;
    ElementType *buffer = new ElementType[SIZEBUFFER];
    LetterNumber toRead = 0;

    SequenceNumber j = 0;
    while ( j < nText )
    {
        AlphabetSymbol currentPile = vectTriple[j].pileN;
        //if (verboseEncode==1)
        // cerr << "\nNew Segment; index text j= " << j << " current SA segment is " << (int)currentPile << endl;
        //cerr << "Pile " << (int)currentPile << endl;
        TmpFilename filenameIn( "sa_", currentPile );
        FILE *InFileSA = fopen( filenameIn, "rb" );
        if ( InFileSA == NULL )
        {
            cerr << "In SA file " << ( int )j << ": Error opening " << endl;
            exit ( EXIT_FAILURE );
        }

        TmpFilename filenameOut( "new_sa_", currentPile );
        FILE *OutFileSA = fopen( filenameOut, "wb" );
        if ( OutFileSA == NULL )
        {
            cerr << "Out SA file " << ( int )j << ": Error opening " << endl;
            exit ( EXIT_FAILURE );
        }
        //cerr << "In File " << filenameIn << endl;
        //cerr << "Out File " << filenameOut << endl;

        //For each new symbol in the same pile
        SequenceNumber k = j;
        LetterNumber cont = 0;
        while ( ( k < nText ) && ( vectTriple[k].pileN == currentPile ) )
        {

            //if (verboseEncode==1)
            // cerr << "k= " << k << " Q[k]= " << (int)vectTriple[k].pileN << " P[k]= " << vectTriple[k].posN << " cont = "<< cont << endl;
            //cerr << "k= " << k << " pileN[k]= " << pileN[k] << " posN[k]= " << posN[k] << endl;
            //So I have to read the k-SA and I have to count the number of the symbols up to the position posN.
            //As PosN starts to the position 1 and I have to insert the new symbol in position posN[k]
            // I have to read posN[k]-1 symbols
            //cont is the number of symbols already read!
            toRead = ( vectTriple[k].posN - 1 ) - cont;
            /*
            if (verboseEncode == 1)
                cerr << "Start: to Read " << toRead << "\n";
            */
            while ( toRead > 0 )            //((numchar!=0) && (toRead > 0)) {
            {
                if ( toRead < SIZEBUFFER ) //The last reading for this sequence
                {
                    numchar = fread( buffer, sizeof( ElementType ), toRead, InFileSA );
                    /*
                    if (verboseEncode == 1)
                        cerr << "number read " << numchar << " to Read " << toRead << "\n";
                    */
                    checkIfEqual( numchar, toRead ); // we should always read/write the same number of characters

                    numcharWrite = fwrite ( buffer, sizeof( ElementType ), numchar , OutFileSA );
                    checkIfEqual( numchar, numcharWrite ); // we should always read/write the same number of characters
                    //cerr << "toread number write " << numcharWrite << "\n";
                }
                else
                {
                    numchar = fread( buffer, sizeof( ElementType ), SIZEBUFFER, InFileSA );
                    //if (verboseEncode == 1)
                    // cerr << "number read " << numchar << "\n";
                    checkIfEqual( numchar, SIZEBUFFER ); // we should always read/write the same number of characters
                    numcharWrite = fwrite ( buffer, sizeof( ElementType ), numchar , OutFileSA );
                    checkIfEqual( numchar , numcharWrite ); // we should always read/write the same number of characters
                    //cerr << "sizebuffer number write " << numcharWrite << "\n";
                }

                cont   += numchar;  //number of read symbols
                toRead -= numchar;
                if ( ( numchar == 0 ) && ( toRead > 0 ) ) //it means that we have read 0 character, but there are still toRead characters to read
                {
                    cerr << "storeSA: sequence number" << ( int )k << " read 0 character, but there are still " << toRead << " characters to read  " << endl;
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
                //cerr << "new number write " << numchar << "\n";
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
            //cerr << "After insert: " << numchar << "\n";
            numcharWrite = fwrite ( buffer, sizeof( ElementType ), numchar , OutFileSA );
            checkIfEqual( numchar, numcharWrite ); // we should always read/write the same number of characters
        }

        fclose( InFileSA );
        fclose( OutFileSA );
        j = k;
    }

    //Renaming new to old
    for ( AlphabetSymbol g = 0 ; g < alphabetSize; g++ )
    {
        TmpFilename filenameIn( "sa_", g );
        TmpFilename filenameOut( "new_sa_", g );
        //cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << endl;
        FILE *OutFileSA = fopen( filenameOut, "rb" );

        if ( OutFileSA != NULL ) //If it exists
        {
            fclose( OutFileSA );
            if ( remove( filenameIn ) != 0 )
                cerr << filenameIn << ": Error deleting file" << endl;
            else if ( safeRename( filenameOut, filenameIn ) )
                cerr << filenameOut << ": Error renaming " << endl;
        }
        /*
        if (verboseEncode == 1) {
            struct stat results;
            if (stat(filenameIn, &results) == 0)
                // The size of the file in bytes is in results.st_size
                //fprintf(tmpFile,"%s\t%u\n", filenameIn, results.st_size);
                cerr << filenameIn <<"\t" << results.st_size << endl;
            else
                //fprintf(tmpFile,"An error occurred %s\n", filenameIn);
                cerr << "An error occurred" << endl;
        }
        */
    }
    //cerr <<  endl;
    //fprintf(tmpFile,"\n");
    //fclose(tmpFile);

    delete [] buffer;
}

void BCRexternalBWT::storeEntirePairSA( const char *fn )
{

    cerr << "\nEntire Pairs SA file (position, number of sequence)" << endl;

    LetterNumber numcharWrite, numcharRead;
    ElementType *buffer = new ElementType[SIZEBUFFER];

    TmpFilename fnSA( fn, ".pairSA" );

    FILE *OutFileSA = fopen( fnSA, "wb" );
    if ( OutFileSA == NULL )
    {
        cerr << "Entire Pairs SA file: Error opening " << fnSA << endl;
        exit ( EXIT_FAILURE );
    }
    /* //it will be useful for varying length reads
        vector <SequenceLength> vectLen;
        vectLen.resize(nText);

        char *fileLen="outFileLen";
        FILE *InFileLen;                  // file of the lengths;
        InFileLen = fopen(fileLen, "rb");
        if (InFileLen==NULL) {
                cerr << "storeEntireSAfromPairSA: could not open file \"" << fileLen << "\"!"<< endl;
                exit (EXIT_FAILURE);
        }

        numcharRead = fread (&vectLen[0], sizeof(SequenceLength), vectLen.size() , InFileLen);
        checkIfEqual(numcharRead , nText); // we should always read the same number of characters

        fclose(InFileLen);
        */
    for ( AlphabetSymbol g = 0 ; g < alphabetSize; g++ )
    {
        TmpFilename filenameIn( "sa_", g );
        FILE *InFileSA = fopen( filenameIn, "rb" );
        if ( InFileSA == NULL )
        {
            cerr << "SA file " << ( int )g << ": Error opening " << endl;
            exit ( EXIT_FAILURE );
        }
        //cerr << "SA file " << (int)g << "= ";

        numcharRead = fread( buffer, sizeof( ElementType ), SIZEBUFFER, InFileSA );

        /* //it will be useful for varying length reads
        //Correction of the length of the sequences.
        for (SequenceNumber num = 0; num < numcharRead; num++) {
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
            for (SequenceNumber num = 0; num < numcharRead; num++) {
                buffer[num].sa = buffer[num].sa - (lengthRead - vectLen[buffer[num].numSeq]);
            }
            */
            numcharWrite = fwrite ( buffer, sizeof( ElementType ), numcharRead , OutFileSA );
            checkIfEqual ( numcharRead , numcharWrite );
        }

        fclose( InFileSA );
        if ( remove( filenameIn ) != 0 )
            cerr << filenameIn << ": Error deleting file" << endl;
    }

    fclose( OutFileSA );

    if ( verboseEncode == 1 )
    {
        OutFileSA = fopen( fnSA, "rb" );
        if ( OutFileSA == NULL )
        {
            cerr << "Entire SA file: Error opening " << endl;
            exit ( EXIT_FAILURE );
        }

        numcharRead = fread( buffer, sizeof( ElementType ), SIZEBUFFER, OutFileSA );
        if ( numcharRead == 0 )
            cerr  << "empty\n";
        else
            for ( LetterNumber g = 0 ; g < numcharRead; g++ )
            {
                cerr  << "(" << ( int )buffer[g].sa << "," << buffer[g].numSeq << ") ";
            }
        while ( numcharRead != 0 )
        {
            numcharRead = fread( buffer, sizeof( ElementType ), SIZEBUFFER, OutFileSA );
            for ( LetterNumber g = 0 ; g < numcharRead; g++ )
            {
                cerr  << "(" << buffer[g].sa << "," << buffer[g].numSeq << ") ";
            }
        }
        cerr << endl;

        fclose( OutFileSA );
    }

    delete [] buffer;
}


void BCRexternalBWT::storeEntireSAfromPairSA( const char *fn )
{
    cerr << "\nSA file from pair SA file" << endl;

    LetterNumber numchar, numcharWrite;
    /* //it will be useful for varying length reads

        vector <SequenceLength> vectSumCumLen;
        vectSumCumLen.resize(nText+1);
        char *fileLen="outFileLen";
        FILE *InFileLen;                  // file of the lengths;
        InFileLen = fopen(fileLen, "rb");
        if (InFileLen==NULL) {
                cerr << "storeEntireSAfromPairSA: could not open file \"" << fileLen << "\"!"<< endl;
                exit (EXIT_FAILURE);
        }
        SequenceLength lenSeq=0;
        numchar = fread (&lenSeq, sizeof(SequenceLength), 1 , InFileLen);
        checkIfEqual( numchar , 1); // we should always read the same number of characters
        vectSumCumLen[0] = 0;
        vectSumCumLen[1] = lenSeq + 1;   //Plus $
        for (SequenceNumber num = 2; num < nText+1; num++) {
            numchar = fread (&lenSeq, sizeof(SequenceLength), 1 , InFileLen);
            checkIfEqual(numchar , 1); // we should always read the same number of characters
            vectSumCumLen[num] = vectSumCumLen[num-1] + lenSeq + 1;  //Plus $
        }
        fclose(InFileLen);
    */
    TmpFilename fnSA    ( fn, ".sa" );
    TmpFilename fnPairSA( fn, ".pairSA" );

    FILE *InFilePairSA = fopen( fnPairSA, "rb" );
    if ( InFilePairSA == NULL )
    {
        cerr << "Entire Pairs SA file: Error opening " << fnPairSA << endl;
        exit ( EXIT_FAILURE );
    }

    FILE *OutFileSA = fopen( fnSA, "wb" );
    if ( OutFileSA == NULL )
    {
        cerr << "Entire SA file: Error opening " << fnSA << endl;
        exit ( EXIT_FAILURE );
    }

    ElementType *buffer = new ElementType[SIZEBUFFER];
    LetterNumber *bufferNChar = new LetterNumber[SIZEBUFFER];

    while ( !feof( InFilePairSA ) )
    {
        numchar = fread( buffer, sizeof( ElementType ), SIZEBUFFER, InFilePairSA );
        //cerr << "number read " << numchar << "\n";
        if ( numchar > 0 )
        {
            for ( LetterNumber i = 0; i < numchar; i++ )
            {
                bufferNChar[i] = ( LetterNumber )( buffer[i].numSeq * ( lengthRead + 1 ) + buffer[i].sa );
                //cerr << buffer[i].numSeq << " " << (int)lengthRead << " " << (int)buffer[i].sa << " --> " << (int)bufferNChar[i] << "\n";
                //bufferNChar[i] = (LetterNumber)(vectSumCumLen[buffer[i].numSeq] + buffer[i].sa);       //it will be useful for varying length reads
                //cerr << "vectSumCumLen["<< buffer[i].numSeq<< "]= " << (int)vectSumCumLen[buffer[i].numSeq] << " + " << (int)buffer[i].sa << " --> " << (int)bufferNChar[i] << "\n";
            }
            numcharWrite = fwrite ( bufferNChar, sizeof( LetterNumber ), numchar, OutFileSA );
            checkIfEqual( numchar , numcharWrite );
            //cerr << "number write " << numcharWrite << "\n";
        }
    }
    fclose( InFilePairSA );
    fclose( OutFileSA );

    if ( verboseEncode == 1 )
    {
        cerr << "\nThe Entire SA. The file is " << fnSA << endl;
        OutFileSA = fopen( fnSA, "rb" );
        if ( OutFileSA == NULL )
        {
            cerr << "Entire SA file: Error opening " << endl;
            exit ( EXIT_FAILURE );
        }

        numchar = fread( bufferNChar, sizeof( LetterNumber ), SIZEBUFFER, OutFileSA );
        if ( numchar == 0 )
            cerr  << "empty\n";
        else
            for ( LetterNumber g = 0 ; g < numchar; g++ )
            {
                cerr  << bufferNChar[g] << " ";
            }
        while ( numchar != 0 )
        {
            numchar = fread( buffer, sizeof( ElementType ), SIZEBUFFER, OutFileSA );
            for ( LetterNumber g = 0 ; g < numchar; g++ )
            {
                cerr  << bufferNChar[g] << " ";
            }
        }
        cerr << endl;

        fclose( OutFileSA );
    }

    delete [] buffer;
    delete [] bufferNChar;
}

void BCRexternalBWT::storeBWTandLCP( uchar const *newSymb )
{
    SequenceLength maxValueLen = lengthRead + 2;
    vector <SequenceLength> minLCPcur;
    vector <bool> minLCPcurFound;
    vector <SequenceLength> minLCPsuc;
    vector <SequenceNumber> minLCPsucText;
    vector <bool> minLCPsucToFind;
    minLCPcur.resize( alphabetSize );     //for each symbol of the alphabet
    minLCPcurFound.resize( alphabetSize );
    minLCPsuc.resize( alphabetSize );
    minLCPsucText.resize( alphabetSize );
    minLCPsucToFind.resize( alphabetSize );

    //I have found the position where I have to insert the chars in the position t of the each text
    //Now I have to update the BWT in each file.
    LetterNumber numchar = 0;
    LetterNumber numcharWrite = 0;
    uchar *buffer = new uchar[SIZEBUFFER];
    SequenceLength *bufferLCP = new SequenceLength[SIZEBUFFER];
    LetterNumber toRead = 0;

    SequenceNumber j = 0;
    while ( j < nText )
    {
        AlphabetSymbol currentPile = vectTriple[j].pileN;
        for ( AlphabetSymbol g = 0 ; g < alphabetSize; g++ )
        {
            minLCPcur[g] = maxValueLen;
            minLCPcurFound[g] = 0;
            minLCPsuc[g] = maxValueLen;
            minLCPsucText[g] = 0;   //denotes the number of the text associated with the symbol g
            minLCPsucToFind[g] = 0;
        }
        assert( currentPile < alphabetSize );
        TmpFilename filenameIn( currentPile );
        FILE *InFileBWT = fopen( filenameIn, "rb" );
        if ( InFileBWT == NULL )
        {
            cerr << "(storeBWTandLCP) In BWT file " << filenameIn << ": Error opening " << endl;
            exit ( EXIT_FAILURE );
        }
        TmpFilename filenameOut( "new_", currentPile );
        FILE *OutFileBWT = fopen( filenameOut, "wb" );
        if ( OutFileBWT == NULL )
        {
            cerr << "(storeBWTandLCP) Out BWT file " << filenameIn << ": Error opening " << endl;
            exit ( EXIT_FAILURE );
        }
        TmpFilename filenameInLCP( "", currentPile, ".lcp" );
        FILE *InFileLCP = fopen( filenameInLCP, "rb" );
        if ( InFileLCP == NULL )
        {
            cerr << "In LCP file " << filenameInLCP << ": Error opening "  << endl;
            exit ( EXIT_FAILURE );
        }
        TmpFilename filenameOutLCP( "new_", currentPile, ".lcp" );
        FILE *OutFileLCP = fopen( filenameOutLCP, "wb" );
        if ( OutFileLCP == NULL )
        {
            cerr << "Out LCP file " << filenameInLCP << ": Error opening " << filenameOutLCP << endl;
            exit ( EXIT_FAILURE );
        }

        //For each new symbol in the same pile
        SequenceNumber k = j;
        LetterNumber cont = 0;
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
                    numchar = fread( bufferLCP, sizeof( SequenceLength ), toRead, InFileLCP );
                    checkIfEqual( numchar , toRead );
                    numcharWrite = fwrite ( bufferLCP, sizeof( SequenceLength ), numchar , OutFileLCP );
                    checkIfEqual( numchar , numcharWrite );
                }
                else
                {
                    numchar = fread( buffer, sizeof( uchar ), SIZEBUFFER, InFileBWT );
                    checkIfEqual( numchar , SIZEBUFFER );
                    numcharWrite = fwrite ( buffer, sizeof( uchar ), numchar , OutFileBWT );
                    checkIfEqual( numchar , numcharWrite );
                    numchar = fread( bufferLCP, sizeof( SequenceLength ), SIZEBUFFER, InFileLCP );
                    checkIfEqual( numchar , SIZEBUFFER );
                    numcharWrite = fwrite ( bufferLCP, sizeof( SequenceLength ), numchar , OutFileLCP );
                    checkIfEqual( numchar , numcharWrite );
                }
                //I must to compute the minimum LCP. It needs to compute the lcpValue for the next iteration
                //cerr << "For each letter in the buffer before of the position where I have to insert the new symbol\n";
                for ( LetterNumber bb = 0 ; bb < numcharWrite; bb++ )
                {
                    //Update the min1 for each letter of the alphabet, for which I have already met the symbol
                    for ( AlphabetSymbol gg = 0 ; gg < alphabetSize; gg++ )
                    {
                        if ( minLCPcurFound[gg] == 1 ) //I already met the symbol gg. So, I must compute the minimum
                            if ( bufferLCP[bb] < minLCPcur[gg] ) //comparison with the last inserted lcp
                                minLCPcur[gg] = bufferLCP[bb];
                    }

                    minLCPcur[whichPile[( int )buffer[bb]]] = maxValueLen; //For each occurrence of buffer[bb], I have to set the min1 (the interval starts from the next symbol)
                    minLCPcurFound[whichPile[( int )buffer[bb]]] = 1; //So I open the LCP interval for buffer[bb] (for the first occurrence of buffer[bb])

                    //First, it needs to check if the symbol buffer[bb] closes a previous interval or it is in the middle or no.
                    //In any case, for each symbol, we have to update the minimum
                    for ( AlphabetSymbol gg = 0 ; gg < alphabetSize; gg++ )
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
                    cerr << "Error storeBWT: sequence number" << ( int )k << " read 0 character, but there are still " << toRead << " characters to read  " << endl;
                    exit ( EXIT_FAILURE );
                }

            }
            //Now I have to insert the new symbol associated with the suffix of the sequence k
            //And I have to update the number of occurrences of each symbol
            //And I have to insert the valueLCP store in lcpCurN + 1 in the previous iteration
            if ( toRead == 0 )
            {
                //cerr << "\nNow I can insert the new symbol and lcp, indeed toRead= " << toRead << endl;
                numchar = fwrite ( &newSymb[vectTriple[k].seqN], sizeof( uchar ), 1, OutFileBWT );
                checkIfEqual( numchar , 1 ); // we should always read/write the same number of characters
                tableOcc_[currentPile].count_[whichPile[( int )newSymb[vectTriple[k].seqN]]]++;
                //tableOcc[currentPile][whichPile[(int)newSymb[vectTriple[k].seqN]]]++;       //update the number of occurrences in BWT of the pileN[k]
                SequenceLength lcpValueNow;
                if ( vectTriple[k].posN == 1 )   //it is the first symbol of the segment. So, the lcp value is 0
                {
                    lcpValueNow = vectTriple[k].getLcpCurN();
                }
                else
                    lcpValueNow = vectTriple[k].getLcpCurN() + 1;
                numchar = fwrite ( &lcpValueNow, sizeof( SequenceLength ), 1, OutFileLCP ); //Insert the lcp for the new symbol
                checkIfEqual( numchar , 1 );
                //cerr << "I insert the symbol= " << newSymb[vectTriple[k].seqN] <<  " and lcp " << lcpValueNow << endl;
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
                for ( AlphabetSymbol gg = 0 ; gg < alphabetSize; gg++ )
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
                //cerr << "minSuc was not the maxValue. The new value for minLCPsuc[" << newSymb[vectTriple[k].seqN] << "] is " << minLCPsuc[alpha[(int)newSymb[vectTriple[k].seqN]]] << "\n";

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
                        cerr << "???? Warning!--Should be the same? triple[" << k << "].lcpSucN(=" << vectTriple[k].getLcpSucN() << ") + 1= " << vectTriple[k].getLcpSucN() + 1 << " == triple[" << k + 1 << "].lcpCurN+1= " << vectTriple[k + 1].getLcpCurN() + 1 << " ";
                        cerr << ", Seq k N. " << vectTriple[k].seqN << " and Seq k+1 N. " << vectTriple[k + 1].seqN << "\n";
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
                        SequenceLength sucLCP = 0;
                        numchar = fread( &sucLCP, sizeof( SequenceLength ), 1, InFileLCP ); //I have to change it
                        checkIfEqual( numchar , 1 ); // we should always read/write the same number of characters

                        //I have to update the lcp of this symbol and I have to copy it into the new bwt segment
                        SequenceLength lcpValueNow = vectTriple[k].getLcpSucN() + 1;
                        numcharWrite = fwrite ( &lcpValueNow , sizeof( SequenceLength ), numchar , OutFileLCP ); //Updated the lcpSuc
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
                        else    			//cerr << "The succSymb is not equal to the new symbol\n";
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
                        for ( AlphabetSymbol gg = 0 ; gg < alphabetSize; gg++ )
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
                        for ( AlphabetSymbol gg = 0 ; gg < alphabetSize; gg++ )
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
            numchar = fread( bufferLCP, sizeof( SequenceLength ), SIZEBUFFER, InFileLCP );
            numcharWrite = fwrite ( bufferLCP, sizeof( SequenceLength ), numchar , OutFileLCP );
            assert( numchar == numcharWrite ); // we should always read/write the same number of characters
            //Compute lcpSucN for the other texts
            //For each symbol in the buffer, we check it it close any interval, while each entry in minLcpSuc is maxValue

            //TBD: TO OPTIMIZE. IT CAN END BEFORE. IT DOES NOT NEED TO READ THE ENTIRE BUFFER

            for ( LetterNumber bb = 0 ; bb < numchar; bb++ )
            {
                //First, I check if the symbol bb closes a previous interval or it is in the middle or no.
                //In any case, for each symbol, we have to update the minimum
                for ( AlphabetSymbol gg = 0 ; gg < alphabetSize; gg++ )
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
        for ( AlphabetSymbol gg = 0 ; gg < alphabetSize; gg++ )
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
        //cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << endl;
        OutFileBWT = fopen( filenameOut, "r" );
        if ( OutFileBWT != NULL ) //If it exists
        {
            fclose( OutFileBWT );
            if ( remove( filenameIn ) != 0 )
                cerr << filenameIn << ": Error deleting file" << endl;
            else if ( safeRename( filenameOut, filenameIn ) )
                cerr << filenameOut << ": Error renaming " << endl;
        }
        //cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << endl;
        OutFileLCP = fopen( filenameOutLCP, "r" );
        if ( OutFileLCP != NULL ) //If it exists
        {
            fclose( OutFileLCP );
            if ( remove( filenameInLCP ) != 0 )
                cerr << filenameInLCP << ": Error deleting file" << endl;
            else if ( safeRename( filenameOutLCP, filenameInLCP ) )
                cerr << filenameOutLCP << ": Error renaming " << endl;
        }

        j = k;
    }

    if ( verboseEncode == 1 )
    {
        cerr << "After the computation of LCP for the next iteration" << endl;
        cerr << "Q  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            cerr << ( int )vectTriple[g].pileN << " ";
        }
        cerr << endl;
        cerr << "P  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            cerr << vectTriple[g].posN  << " ";
        }
        cerr << endl;
        cerr << "N  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            cerr << vectTriple[g].seqN  << " ";
        }
        cerr << endl;
        cerr << "C  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            cerr << vectTriple[g].getLcpCurN()  << " ";
        }
        cerr << endl;
        cerr << "S  ";
        for ( SequenceNumber g = 0 ; g < nText; g++ )
        {
            cerr << vectTriple[g].getLcpSucN()  << " ";
        }
        cerr << endl;
    }

    //Renaming new to old
    for ( AlphabetSymbol g = 0 ; g < alphabetSize; g++ )
    {
        TmpFilename filenameOut( "new_", g );
        TmpFilename filenameIn( g );
        //cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << endl;
        FILE *OutFileBWT = fopen( filenameOut, "rb" );
        if ( OutFileBWT != NULL ) //If it exists
        {
            fclose( OutFileBWT );
            if ( remove( filenameIn ) != 0 )
                cerr << filenameIn << ": Error deleting file" << endl;
            else if ( safeRename( filenameOut, filenameIn ) )
                cerr << filenameOut << ": Error renaming " << endl;
        }

        TmpFilename filenameOutLCP( "new_", g, ".lcp" );
        TmpFilename filenameInLCP( "", g, ".lcp" );
        //cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << endl;
        FILE *OutFileLCP = fopen( filenameOutLCP, "rb" );
        if ( OutFileLCP != NULL ) //If it exists
        {
            fclose( OutFileLCP );
            if ( remove( filenameInLCP ) != 0 )
                cerr << filenameInLCP << ": Error deleting file" << endl;
            else if ( safeRename( filenameOutLCP, filenameInLCP ) )
                cerr << filenameOutLCP << ": Error renaming " << endl;
        }
    }

    delete [] buffer;
    delete [] bufferLCP;
}

void BCRexternalBWT::storeEntireLCP( const string &fn )
{
    assert( false && "TODO" );
}

void BCRexternalBWT::pauseBetweenCyclesIfNeeded()
{
    if ( bwtParams_->getValue( PARAMETER_PAUSE_BETWEEN_CYCLES ) == true )
        pauseBetweenCycles();
}

void BCRexternalBWT::writeEndPosFile( const uint8_t subSequenceNum, const bool lastFile )
{
    Logger_if( LOG_SHOW_IF_VERBOSE )
    {
        Logger::out() << "Storing 'outFileEndPos', time now: " << timer.timeNow();
        Logger::out() << "Storing 'outFileEndPos', usage: " << timer << endl;

        Logger::out() << "Stores the 'end positions' of the $!" << endl;
    }
    LetterNumber numchar;
    FILE *inFileEndPos = NULL;
    FILE *OutFileEndPos = NULL;
    static bool firstTime = true;
    char inBuf[ sizeof( SequenceNumber ) + sizeof( uint8_t ) ];
    if ( !firstTime )
    {
        Filename previousFileEndPos( bwtParams_->getStringValue( "output filename" ).c_str(), "-end-pos-intermediate" );
        inFileEndPos = fopen( previousFileEndPos, "rb" );
        SequenceNumber skipNumText;
        uint8_t skipSubSequenceCount, skipHasRevComp;
        fread ( &skipNumText, sizeof( SequenceNumber ), 1 , inFileEndPos );
        fread ( &skipSubSequenceCount, sizeof( uint8_t ), 1 , inFileEndPos );
        fread ( &skipHasRevComp, sizeof( uint8_t ), 1 , inFileEndPos );
    }
    Filename fileEndPos( bwtParams_->getStringValue( "output filename" ).c_str(), lastFile ? "-end-pos" : "-end-pos-intermediate" );
    OutFileEndPos = fopen( fileEndPos, "wb" );
    if ( OutFileEndPos == NULL )
    {
        cerr << "Error opening \"" << fileEndPos << "\" file" << endl;
        exit ( EXIT_FAILURE );
    }

    // Writing -end-pos file header
    uint8_t hasRevComp = ( bwtParams_->getValue( PARAMETER_ADD_REV_COMP ) == 0 ) ? 0 : 1;
    SequenceNumber seqCount = nText;
    if ( hasRevComp )
    {
        assert( seqCount % 2 == 0 );
        seqCount /= 2;
    }
    uint8_t subSequenceCount = 1;

    if ( bwtParams_->getValue( PARAMETER_PAIRED_READS_INPUT ) == PAIRED_READS_INPUT_ALL1ALL2 )
    {
        assert( firstTime ); // We can't have both paired-end ways (--sub-sequence-length and --paired-end-input) at the same time

        assert( seqCount % 2 == 0 );
        seqCount /= 2;
        subSequenceCount = 2;
    }
    else
    {
        subSequenceCount = firstTime ? 1 : 2;
    }
    numchar = fwrite ( &seqCount, sizeof( SequenceNumber ), 1 , OutFileEndPos );
    assert( numchar == 1 );
    numchar = fwrite ( &subSequenceCount, sizeof( uint8_t ), 1 , OutFileEndPos );
    assert( numchar == 1 );
    numchar = fwrite ( &hasRevComp, sizeof( uint8_t ), 1 , OutFileEndPos );
    assert( numchar == 1 );


    LetterNumber cont = 0;
    LetterCount counters;
    int currentPile = 0;
    unique_ptr<BwtReaderBase> pReader;
    for ( SequenceNumber i = 0; i < nText; i++ )
    {
        while ( currentPile != vectTriple[i].pileN )
        {
            ++currentPile;
            // finish read&counting current bwt file
            counters.clear();
            if ( pReader )
            {
                pReader->readAndCount( counters );
            }

            // how many '$' signs found? transfer as many end-pos-intermediate entries
            Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "Transfer " << counters.count_[0] << " entries" << endl;
            for ( unsigned int j = 0; j < counters.count_[0]; ++j )
            {
                fread( inBuf, sizeof( SequenceNumber ) + sizeof( uint8_t ), 1, inFileEndPos );
                fwrite( inBuf, sizeof( SequenceNumber ) + sizeof( uint8_t ), 1, OutFileEndPos );
            }

            // open next bwt file
            TmpFilename filenameIn( "", currentPile, "" );
            pReader.reset( instantiateBwtReaderForLastCycle( filenameIn ) );
            cont = 0;
        }

        // read&count current bwt file until position vectTriple[i].posN
        LetterNumber toRead = vectTriple[i].posN - cont - 1;
        counters.clear();
        LetterNumber numberRead = ( *pReader ).readAndCount( counters, toRead );
        if ( toRead != numberRead )
        {
            cerr << "ERROR: toRead != numberRead" << endl;
            cerr << "  toRead = " << toRead << endl;
            cerr << "  numberRead=" << numberRead << endl;
            cerr << "  cont=" << cont << endl;
            cerr << "  i=" << i << endl;
            cerr << "  vectTriple[i].posN=" << vectTriple[i].posN << endl;
            assert( false );
        }
        cont += toRead;

        // how many '$' signs found? transfer as many end-pos-intermediate entries
        Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "Transfer(2) " << counters.count_[0] << " entries" << endl;
        for ( unsigned int j = 0; j < counters.count_[0]; ++j )
        {
            fread( inBuf, sizeof( SequenceNumber ) + sizeof( uint8_t ), 1, inFileEndPos );
            fwrite( inBuf, sizeof( SequenceNumber ) + sizeof( uint8_t ), 1, OutFileEndPos );
        }

        // read&count 1 char: check that it's a '$' sign
        counters.clear();
        numberRead = ( *pReader ).readAndCount( counters, 1 );
        assert( counters.count_[0] == 1 );
        ++cont;

        // insert entry
        Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << "Triple: " << vectTriple[i].seqN << " " << vectTriple[i].posN << " " << ( int )vectTriple[i].pileN << endl;
        uint8_t extendedSubSeqNum = vectTriple[i].seqN / seqCount;
        SequenceNumber seqNumWithoutRevCompOrPair = vectTriple[i].seqN % seqCount;

        numchar = fwrite ( &seqNumWithoutRevCompOrPair, sizeof( SequenceNumber ), 1 , OutFileEndPos );
        assert( numchar == 1 );
        numchar = fwrite ( &extendedSubSeqNum, sizeof( uint8_t ), 1 , OutFileEndPos );
        assert( numchar == 1 );
    }

    fclose( OutFileEndPos );
    Logger::out() << "'end positions' stored!" << endl;

    if ( !firstTime )
    {
        Filename previousFileEndPos( bwtParams_->getStringValue( "output filename" ).c_str(), "-end-pos-intermediate" );
        if ( remove( previousFileEndPos ) != 0 )
            cerr << "Error deleting file " << previousFileEndPos << endl;
    }
    firstTime = false;
}
