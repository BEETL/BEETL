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

#include "TransposeFasta.hh"

#include "Filename.hh"
#include "SeqReader.hh"
#include "Tools.hh"
#include "libzoo/util/Logger.hh"
#include "libzoo/util/TemporaryFilesManager.hh"

#include <cassert>
#include <cstdlib>

using namespace std;


TransposeFasta::TransposeFasta()
    : pReader_( NULL )
    , cycleNum_( 0 )
    , processQualities_( false )
{
    for ( int i( 0 ); i < 256; i++ ) freq[i] = 0;
}

void TransposeFasta::init( SeqReaderFile *pReader, const bool processQualities )
{
    pReader_ = pReader;
    cycleNum_ = pReader->length();
    outputFiles_.resize( pReader->length() );
    buf_.resize( pReader->length(), vector<uchar>( BUFFERSIZE ) );
    processQualities_ = processQualities;

    cerr << "Constructing TransposeFasta, found read length of "
         << cycleNum_ << endl;

    if ( processQualities_ && pReader_->thisQual()[0] == '\0' )
    {
        // If the first entry of the file (which can be fastq or any other format (raw/fasta/etc)) doesn't contain any quality info
        // , deactivate qualities processing
        processQualities_ = false;
    }
}

TransposeFasta::~TransposeFasta()
{

}


bool TransposeFasta::convert( /*const string &input,*/ const string &output, bool generatedFilesAreTemporary )
{
    vector<vector<uchar> > bufQual;
    vector<FILE *> outputFilesQual;
    if ( processQualities_ )
    {
        bufQual.resize( pReader_->length(), vector<uchar>( BUFFERSIZE ) );
        outputFilesQual.resize( pReader_->length() );
    }

    //TO DO
    lengthRead = cycleNum_;
    //The distribution of characters is useful
    //for alpha[256] -->Corresponding between the alphabet, the piles and tableOcc
    //and to know sizeAlpha
    //We supposed that the symbols in the input file are the following
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

    // create output files
    for ( SequenceLength i = 0; i < cycleNum_; i++ )
    {
        Filename fn( output, i, "" );
        outputFiles_[i] = fopen( fn, "w" );
        if ( outputFiles_[i] == NULL )
        {
            cerr << "Error: couldn't open output file " << fn << endl;
            if ( i > 0 )
            {
                cerr << "  You may have reached the maximum number of opened files (see `ulimit -n`) or the maximum number of files allowed in one directory, as we create one file per cycle (and a second one if qualities are present)" << endl;
                exit ( -1 );
            }
        }
        if ( generatedFilesAreTemporary )
            TemporaryFilesManager::get().addFilename( fn );
        if ( processQualities_ )
        {
            Filename fnQual( output + "qual.", i, "" );
            outputFilesQual[i] = fopen( fnQual, "w" );
            if ( outputFilesQual[i] == NULL )
            {
                cerr << "Error: couldn't open output file " << fnQual << endl;
                if ( i > 0 )
                {
                    cerr << "  You may have reached the maximum number of opened files (see `ulimit -n`) or the maximum number of files allowed in one directory, as we create one file per cycle (and a second one if qualities are present)" << endl;
                    exit ( -1 );
                }
            }
            if ( generatedFilesAreTemporary )
                TemporaryFilesManager::get().addFilename( fnQual );
        }
    }


    // looping through the input file, add the characters to the buffer, print buffer when it's full
    //    unsigned int num_read = 0;
    unsigned int num_write = 0;
    unsigned int charsBuffered = 0;

    //******************************buf[cycleNum_+1];  ********* is cycleNum_ right?
    lengthTexts = 0;
    nSeq = 0;
    //    num_read = fread(buf,sizeof(uchar),cycleNum_,ifile);

    //    fgets ( buf,1024, ifile ); %%%%%
    //    while( !feof(ifile) ) %%%%%
    while ( pReader_->allRead() == false )
    {
        //cerr << "current line : " << buf << endl;

        if ( charsBuffered == BUFFERSIZE )
        {
            // write buffers to the files, clear buffers
            #pragma omp parallel for num_threads(4)
            for ( SequenceLength i = 0; i < cycleNum_; i++ )
            {
                //cerr << "writing to " << i << " : " << buf_[i] << endl;
                size_t num_write_bases = fwrite ( buf_[i].data(), sizeof( char ), charsBuffered, outputFiles_[i] );
                checkIfEqual( num_write_bases, charsBuffered ); // we should always read/write the same number of characters
                if ( processQualities_ )
                {
                    size_t num_write_qual = fwrite ( bufQual[i].data(), sizeof( char ), charsBuffered, outputFilesQual[i] );
                    checkIfEqual( num_write_bases, num_write_qual );
                }
            }
            lengthTexts += ( num_write * cycleNum_ );


            charsBuffered = 0;
        }

        for ( SequenceLength i = 0; i < cycleNum_; i++ )
        {
            buf_[i][charsBuffered] = pReader_->thisSeq()[i];

            if ( processQualities_ )
            {
                bufQual[i][charsBuffered] = pReader_->thisQual()[i];
            }
        }
        // increase the counter of chars buffered
        charsBuffered++;
        nSeq++;


#ifdef XXX
        // process the input
        if ( buf[0] != '>' )
        {
            // add the characters
            for ( SequenceLength i = 0; i < cycleNum_; i++ )
            {
                buf_[i][charsBuffered] = buf[i];
            }

            // increase the counter of chars buffered
            charsBuffered++;
            nSeq++;
        }
#endif                //else


        //num_read = fread(buf,sizeof(uchar),cycleNum_,ifile);
        //        fgets ( buf, 1024, ifile );
        pReader_->readNext();
    }

    // write the rest
    for ( SequenceLength i = 0; i < cycleNum_; i++ )
    {
        num_write = fwrite ( buf_[i].data(), sizeof( uchar ), charsBuffered, outputFiles_[i] );
        lengthTexts += num_write;
        if ( processQualities_ )
        {
            size_t num_write_qual = fwrite ( bufQual[i].data(), sizeof( uchar ), charsBuffered, outputFilesQual[i] );
            checkIfEqual( num_write, num_write_qual );
        }
    }
    checkIfEqual( num_write, charsBuffered );



    // closing all the output file streams
    for ( SequenceLength i = 0; i < cycleNum_; i++ )
    {
        fclose( outputFiles_[i] );
        if ( processQualities_ )
        {
            fclose( outputFilesQual[i] );
        }
    }

    std::cout << "Number of sequences reading/writing: " << nSeq << "\n";
    std::cout << "Number of characters reading/writing: " << lengthTexts << "\n";

    //    delete pReader;
    return true;
}

bool TransposeFasta::inputCycFile( const string &cycPrefix )
{
    //TO DO
    //The distribution of characters is useful
    //for alpha[256] -->Corresponding between the alphabet, the piles and tableOcc
    //and to know sizeAlpha

    //1) Alphabet
    //We supposed that the symbols in the input file are the following
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

    //2) Number of sequences
    string cyc1Filename = cycPrefix + "1";
    FILE *f = fopen( cyc1Filename.c_str(), "rb" );
    if ( !f )
    {
        cerr << "ERROR: Cycle file " << cyc1Filename << " not found!" << endl;
        exit( -1 );
    }
    fseek( f, 0, SEEK_END );
    nSeq = ftell( f );
    if (nSeq != ftell( f))
    {
        Logger::error() << "Error: Too many sequences. This version of BEETL was compiled for a maximum of " << maxSequenceNumber << " sequences, but this input has " << ftell(f) << " sequences. You can increase this limit by changing the type definition of 'SequenceNumber' in Types.hh and recompiling BEETL." << endl;
        exit( -1 );
    }
    fclose( f );

    //3) Length of the longest sequence
    for ( lengthRead = 1; ; ++lengthRead )
    {
        Filename cycFilename( cycPrefix, lengthRead, "" );
        FILE *f = fopen( cycFilename, "rb" );
        if ( f )
            fclose( f );
        else
            break;
    }

    //4) qualities detection
    string qual1Filename = cycPrefix + "qual.1";
    f = fopen( qual1Filename.c_str(), "rb" );
    if ( f )
    {
        processQualities_ = true;
        fclose( f );
    }
    else
        processQualities_ = false;

    //5) Total Length
    lengthTexts = lengthRead * nSeq;

    // Report
    Logger_if( LOG_SHOW_IF_VERBOSE )
    {
        Logger::out() << "****processing qualities: " << processQualities_ << "\n";
        Logger::out() << "****number of sequences: " << nSeq << "\n";
        Logger::out() << "****max length of each sequence: " << lengthRead << "\n";
        Logger::out() << "****lengthTot: " << lengthTexts << "\n";
    }

    return 1;
}


bool TransposeFasta::convertFromCycFileToFastaOrFastq( const string &fileInputPrefix, const string &fileOutput, bool generatedFilesAreTemporary, SequenceExtractor *sequenceExtractor )
{
    bool outputIsFastq = hasSuffix( fileOutput, ".fastq" );
    vector <FILE *> inFilesCyc;
    vector <FILE *> inFilesCycQual;
    //Open all cyc files
    for ( int i = 0; ; ++i )
    {
        Filename fn( fileInputPrefix, i, "" );
        FILE *f = fopen( fn, "rb" );
        if ( !f ) break;
        inFilesCyc.push_back( f );

        if ( outputIsFastq )
        {
            Filename fnQual( fileInputPrefix, i, ".qual" );
            inFilesCycQual.push_back( fopen( fnQual, "rb" ) );
            if ( inFilesCycQual[i] == NULL )
            {
                std::cerr << "TransposeFasta: could not open file " << fnQual << std::endl;
                exit ( EXIT_FAILURE );
            }
        }
    }
    if ( inFilesCyc.empty() )
    {
        std::cerr << "TransposeFasta: could not open file " << fileInputPrefix << "0" << std::endl;
        exit ( EXIT_FAILURE );
    }
    SequenceLength lengthRead = inFilesCyc.size();
    fseek( inFilesCyc[0], 0, SEEK_END );
    SequenceNumber nSeq = ftell( inFilesCyc[0] );
    fseek( inFilesCyc[0], 0, SEEK_SET );

    ofstream outFile ( fileOutput.c_str() );
    if ( outFile.is_open() == false )
    {
        std::cerr << "Error opening \"" << fileOutput << "\" file" << std::endl;
        exit ( 1 );
    }

    //I must read a char for each sequence. The chars at the position i corresponds to the chars of the sequence i.
    char symbol;
    string sequence = "";
    // buf to accelerate SequenceExtractor usage
    const int SEQ_EXTRACTION_BUF_SIZE = 1024;
    char seqExtractionBuf[SEQ_EXTRACTION_BUF_SIZE];
    int seqCountToSkip = 0;
    for ( SequenceNumber j = 0; j < nSeq; j++ )
    {
        bool extractThisSeq = !sequenceExtractor || sequenceExtractor->doWeExtractNextSequence();

        if ( !extractThisSeq )
        {
            ++seqCountToSkip;
            continue;
        }
        else
        {
            while ( seqCountToSkip > 0 )
            {
                size_t skip = min( seqCountToSkip, SEQ_EXTRACTION_BUF_SIZE );
                for ( SequenceLength i = 0; i < lengthRead; i++ )
                {
                    assert( fread ( seqExtractionBuf, sizeof( char ), skip, inFilesCyc[i] ) == skip );
                    if ( outputIsFastq && inFilesCycQual.size() >= lengthRead )
                        assert( fread ( seqExtractionBuf, sizeof( char ), skip, inFilesCycQual[i] ) == skip );
                }
                seqCountToSkip -= skip;
            }
        }

        if ( outputIsFastq )
            outFile << "@Read"  << j << std::endl;
        else
            outFile << "> Read "  << j << std::endl;
        for ( SequenceLength i = 0; i < lengthRead; i++ )
        {
            assert( fread ( &symbol, sizeof( char ), 1, inFilesCyc[i] ) == 1 );
            sequence.append ( 1, symbol );
        }
        outFile << sequence << std::endl;
        Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << sequence << std::endl;
        sequence.clear();

        if ( outputIsFastq )
        {
            outFile << "+" << std::endl;
            if ( outputIsFastq && inFilesCycQual.size() >= lengthRead )
            {
                for ( SequenceLength i = 0; i < lengthRead; i++ )
                {
                    assert( fread ( &symbol, sizeof( char ), 1, inFilesCycQual[i] ) == 1 );
                    sequence.append ( 1, symbol );
                }
                outFile << sequence << std::endl;
                sequence.clear();
            }
            else
                outFile << "<qualities not available>" << std::endl;
        }
    }


    outFile.close();


    //Close all cyc files
    for ( SequenceLength i = 0; i < lengthRead; i++ )
    {
        fclose( inFilesCyc[i] );
    }

    return 1;
}
