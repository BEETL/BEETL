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

#include "TransposeFasta.hh"

#include "Logger.hh"
#include "SeqReader.hh"
#include "TemporaryFilesManager.hh"
#include "Tools.hh"

#include <cassert>
#include <cstdlib>

using namespace std;


TransposeFasta::TransposeFasta()
{
    for ( int i( 0 ); i<256; i++ ) freq[i]=0;
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


bool TransposeFasta::convert( const string &input, const string &output )
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
    freq[int( TERMINATE_CHAR )]=1;
    freq[int( 'A' )]=1;
    freq[int( 'C' )]=1;
    freq[int( 'G' )]=1;
    freq[int( 'N' )]=1;
    freq[int( 'T' )]=1;
    //GIOVANNA: ADDED THE SYMBOL Z IN THE ALPHABET, SO sizeAlpha = alphabetSize
    freq[int( 'Z' )]=1;

    FILE *ifile;
    ifile = fopen( input.c_str(), "rb" );

    //        SeqReaderFile* pReader(SeqReaderFile::getReader(fopen(input.c_str(),"rb")));

    //        cycleNum_=pReader->getLength();
    //        cerr << "Deduced read length of " << cycleNum_ << endl;


    if ( ifile == NULL )
    {
        cerr << "TransposeFasta: could not open file " << input << " !" << endl;
    }

    // create output files
    for ( dataTypelenSeq i=0; i<cycleNum_; i++ )
    {
        stringstream fn;
        fn << output <<  i << ".txt";
        outputFiles_[i] = fopen( fn.str().c_str(),"w" );
        TemporaryFilesManager::get().addFilename( fn.str() );
        if ( processQualities_ )
        {
            stringstream fnQual;
            fnQual << output << "qual." << i << ".txt";
            outputFilesQual[i] = fopen( fnQual.str().c_str(),"w" );
            TemporaryFilesManager::get().addFilename( fnQual.str() );
        }
    }


    // looping through the input file, add the characters to the buffer, print buffer when it's full
    //    unsigned int num_read = 0;
    unsigned int num_write = 0;
    unsigned int charsBuffered = 0;

    //******************************buf[cycleNum_+1];  ********* is cycleNum_ right?
    char buf[cycleNum_+1];
    lengthTexts = 0;


    for ( dataTypelenSeq i=0; i<cycleNum_; i++ )
    {
        buf[i] = '\0';
    }

    nSeq = 0;
    //    num_read = fread(buf,sizeof(uchar),cycleNum_,ifile);

    //    fgets ( buf,1024, ifile ); %%%%%
    //    while( !feof(ifile) ) %%%%%
    while ( pReader_->allRead()==false )
    {
        //cerr << "current line : " << buf << endl;

        if ( charsBuffered == BUFFERSIZE )
        {
            // write buffers to the files, clear buffers
            for ( dataTypelenSeq i=0; i<cycleNum_; i++ )
            {
                //cerr << "writing to " << i << " : " << buf_[i] << endl;
                num_write = fwrite ( ( void * )( &buf_[i][0] ),sizeof( char ),charsBuffered,outputFiles_[i] );
                lengthTexts += num_write;
                if ( processQualities_ )
                {
                    size_t num_write_qual = fwrite ( ( void * )( &bufQual[i][0] ),sizeof( char ),charsBuffered,outputFilesQual[i] );
                    checkIfEqual( num_write, num_write_qual );
                }
            }
            checkIfEqual( num_write,charsBuffered ); // we should always read/write the same number of characters


            charsBuffered=0;
        }

        for ( dataTypelenSeq i=0; i<cycleNum_; i++ )
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
            for ( dataTypelenSeq i=0; i<cycleNum_; i++ )
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
    for ( dataTypelenSeq i=0; i<cycleNum_; i++ )
    {
        num_write = fwrite ( ( void * )( &buf_[i][0] ),sizeof( uchar ),charsBuffered,outputFiles_[i] );
        lengthTexts += num_write;
        if ( processQualities_ )
        {
            size_t num_write_qual = fwrite ( ( void * )( &bufQual[i][0] ),sizeof( uchar ),charsBuffered,outputFilesQual[i] );
            checkIfEqual( num_write, num_write_qual );
        }
    }
    checkIfEqual( num_write,charsBuffered );



    // closing all the output file streams
    for ( dataTypelenSeq i=0; i<cycleNum_; i++ )
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
    freq[int( TERMINATE_CHAR )]=1;
    freq[int( 'A' )]=1;
    freq[int( 'C' )]=1;
    freq[int( 'G' )]=1;
    freq[int( 'N' )]=1;
    freq[int( 'T' )]=1;
    //GIOVANNA: ADDED THE SYMBOL Z IN THE ALPHABET, SO sizeAlpha = alphabetSize
    freq[int( 'Z' )]=1;

    //2) Number of sequences
    string cyc1Filename = cycPrefix + "1.txt";
    FILE *f = fopen( cyc1Filename.c_str(), "rb" );
    if ( !f )
    {
        cerr << "ERROR: Cycle file " << cyc1Filename << " not found!" << endl;
        exit( -1 );
    }
    fseek( f, 0, SEEK_END );
    nSeq = ftell( f );
    fclose( f );

    //3) Length of the longest sequence
    for ( lengthRead=1; ; ++lengthRead )
    {
        stringstream cycFilename;
        cycFilename << cycPrefix << lengthRead << ".txt";
        FILE *f = fopen( cycFilename.str().c_str(), "rb" );
        if ( f )
            fclose( f );
        else
            break;
    }

    //4) qualities detection
    string qual1Filename = cycPrefix + "qual.1.txt";
    f = fopen( qual1Filename.c_str(), "rb" );
    if ( f )
    {
        processQualities_ = true;
        fclose( f );
    }
    else
        processQualities_ = false;
    Logger::out( LOG_SHOW_IF_VERBOSE ) << "****processing qualities: " << processQualities_ << "\n";


    Logger::out( LOG_SHOW_IF_VERBOSE ) << "****number of sequences: " << nSeq << "\n";
    Logger::out( LOG_SHOW_IF_VERBOSE ) << "****max length of each sequence: " << lengthRead << "\n";

    //5) Total Length
    lengthTexts = lengthRead * nSeq;
    Logger::out( LOG_SHOW_IF_VERBOSE ) << "****lengthTot: " << lengthTexts << "\n";
    return 1;
}


bool TransposeFasta::convertFromCycFileToFastaOrFastq( const string &fileInputPrefix, const string &fileOutput )
{
    bool outputIsFastq = hasSuffix( fileOutput, ".fastq" );
    vector <FILE *> inFilesCyc;
    vector <FILE *> inFilesCycQual;
    //Open all cyc files
    for ( int i=0; ; ++i )
    {
        stringstream fn;
        fn << fileInputPrefix <<  ( int )i << ".txt";
        FILE *f = fopen( fn.str().c_str(),"rb" );
        if ( !f ) break;
        inFilesCyc.push_back( f );

        if ( outputIsFastq )
        {
            stringstream fnQual;
            fnQual << fileInputPrefix <<  ( int )i << ".qual.txt";
            inFilesCycQual.push_back( fopen( fnQual.str().c_str(),"rb" ) );
            if ( inFilesCycQual[i] == NULL )
            {
                std::cerr << "TransposeFasta: could not open file "  <<  fn.str().c_str() << std::endl;
                exit ( EXIT_FAILURE );
            }
        }
    }
    if ( inFilesCyc.empty() )
    {
        std::cerr << "TransposeFasta: could not open file " << fileInputPrefix << "0.txt" << std::endl;
        exit ( EXIT_FAILURE );
    }
    dataTypelenSeq lengthRead = inFilesCyc.size();
    fseek( inFilesCyc[0], 0, SEEK_END );
    dataTypeNSeq nSeq = ftell( inFilesCyc[0] );
    fseek( inFilesCyc[0], 0, SEEK_SET );

    ofstream outFile ( fileOutput.c_str() );
    if ( outFile.is_open() == false )
    {
        std::cerr << "Error opening \"" << fileOutput << "\" file"<< std::endl;
        exit ( 1 );
    }

    //I must read a char for each sequence. The chars at the position i corresponds to the chars of the sequence i.
    dataTypeNChar num_read;
    char symbol;
    string sequence = "";
    for ( dataTypeNSeq j=0; j<nSeq; j++ )
    {
        if ( outputIsFastq )
            outFile << "@Read"  << j << std::endl;
        else
            outFile << "> Read "  << j << std::endl;
        for ( dataTypelenSeq i=0; i<lengthRead; i++ )
        {
            num_read = fread ( &symbol, sizeof( char ), 1, inFilesCyc[i] );
            sequence.append ( 1, symbol );
        }
        outFile << sequence << std::endl;
        Logger::out( LOG_FOR_DEBUGGING ) << sequence << std::endl;
        sequence.clear();

        if ( outputIsFastq )
        {
            outFile << "+" << std::endl;
            if ( outputIsFastq && inFilesCycQual.size() >= lengthRead )
            {
                for ( dataTypelenSeq i=0; i<lengthRead; i++ )
                {
                    num_read = fread ( &symbol, sizeof( char ), 1, inFilesCycQual[i] );
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
    for ( dataTypelenSeq i=0; i<lengthRead; i++ )
    {
        fclose( inFilesCyc[i] );
    }

    return 1;
}
