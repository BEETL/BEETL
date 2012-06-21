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

#include "TransposeFasta.h"
#include "Tools.h"
#include <stdlib.h>
#include <assert.h>
#include "SeqReader.hh"



using std::cout;
using std::cerr;
using std::endl;

TransposeFasta::TransposeFasta( SeqReaderFile* pReader ) :
  pReader_(pReader),
  cycleNum_(pReader->length()),
  outputFiles_(pReader->length()),
  buf_( pReader->length(), vector<uchar>(BUFFERSIZE))
{
  cerr << "Constructing TransposeFasta, found read length of " 
       << cycleNum_ << endl;
  for (int i(0);i<256;i++) freq[i]=0;
}



TransposeFasta::TransposeFasta()
{
  for (int i(0);i<256;i++) freq[i]=0;
}


TransposeFasta::~TransposeFasta()
{

}


bool TransposeFasta::convert( const string& input,const string& output )
{
	//TO DO
	lengthRead = cycleNum_;
	//The distribution of characters is useful
	//for alpha[256] -->Corresponding between the alphabet, the piles and tableOcc
	//and to know sizeAlpha
	//We supposed that the symbols in the input file are the following
	freq[int(TERMINATE_CHAR)]=1;
	freq[int('A')]=1;
	freq[int('C')]=1;
	freq[int('G')]=1;
	freq[int('N')]=1;
	freq[int('T')]=1;
	//GIOVANNA: ADDED THE SYMBOL Z IN THE ALPHABET, SO sizeAlpha = alphabetSize
	freq[int('Z')]=1;	

    FILE* ifile;
	ifile = fopen(input.c_str(), "rb");

	//	SeqReaderFile* pReader(SeqReaderFile::getReader(fopen(input.c_str(),"rb")));

	//	cycleNum_=pReader->getLength();
	//	cerr << "Deduced read length of " << cycleNum_ << endl;

    
    if( ifile == NULL ) { cerr << "TrasposeFasta: could not open file " << input << " !" << endl; }

    // create output files
    for(dataTypelenSeq i=0;i<cycleNum_;i++ )
    {
        std::stringstream fn;
        fn << output <<  i << ".txt";
        outputFiles_[i] = fopen( fn.str().c_str(),"w" );
    }


    // looping through the input file, add the characters to the buffer, print buffer when it's full
//    unsigned int num_read = 0;
    unsigned int num_write = 0;    
    unsigned int charsBuffered = 0;

    	//******************************buf[cycleNum_+1];  ********* is cycleNum_ right?
	char buf[cycleNum_+1];
	lengthTexts = 0;


    for(dataTypelenSeq i=0;i<cycleNum_;i++ )
    {
        buf[i] = '\0';
    }

    nSeq = 0;
//    num_read = fread(buf,sizeof(uchar),cycleNum_,ifile);

//    fgets ( buf,1024, ifile ); %%%%%
//    while( !feof(ifile) ) %%%%%
    while( pReader_->allRead()==false )
    {
        //cerr << "current line : " << buf << endl;
        
        if( charsBuffered == BUFFERSIZE )
        {
            // write buffers to the files, clear buffers
            for(dataTypelenSeq i=0;i<cycleNum_;i++ )
            {
                //cerr << "writing to " << i << " : " << buf_[i] << endl;
	      num_write = fwrite ( (void*)(&buf_[i][0]),sizeof(char),charsBuffered,outputFiles_[i] );
				lengthTexts += num_write;
            }
            checkIfEqual( num_write,charsBuffered ); // we should always read/write the same number of characters


            charsBuffered=0;
        }

	for(dataTypelenSeq i=0;i<cycleNum_;i++ )
	{
	  buf_[i][charsBuffered] = pReader_->thisSeq()[i];
            // increase the counter of chars buffered
	}
	  charsBuffered++;        
	  nSeq++;


#ifdef XXX
        // process the input
        if( buf[0] != '>' ) 
        {
            // add the characters
            for(dataTypelenSeq i=0;i<cycleNum_;i++ )
            {
                buf_[i][charsBuffered] = buf[i];
            }
            
            // increase the counter of chars buffered
            charsBuffered++;        
	    			nSeq++;
        }
#endif		//else 


        //num_read = fread(buf,sizeof(uchar),cycleNum_,ifile);        
	//        fgets ( buf, 1024, ifile );
	pReader_->readNext();
    }

    // write the rest
    for(dataTypelenSeq i=0;i<cycleNum_;i++ )
    {
      num_write = fwrite ( (void*)(&buf_[i][0]),sizeof(uchar),charsBuffered,outputFiles_[i] );
		lengthTexts += num_write;
    }
    checkIfEqual( num_write,charsBuffered );



    // closing all the output file streams
    for(dataTypelenSeq i=0;i<cycleNum_;i++ )
    {
        fclose( outputFiles_[i] );
    }

    std::cout << "Number of sequences reading/writing: " << nSeq << "\n";
    std::cout << "Number of characters reading/writing: " << lengthTexts << "\n";

    //    delete pReader;
    return true;
}

bool TransposeFasta::inputCycFile() {
  		//TO DO
	//lengthRead = cycleNum_;
	//The distribution of characters is useful
	//for alpha[256] -->Corresponding between the alphabet, the piles and tableOcc
	//and to know sizeAlpha

	//TDB

	// Setting useful to work:
	//1) Alphabet
	//We supposed that the symbols in the input file are the following
	freq[int(TERMINATE_CHAR)]=1;
	freq[int('A')]=1;
	freq[int('C')]=1;
	freq[int('G')]=1;
	freq[int('N')]=1;
	freq[int('T')]=1;
	//GIOVANNA: ADDED THE SYMBOL Z IN THE ALPHABET, SO sizeAlpha = alphabetSize
	freq[int('Z')]=1;	

	//2) Number of sequences
	nSeq = 100;
	//3) Length of the longest sequence
	lengthRead = 100;

	std::cerr << "****number of sequences: " << nSeq << "\n";
	std::cerr << "****max length of each sequence: " << lengthRead << "\n";
	//4) Total Length 
	lengthTexts = lengthRead * nSeq;
	std::cerr << "****lengthTot: " << lengthTexts << "\n";
    return 1;
} 

bool TransposeFasta::convertFromCycFileToFasta(const string& fileOutput, dataTypeNSeq nSeq, dataTypelenSeq lengthRead) {	
	vector <FILE*> inFilesCyc_; 
	inFilesCyc_.resize(lengthRead);    //One for each symbol of the read.
    const char* fileOut = "cyc.";
	//Open all cyc files
	for(dataTypelenSeq i=0;i<lengthRead;i++ ) {
        std::stringstream fn;
        fn << fileOut <<  (int)i << ".txt";
        inFilesCyc_[i] = fopen( fn.str().c_str(),"rb" );
        if (inFilesCyc_[i] == NULL) {
                std::cerr << "TrasposeFasta: could not open file "  <<  fn.str().c_str() << std::endl;
				exit (EXIT_FAILURE);
		}
    }
//	FILE *outFile = fopen(fileOutput.c_str(), "rb");
//	if (outFile==NULL) 

	ofstream outFile (fileOutput.c_str());
	if (outFile.is_open() == false)	{
		 std::cerr << "Error opening \"" << fileOutput << "\" file"<< std::endl;
		exit (1);
	}
	
	//I must read a char for each sequence. The chars at the position i corresponds to the chars of the sequence i.
	dataTypeNChar num_read;
	char symbol;
	string sequence = "";
	for(dataTypeNSeq j=0;j<nSeq;j++ ) {
		outFile << "> Read "  << j << std::endl;
		for(dataTypelenSeq i=0;i<lengthRead;i++ ) {
			num_read = fread (&symbol, sizeof(char), 1, inFilesCyc_[i] );
			sequence.append (1, symbol);
		}
		outFile << sequence << std::endl;
		if (verboseDecode == 1) 
			std::cerr << sequence << std::endl;
		sequence = "";
	}

	
	outFile.close();
		

	//Close all cyc files
	for(dataTypelenSeq i=0;i<lengthRead;i++){
        fclose(inFilesCyc_[i]);
    }

	return 1;
}
