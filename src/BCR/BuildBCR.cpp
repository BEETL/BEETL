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

#include <string.h>
#include <sys/stat.h>
#include <stdlib.h>
#include "BCRexternalBWT.h"
#include "BWTCollection.h"
#include "TransposeFasta.h"
#include "LetterCount.hh"
#include "BwtReader.hh"
#include "BwtWriter.hh"



//#ifdef COMPRESS_BWT
//typedef BwtReaderRunLength BwtReader;
//typedef BwtWriterRunLength BwtWriter;
//#else
//typedef BwtReaderASCII BwtReader;
//typedef BwtWriterASCII BwtWriter;
//#endif


#define SIZEBUFFER 1024 // TBD get rid of this



int BCRexternalBWT::buildBCR(char const * file1, char const * fileOut)
{
  TransposeFasta trasp;
  int controllo = trasp.convert( file1, fileOut);
	
  if (controllo == false) {  //Error in the reading
    std::cerr << "Reading Error \n";
    exit (EXIT_FAILURE);
  }
	
  nText = trasp.nSeq;
  lengthRead = trasp.lengthRead;
  lengthTot = trasp.lengthTexts;

  //#ifdef REPLACE_TABLEOCC
  sizeAlpha=0;
  for (dataTypedimAlpha i = 0; i < 255; ++i)
    if (trasp.freq[i] > 0) {
      alpha[i] = sizeAlpha;
      sizeAlpha++;
    }
  std::cerr << "We supposed that the symbols in the input file are:\n";
  for (dataTypedimAlpha i = 0; i < 255; ++i)
    if (trasp.freq[i] > 0) 
      std::cerr << i << " " << trasp.freq[i] << " " << (int)alpha[i] << "\n";
  //#endif

  lengthTot_plus_eof = lengthTot+nText;

  std::cerr << "sizeof(type size of alpha): " << sizeof(dataTypedimAlpha) << "\n";
  std::cerr << "sizeof(type of #sequences): " << sizeof(dataTypeNSeq) << "\n";
  std::cerr << "sizeof(type of #characters): " << sizeof(dataTypeNChar) << "\n";
	
  std::cerr << "\nNumber of sequences: " << nText << "\n";
  std::cerr << "Length of each sequence: " << lengthRead << "\n\n";
  std::cerr << "Total length (without $): " << lengthTot << "\n";
  std::cerr << "Total length (with $): " << lengthTot_plus_eof << "\n";

  dataTypeNChar numchar;
  char *filename;
  filename = new char[strlen(fileOut)+sizeof(dataTypelenSeq)*8];
	
  std::cerr << "Partial File name for input: " << fileOut <<" \n\n";

  static FILE *InFileInputText;
	
  uchar *newSymb = new uchar[nText];
  vectTriple.resize(nText); 

#ifdef REPLACE_TABLEOCC
  tableOcc = new dataTypeNChar*[sizeAlpha]; 
  for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {   //Counting for each pile: $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
    tableOcc[j] = new dataTypeNChar[sizeAlpha];
  }
  for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) 
    for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++) 
      tableOcc[j][h]=0;
#endif
  tableOcc_.clear();

  //lengthTot = 0;  //Counts the number of symbols

  std::cerr << "\nFirst symbols: "<< "j= "<< 1 <<" - symbols in position " << lengthRead-1 << "\n";

  numchar=sprintf (filename, "%s%u.txt", fileOut, lengthRead-1);
  InFileInputText = fopen(filename, "rb");
  // fileOut opened?
  if (InFileInputText==NULL) {
    // NO, abort program
    std::cerr << filename <<" : Error opening " << std::endl;
    exit (EXIT_FAILURE);
  }
  dataTypeNChar num = fread(newSymb,sizeof(uchar),nText,InFileInputText);
  checkIfEqual( num,nText); // we should always read the same number of characters
  //lengthTot += num;   //Increment the number of chars
  fclose(InFileInputText);
  /*
    std::cerr << filename  << std::endl;
    std::cerr << "number read " << num << "\n";
    for (dataTypeNSeq j = 0 ; j < nText; j++) 
    std::cerr << newSymb[j];
    std::cerr << ".\n";
  */
  InsertFirstsymbols(newSymb);  
  
  //maxLengthRead-2
  for (int t = lengthRead-2 ; t >= 0; t--) {
    std::cerr << "j= "<< lengthRead - t <<" - symbols in position " << t << "\n";
    //cout << "Starting iteration " << t << ", time now: " << timer.timeNow();
    //cout << "Starting iteration " << t << ", usage: " << timer << endl;

    //To insert the symbol from position m-3 to position 1
    //The last inserted symbol is in position i+1 (or it is newSymb[j]), 
    //the next symbol (to insert) is in position i
		
    numchar=sprintf (filename, "%s%u.txt", fileOut, t);
    InFileInputText = fopen(filename, "rb");
    // fileOut opened?
    if (InFileInputText==NULL) {
      // NO, abort program
      std::cerr << filename <<" : Error opening " << std::endl;
      exit (EXIT_FAILURE);
    }
    num = fread(newSymb,sizeof(uchar),nText,InFileInputText);
    checkIfEqual( num,nText); // we should always read the same number of characters
		
    /*
      std::cerr << filename  << std::endl;
      std::cerr << "read number is " << num << "\n";
      for (dataTypeNSeq j = 0 ; j < nText; j++) 
      std::cerr << newSymb[j];
      std::cerr << ".\n";
    */
    //lengthTot += num;   //Increment the number of chars
    fclose(InFileInputText);
    InsertNsymbols(newSymb, t);
    //    exit (0); // %%%%%
  }

  //The last inserted symbol is in position 0 (or it is newSymb[j]), 
  //the next symbol (to insert) is in position m-1, that is, I have to inserted the symbols $ 
  std::cerr << "\n*******Inserting $ symbol\n\n";
  //cout << "Starting iteration " << lengthRead-1 << ", time now: " << timer.timeNow();
  //cout << "Starting iteration " << lengthRead-1 << ", usage: " << timer << endl;
  for (dataTypeNSeq j = 0 ; j < nText; j++) {
    newSymb[j] = '$';
    //std::cerr << newSymb[j];
  }
  //std::cerr << "\n";

  InsertNsymbols(newSymb, lengthRead); 

  //cout << "All new characters inserted, usage: " << timer << endl;


  // to delete those:
  delete [] newSymb;
  //	vectTriple.~vector<sortElement>();
  /*
    delete [] seqN;
    delete [] pileN;
    delete [] posN;
  */
  //  std::cerr << std::endl;
  //std::cerr << "The input text is long " << lengthTot << std::endl;

  dataTypeNChar numCharInTable = 0;
  for (dataTypedimAlpha r = 0; r < alphabetSize; r++) {
    for (dataTypedimAlpha t = 0; t < alphabetSize; t++) {
      numCharInTable += tableOcc_[r].count_[t];
    }
  } 
	std::cerr << "In tableOcc_, there are " << numCharInTable << " letters" << std::endl;

#ifdef REPLACING_TABLEOCC
	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) { 
	  delete [] tableOcc[j];
	  tableOcc[j] = NULL;
	}
	delete [] tableOcc;
#endif
	
	return true;
} // ~buildBCR

void BCRexternalBWT::InsertFirstsymbols(uchar const * newSymb)
{
	for (dataTypeNSeq j = 0 ; j < nText; j++) {
		vectTriple[j].posN = j+1;  // position of the suffix (1-based)
		vectTriple[j].seqN = j;	  // number of sequence
		vectTriple[j].pileN = 0;    //The first symbols are in $-pile
	}
	std::cerr << "First step" << std::endl;
	if (verboseEncode == 1) {
		std::cerr << "Q  ";
		for (dataTypeNSeq g = 0 ; g < nText; g++) {  
			std::cerr << (int)vectTriple[g].pileN << " ";
		}
		std::cerr << std::endl;
		std::cerr << "P  ";
		for (dataTypeNSeq g = 0 ; g < nText; g++) {  
			std::cerr << vectTriple[g].posN  << " ";
		}
		std::cerr << std::endl;
		std::cerr << "N  ";
		for (dataTypeNSeq g = 0 ; g < nText; g++) {  
			std::cerr << vectTriple[g].seqN  << " ";
		}
		std::cerr << std::endl;
	}


	static FILE *OutFileBWT;                  // output and input file BWT;
	char *filenameOut = new char[12];
	char *filename = new char[8];
	dataTypeNChar numchar;
	const char *ext = ".aux";
	numchar=sprintf (filename, "bwt_%d",0);
	
	numchar=sprintf (filenameOut,"%s%s",filename,ext);
	
	OutFileBWT = fopen(filenameOut, "w");
	// fileOut aperto?
	if (OutFileBWT==NULL) {
	// NO, abort program
		std::cerr << "BWT file $: Error opening " << std::endl;
		exit (EXIT_FAILURE);
	}
	for (dataTypeNSeq j = 0 ; j < nText; j++) {
		tableOcc_[0].count_[whichPile[(int)newSymb[j]]]++;       //counting the number of occurrences in BWT of the $-pile
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
	dataTypeNChar num = fwrite (newSymb, sizeof(uchar), nText , OutFileBWT); 
        checkIfEqual( num,nText);// we should always read the same number of characters

	//if (num != nText)
	//	std::cerr << "the written characters is not equal to number of the texts" << num << " and "<< nText <<"\n";
	fclose(OutFileBWT);

	//Creates one file for each letter in the alphabet. From 1 to alphabetSize-1
	for (dataTypedimAlpha i = 1; i < alphabetSize-1; i++) {
		numchar=sprintf (filename, "bwt_%d", i);
		numchar=sprintf (filenameOut,"%s%s",filename,ext);
		
		OutFileBWT = fopen(filenameOut, "w");
		// fileOut aperto?
		if (OutFileBWT==NULL) {
		// NO, abort program
			std::cerr << "BWT file " << (int)i <<" : Error opening " << std::endl;
			exit (EXIT_FAILURE);
		}
		fclose(OutFileBWT);
	}

	delete [] filenameOut;
	delete [] filename;
}


void BCRexternalBWT::InsertNsymbols(uchar const * newSymb, dataTypelenSeq posSymb)
{
  #ifdef DEBUG
  ulong zz(0); // TEMP
  #endif

  //<<<<<<< BuildBCR.cpp
  //they are not first symbols
  static FILE *InFileBWT;                  // output and input file BWT;
  char *filenameIn = new char[12];
  char *filename = new char[8];
  const char *ext = ".aux";
	
  //std::cerr << "Compute new posN" << std::endl;

  LetterCount counters;
  
  dataTypeNChar numchar=0;
  dataTypeNChar toRead = 0;
  //Find the positions of the new symbols
  dataTypeNSeq j = 0;

  BwtReaderBase* pReader(NULL);

  counters.clear();
	
  while (j < nText) {			
    dataTypedimAlpha currentPile = vectTriple[j].pileN;
    assert (currentPile<alphabetSize-1);
    numchar=sprintf (filename, "bwt_%d", currentPile);
    numchar=sprintf (filenameIn,"%s%s",filename,ext);
    //printf("===Current BWT-partial= %d\n",currentPile);
    if (pReader!=NULL) delete pReader;
    //    pReader= new BwtReader(filenameIn);
    if ((outputCompression_==compressionASCII)||(posSymb==98))
      pReader = new BwtReaderASCII(filenameIn);
    else if (outputCompression_==compressionRunLength)
      pReader = new BwtReaderRunLength(filenameIn);
    else if (outputCompression_==compressionHuffman)
      pReader = new BwtReaderHuffman(filenameIn);
    else 
    { 
      cerr << "!! BuildBcr.cpp - unknown compression type specified: " 
	   << outputCompression_  << endl;
      exit (-1);
    } // ~else
    assert(pReader!=NULL);


    //    BwtReader reader(filenameIn);

    dataTypeNSeq k=j;
    //For each pile, we have a different counter of characters
    for (dataTypedimAlpha i = 0 ; i < alphabetSize; i++)
      counters.count_[i]=0;
    dataTypeNChar cont = 0;   //number of the read symbols
    uchar foundSymbol;
    dataTypeNChar numberRead=0;
    while ((k< nText) && (vectTriple[k].pileN == currentPile)) {
      if (verboseEncode == 1)
	std::cerr << "j-1: Q["<<k<<"]=" << (int)vectTriple[k].pileN << " P["<<k<<"]=" << (dataTypeNChar)vectTriple[k].posN << " N["<<k<<"]=" << (dataTypeNSeq)vectTriple[k].seqN << "\t";
      
      //std::cerr << "--k= " << k << " pileN[k]= " << pileN[k] << " posN[k]= " << posN[k] <<  " seqN[k]= " << seqN[k] << std::endl;
      //For any character (of differents sequences) in the same pile
      foundSymbol = '\0';

      //cont is the number of symbols already read!
      //      toRead = vectTriple[k].posN - cont -1;

      

      toRead = vectTriple[k].posN - cont;
      //      cout << "toRead: " << toRead << endl;
      if (toRead>0)
      {
	if (toRead>1)
	{
	  numberRead = (*pReader).readAndCount( counters, toRead-1 );
	  assert (toRead-1 == numberRead);
	}
	assert((*pReader)( (char*)&foundSymbol, 1)==1);
	if(whichPile[(int)foundSymbol]<alphabetSize-1) {}
	else { cout << (int)foundSymbol << " " << foundSymbol << endl; assert(1==0); }

	counters.count_[whichPile[(int)foundSymbol]]++;
	cont += toRead;
      }

#ifdef DEBUG
      cout << zz << " " << toRead << " " << foundSymbol  << " "; counters.print();
#endif

      //std::cerr << "toRead " << toRead << "Found Symbol is " << foundSymbol << "\n";


      //I have to update the value in vectTriple[k].posN, it must contain the position of the new symbol
      //#ifdef XXX
      vectTriple[k].posN = counters.count_[whichPile[(int)foundSymbol]];   
      //#endif
      //std::cerr << "--New posN[k]=" << (int)posN[k] <<std::endl;
      if (verboseEncode == 1)				
	std::cerr << "\nInit New P["<< k <<"]= " << vectTriple[k].posN <<std::endl;

      for (dataTypedimAlpha g = 0 ; g < currentPile; g++) {  //I have to count in each pile g= 0... (currentPile-1)-pile
	vectTriple[k].posN = vectTriple[k].posN + tableOcc_[g].count_[whichPile[(int)foundSymbol]];
	//std::cerr << "--New posN[k]=" << (int)posN[k] << " tableOcc[g][whichPile[(int)symbol]] " << tableOcc[g][whichPile[(int)symbol]] <<std::endl;
	if (verboseEncode == 1) {				
	  std::cerr << "g= " << (int)g << " symbol= " << (int)foundSymbol << " whichPile[symbol]= "
		    << (int)whichPile[(int)foundSymbol] <<std::endl;
	  std::cerr << "Add New posN[k]=" << vectTriple[k].posN << " tableOcc[g][whichPile[(int)symbol]] " 
		    << tableOcc_[g].count_[whichPile[(int)foundSymbol]] <<std::endl;
	}
	
      }
      //I have to insert the new symbol in the symbol-pile
      assert(whichPile[(int)foundSymbol]<alphabetSize-1);
      vectTriple[k].pileN=whichPile[(int)foundSymbol];
      //std::cerr << "New posN[k]=" << (int)posN[k] << " New pileN[k]=" << (int)pileN[k] << std::endl;
      if (verboseEncode == 1)
	std::cerr << "j  : Q[q]=" << (int)vectTriple[k].pileN << " P[q]=" << (dataTypeNChar)vectTriple[k].posN <<  " N[q]=" << (dataTypeNSeq)vectTriple[k].seqN << std::endl;
      
      k++;
    }

    //    fclose(InFileBWT);
    j=k;
  } // ~while j
  
  if (verboseEncode==1) {
    uchar *buffer = new uchar[SIZEBUFFER];
    dataTypedimAlpha mmm=0;
    while (mmm < alphabetSize) {			
      numchar=sprintf (filename, "bwt_%d", mmm);
      numchar=sprintf (filenameIn,"%s%s",filename,ext);
      //printf("===currentPile= %d\n",mmm);
      InFileBWT = fopen(filenameIn, "r");
      for (dataTypeNSeq g = 0 ; g < SIZEBUFFER; g++) 
	buffer[g] = '\0';
      numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,InFileBWT);
      std::cerr << "B[" << mmm << "]:\t";
      if (numchar==0)
	std::cerr  << "empty\n";
      else
	std::cerr  << buffer << "\n";
      fclose(InFileBWT);
			mmm++;
    }
    delete [] buffer;
  } // ~if verboseEncode


#ifdef XXX
  delete [] counters;
#endif
  delete [] filenameIn;
  delete [] filename;
	
  if (verboseEncode==1) {
    std::cerr << "NewSymbols " ;
    for (dataTypeNSeq g = 0 ; g < nText; g++) {  
      std::cerr << newSymb[g] << " ";
    }
    std::cerr << std::endl;
    std::cerr << "Before Sorting" << std::endl;
    std::cerr << "Q  ";
    for (dataTypeNSeq g = 0 ; g < nText; g++) {  
      std::cerr << (int)vectTriple[g].pileN << " ";
    }
    std::cerr << std::endl;
    std::cerr << "P  ";
    for (dataTypeNSeq g = 0 ; g < nText; g++) {  
      std::cerr << vectTriple[g].posN  << " ";
    }
    std::cerr << std::endl;
    std::cerr << "N  ";
    for (dataTypeNSeq g = 0 ; g < nText; g++) {  
      std::cerr << vectTriple[g].seqN  << " ";
    }
    std::cerr << std::endl;
  } // ~if verboseEncode
	
  quickSort(vectTriple);
		
  if (verboseEncode==1) {
    std::cerr << "After Sorting" << std::endl;
    std::cerr << "Q  ";
    for (dataTypeNSeq g = 0 ; g < nText; g++) {  
      std::cerr << (int)vectTriple[g].pileN << " ";
    }
    std::cerr << std::endl;
    std::cerr << "P  ";
    for (dataTypeNSeq g = 0 ; g < nText; g++) {  
      std::cerr << vectTriple[g].posN  << " ";
    }
    std::cerr << std::endl;
    std::cerr << "N  ";
    for (dataTypeNSeq g = 0 ; g < nText; g++) {  
      std::cerr << vectTriple[g].seqN  << " ";
    }	
    std::cerr << std::endl;
  } // ~if verboseEncode

  storeBWT(newSymb);
  //std::cerr << "End storing BWT" << std::endl;
	
  if (posSymb == lengthRead) { //We are storing the last simbols of sequences: $
    std::cerr << "Stores the 'end positions' of the $!"<< std::endl;
    const char *fileEndPos="outFileEndPos.bwt";
    static FILE *OutFileEndPos;                  // output file of the end positions;
    OutFileEndPos = fopen(fileEndPos, "w");
    // fileOut aperto?
    if (OutFileEndPos==NULL) {
      // NO, abort program
      std::cerr << "Error opening \"" << fileEndPos << "\" file"<< std::endl;
      exit (EXIT_FAILURE);
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
		
    //<<<<<<< BuildBCR.cpp
    numchar = fwrite (&nText, sizeof(dataTypeNChar), 1 , OutFileEndPos);
    assert( numchar == 1); // we should always read the same number of characters
    
    for (dataTypeNSeq i = 0; i < nText; i++) {
      if (verboseEncode == 1)
	std::cerr << "Triple: " << vectTriple[i].seqN << " " << vectTriple[i].posN << " " << (int)vectTriple[i].pileN << "\n" << std::endl;
      numchar = fwrite (&vectTriple[i].seqN, sizeof(dataTypeNSeq), 1 , OutFileEndPos); 
      assert( numchar == 1); // we should always read the same number of characters
      numchar = fwrite (&vectTriple[i].posN, sizeof(dataTypeNChar), 1 , OutFileEndPos); //here vectTriple[i].posN is the relative position of $ in the partial BWT
      assert( numchar == 1); // we should always read the same number of characters
      numchar = fwrite (&vectTriple[i].pileN, sizeof(dataTypedimAlpha), 1 , OutFileEndPos); 
      assert( numchar == 1); // we should always read the same number of characters
    }
    
    fclose(OutFileEndPos);		
    std::cerr << "The'end positions' stored!"<< std::endl;
  } // ~if (posSymb==lengthRead)

  //  delete pReader;
} // ~InsertNsymbols

void BCRexternalBWT::storeBWT(uchar const * newSymb) {
	
  //I have found the position where I have to insert the chars in the position t of the each text
  //Now I have to update the BWT in each file.
  static FILE *OutFileBWT;                  // output and input file BWT;
  char *filenameOut = new char[16]; 
  char *filenameIn = new char[12];
  char *filename = new char[8];
  const char *ext = ".aux";
	
  dataTypeNChar numchar=0;
  //  dataTypeNChar numcharWrite=0;
#ifdef OLD
  uchar *buffer = new uchar[SIZEBUFFER];
#endif
  dataTypeNChar toRead = 0;
  
  dataTypeNSeq j;
  dataTypedimAlpha currentPile;
  uchar symbol='\0';

  BwtReaderBase* pReader(NULL);
  BwtWriterBase* pWriter(NULL);


  j=0;
  while (j < nText) {	
    currentPile = vectTriple[j].pileN;
    if (verboseEncode==1)		
      std::cerr << "index j= " << j << " current BWT segment " << (int)currentPile << std::endl;

    //std::cerr << "Pile " << (int)currentPile << std::endl;
    numchar=sprintf (filename, "bwt_%d", currentPile);
    numchar=sprintf (filenameIn,"%s%s",filename,ext);
    numchar=sprintf (filenameOut,"new_%s%s",filename,ext);

    if (pReader!=NULL) delete pReader;
    if (pWriter!=NULL) delete pWriter;

    if (currentPile==0)
    {
      pReader = new BwtReaderASCII(filenameIn);
      pWriter = new BwtWriterASCII(filenameOut);
    }
    else
    {
    if (outputCompression_==compressionASCII)
    {
      // if (currentPile!=0)
	pReader = new BwtReaderASCII(filenameIn);
      pWriter = new BwtWriterASCII(filenameOut);
    }
    else if (outputCompression_==compressionRunLength)
    {
      //      if (currentPile!=0)
	pReader = new BwtReaderRunLength(filenameIn);
      pWriter = new BwtWriterRunLength(filenameOut);

    }
    else if (outputCompression_==compressionHuffman)
    {
      //      if (currentPile!=0)
	pReader = new BwtReaderHuffman(filenameIn);
      pWriter = new BwtWriterHuffman(filenameOut);

    }
    else 
    { 
      cerr << "!! BuildBcr.cpp - unknown compression type specified: " 
	   << outputCompression_  << endl;
      exit (-1);
    }
    }

    //    pReader=new BwtReader(filenameIn);
    //  pWriter= new BwtWriter(filenameOut);

    assert(pReader!=NULL);
    assert(pWriter!=NULL);


    //For each new symbol in the same pile
    dataTypeNSeq k=j;
    dataTypeNChar cont = 0;
    while ((k< nText) && (vectTriple[k].pileN == currentPile)) {
      if (verboseEncode==1)
	std::cerr << "k= " << k << " Q[k]= " << (int)vectTriple[k].pileN << " P[k]= " << vectTriple[k].posN << " cont = "<< cont << std::endl;
      //std::cerr << "k= " << k << " pileN[k]= " << pileN[k] << " posN[k]= " << posN[k] << std::endl;
      //So I have to read the k-BWT and I have to count the number of the symbols up to the position posN.
      symbol = '\0';
      //As PosN starts to the position 1 and I have to insert the new symbol in position posN[k]
      // I have to read posN[k]-1 symbols
      //cont is the number of symbols already read!
      toRead = (vectTriple[k].posN-1) - cont;
      if (verboseEncode == 1)					
	std::cerr << "Start: to Read " << toRead << "\n";

      (*pReader).readAndSend((*pWriter), toRead);

      cont+=toRead;

      //      toRead=0;


      //Now I have to insert the new symbol associated with the suffix of the sequence k
      //And I have to update the number of occurrences of each symbol
      //   if (toRead==0) {

	(*pWriter)((char*)&newSymb[vectTriple[k].seqN], 1);

	tableOcc_[currentPile].count_[whichPile[(int)newSymb[vectTriple[k].seqN]]]++;       //update the number of occurrences in BWT of the pileN[k]
	//std::cerr << "new number write " << numchar << "\n";
	cont++;    //number of read symbols
	//toRead--;
	//      }
      k++;   //  I changed the number of the sequence. New iteration.
    }
    
    //it means that posN[k]<>currentPile, so I have to change BWT-file
		//But before, I have to copy the remainder symbols from the old BWT to new BWT
    (*pReader).readAndSend((*pWriter));

    j=k;
  }

  delete pReader; pReader=NULL; // %%%
  delete pWriter; pWriter=NULL; // %%%
	
  //static FILE *tmpFile;
  //tmpFile = fopen("sizeFileBwt.txt", "a");
  
  //Renaming new to old
  for (dataTypedimAlpha g = 0 ; g < alphabetSize-1; g++) {  
    numchar=sprintf (filename, "bwt_%d", g);
    numchar=sprintf (filenameIn,"%s%s",filename,ext);
    numchar=sprintf (filenameOut,"new_%s%s",filename,ext);
    //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
    OutFileBWT = fopen(filenameOut, "r");
    // fileOut aperto?
    
    if (OutFileBWT!=NULL) { //If it exists
      fclose(OutFileBWT);
      if (remove(filenameIn)!=0) 
	std::cerr << filenameIn <<": Error deleting file" << std::endl;
      else
	if(rename(filenameOut,filenameIn))
	  std::cerr << filenameOut <<": Error renaming " << std::endl;
    }
    
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
  }
  std::cerr <<  std::endl;
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

void BCRexternalBWT::storeEntireBWT( const char* fn ) {

	static FILE *OutFileBWT, *InFileBWT;                  // output and input file BWT;
	char *filenameIn = new char[12];
	char *filename = new char[8];
	const char *ext = ".aux";
	dataTypeNChar numchar=0;
	dataTypeNChar numcharWrite=0;

	uchar *buffer = new uchar[SIZEBUFFER];

	dataTypeNChar *freqOut = new dataTypeNChar [256];
	for (unsigned i = 0; i < 255; ++i)
        freqOut[i] = 0;
	
	OutFileBWT = fopen(fn, "w");
			// fileOut aperto?
	if (OutFileBWT==NULL) {
				// NO, abort program
		std::cerr << "Entire BWT file: Error opening " << std::endl;
		exit (EXIT_FAILURE);
	}

	if (verboseEncode==1) {
		std::cerr << "\nThe last BWT-segment:"<< std::endl;
		uint mmm=0;
		while (mmm < alphabetSize) {			
			numchar=sprintf (filename, "bwt_%d", mmm);
			numchar=sprintf (filenameIn,"%s%s",filename,ext);
			//printf("===Current BWT-partial= %d\n",mmm);
			InFileBWT = fopen(filenameIn, "r");
			for (dataTypeNSeq g = 0 ; g < SIZEBUFFER; g++) 
				buffer[g] = '\0';
			numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,InFileBWT);
			std::cerr << "B[" << mmm << "]:\t";
			if (numchar==0)
				std::cerr  << "empty\n";
			else
				std::cerr  << buffer << "\n";
			fclose(InFileBWT);
			mmm++;
		}
	}


	std::cerr << "Entire BWT file" << std::endl;
	
	std::cerr << "Compute the distribution of chars \n";

	for (dataTypedimAlpha g = 0 ; g < alphabetSize; g++) {  
		numchar=sprintf (filename, "bwt_%d", g);
		numchar=sprintf (filenameIn,"%s%s",filename,ext);
		InFileBWT = fopen(filenameIn, "r");
			// fileOut aperto?
		if (InFileBWT==NULL) {
				// NO, abort program
			std::cerr << "BWT file " << (int)g <<": Error opening " << std::endl;
			exit (EXIT_FAILURE);
		}
		//std::cerr << "BWT file " << (int)g << "= ";
		while (numchar!=0) {
			numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,InFileBWT);
			//std::cerr << "number read " << numchar << "\n";			
			numcharWrite = fwrite (buffer, sizeof(uchar), numchar , OutFileBWT); 
			checkIfEqual(numchar,numcharWrite);  // we should always read/write the same number of characters

			for (unsigned j = 0 ; j < numchar; j++) 
				freqOut[(int)(buffer[j])]++;
		}
		
		fclose(InFileBWT);	


		//Delete or renome the partial BWT bwt_%d????
		if (deletePartialBWT == 1)          
		{
			if (remove(filenameIn)!=0) 
				std::cerr << "storeEntireBWT: Error deleting file" << std::endl;
		} 
		else //renome the aux bwt file
		{
			int lung = strlen(fn) + strlen(filenameIn)+1;
			char *newfilename = new char[lung];
			numchar=sprintf (newfilename,"%s%s",fn,filenameIn);
			//std::cerr  << newfilename << " " << filenameIn << std::endl;
			if(rename(filenameIn, newfilename))
					std::cerr  <<"storeEntireBWT: Error renaming file" << std::endl;
		}
	}
	fclose(OutFileBWT);

	if (verboseEncode==1) {
		std::cerr << "\nThe Entire BWT:"<< std::endl;
		OutFileBWT = fopen(fn, "r");		
		if (OutFileBWT==NULL) {
			// NO, abort program
			std::cerr << "Entire BWT file: Error opening " << std::endl;
			exit (EXIT_FAILURE);
		}

		for (dataTypeNSeq g = 0 ; g < SIZEBUFFER; g++) 
			buffer[g] = '\0';
		numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,OutFileBWT);
		if (numchar==0)
			std::cerr  << "empty\n";
		else
			std::cerr  << buffer << "\n";
		fclose(OutFileBWT);
	}
	delete [] buffer;


	
	freqOut[256]=0;
		std::cerr << "Distribution in BWT\n";
	for (dataTypedimAlpha i = 0; i < 255; ++i)
		if (freqOut[i] > 0) 
			std::cerr << i << " " << freqOut[i] << "\n";
	delete [] freqOut;
	
	delete [] filenameIn;
	delete [] filename;
}
