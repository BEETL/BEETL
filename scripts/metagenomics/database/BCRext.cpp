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
 	
 	
 	#include <iostream>
 	#include <fstream>
 	#include <string>
 	#include <cassert>
 	#include <cstring>
 	#include <vector>
 	#include "Types.hh"
 	#include "Config.hh"
 	#include "Timer.hh"
 	
 	#include "LetterCount.hh"
 	#include "ReadBuffer.hh"
 	#include "BwtWriter.hh"
 	#include "BwtReader.hh"
 	
 	#include "BCRext.hh"
 	
 	
 	
 	#ifdef USE_POSIX_FILE_OPTIMIZATIONS
 	#define _XOPEN_SOURCE 600
 	#endif
 	#include <fcntl.h>
 	
 	
 	static const char* id="@(#) $Id: BCRext.cpp,v 1.13 2011/10/13 09:16:02 acox Exp $";
 	
 	
 	using namespace std;
 	
 	//const int maxSeqSize(1023);
 	//const int bwtBufferSize(16384); // 1<<20=1048576
 	//typedef unsigned int uint;
 	
 	// below limits to 4 billion reads max - change to unsigned long for more
 	//typedef unsigned int SequenceNumberType;
 	
 	// Should work for BWT of up to 2^64 characters in size
 	//typedef unsigned long LetterCountType;
 	
 	
 	
 	
 	//typedef ReadBufferASCII ReadBuffer;
 	#ifdef USE_4_BITS_PER_BASE
 	typedef ReadBuffer4Bits ReadBuffer;
 	#else
 	#ifdef USE_PREFIX_ONLY
 	typedef ReadBuffer4Bits ReadBuffer;
 	#else
 	typedef ReadBufferASCII ReadBuffer;
 	#endif
 	#endif
 	
 	
 	#ifdef COMPRESS_BWT
 	typedef BwtReaderRunLength BwtReader;
 	typedef BwtWriterRunLength BwtWriter;
 	#else
 	typedef BwtReaderASCII BwtReader;
 	typedef BwtWriterASCII BwtWriter;
 	#endif
 	
 	
 	const string tmp1("temp1-XX");
 	const string tmp2("temp2-XX");
 	//const int tmp1Size(strlen(tmp1));
 	//const int tmp2Size(strlen(tmp2));
 	
 	const string fileStem("bwt");
 	const string fileStemTemp("bwt_tmp");
 	
 	
 	void getFileName( const string& stem, const char code, const int pile,
 	string& fileName )
 	{
 	fileName=stem;
 	fileName+='-';
 	fileName+=code;
 	fileName+='0';
 	assert(pile<=9);
 	fileName+=(char)(48+pile);
 	// cerr << "Made file name " << fileName << endl;
 	}
 	
 	
 	
 	int main( int numArgs, char* args[])
 	{
 	Timer timer;
 	
 	string prefix = (string)"[" + (string)args[0] + (string)"]: ";
 	cerr << prefix << "time now is " << timer.timeNow();
 	cerr << prefix << "command line is";
 	for (int i(0);i<numArgs;i++) cerr << " " << args[i];
 	cerr << endl;
 	cerr << prefix << "software version is " << id << endl;
 	
 	
 	
 	string tmpIn=fileStem;
 	string tmpOut=fileStemTemp;
 	string tmpSwap;
 	string fileName;
 	
 	// output streams for sequences - 1 per pile
 	vector <FILE*> outSeq(alphabetSize);
 	// output streams for positions of suffixes in pile - 1 per pile
 	vector <FILE*> outPtr(alphabetSize);
 	// output stream for updated BWTs - 1 per pile
 	vector <BwtWriterBase*> outBwt(alphabetSize);
 	// output streams for original read numberings - 1 per pile
 	vector <FILE*> outNum(alphabetSize);
 	
 	// input streams for BWT from previous iter - 1 per pile
 	// this one used to compute the counts to do backward search
 	vector <BwtReader*> inBwt(alphabetSize);
 	
 	// input streams for BWT from previous iter - 1 per pile
 	// this one used to read the BWT chunk by chunk to allow new
 	// chars to be interspersed
 	vector <BwtReader*> inBwt2(alphabetSize);
 	
 	FILE* inSeq;
 	// FILE* inPtr;
 	// FILE* inNum;
 	FILE* outDollarBwt;
 	
 	// string thisSeq;
 	char seqBuf[maxSeqSize+1];
 	// char* seqBuff;
 	
 	vector<char> bwtBuf;
 	// extra byte accounts for a fact that inserted BWT characters
 	// are appended to the buffer after the interspersed characters so as to
 	// make a single write to file
 	bwtBuf.resize(bwtBufferSize+1);
 	
 	bool compressedOutput(false);
 	int nextArg(1);
 	
 	if (numArgs!=2)
 	{
 	if ((numArgs==3)&&(strcmp(args[1],"-c")==0))
 	{
 	cerr << prefix << "Will output compressed BWT strings" << endl;
 	nextArg++;
 	compressedOutput=true;
 	} // ~if
 	else
 	{
 	cerr << "Usage #1: " << args[0] << " [-c] inputFile" << endl
 	<< "Usage #2: cat inputFile | " << args[0] << " [-c] -" << endl
 	<< "-c: output BWT files in compressed form" << endl
 	<< "inputFile: input sequences - ASCII, carriage return after each"
 	<< endl;
 	assert(1==0);
 	} // ~else
 	} // ~if
 	
 	assert(whichPile[(int)'$']==0);
 	assert(whichPile[(int)'A']==1);
 	assert(whichPile[(int)'C']==2);
 	assert(whichPile[(int)'G']==3);
 	// next two assertions handle both possible orderings of N and T
 	assert(whichPile[(int)alphabet[4]]==4);
 	assert(whichPile[(int)alphabet[5]]==5);
 	// assert(whichPile['T']==3);
 	// assert(whichPile['N']==4);
 	assert(whichPile[(int)notInAlphabet]==-1);
 	cerr << prefix << "Using alphabet = " << alphabet << ", size = " << alphabetSize
 	<< endl;
 	// cerr << BUFSIZ << endl;
 	assert(alphabetSize<10); // otherwise need to mod file names etc
 	assert(sizeof(LetterCountType)==8);
 	
 	
 	SequenceNumberType seqNum(0);
 	LetterCountType seqPtr;
 	int thisPile, lastPile;
 	LetterCountType posInPile;
 	// char inChar;
 	int charsToGrab;// charsLeft, charsThisBatch; // TBD make these long?
 	
 	
 	if (strcmp(args[nextArg],"-")==0)
 	{
 	cerr << prefix << "Will read sequences from standard input" << endl;
 	inSeq=stdin;
 	} // ~if
 	else
 	{
 	cerr << prefix
 	<< "Will read sequences from file " << args[nextArg] << endl;
 	inSeq=fopen(args[nextArg],"r");
 	assert(inSeq!=NULL);
 	} // ~else
 	
 	// read first sequence to determine read size
 	assert( fgets( seqBuf, maxSeqSize, inSeq)!=NULL);
 	
 	const int seqSize(strlen(seqBuf)-1); // -1 compensates for \n at end
 	cerr << prefix << "Assuming all sequences are of length " << seqSize << endl;
 	// inFile.seekg(0,ios::beg);
 	// rewind(inSeq);
 	
 	if ((seqSize%2)==1) // if odd
 	{
 	// cout << "ODD" << endl;
 	tmpIn=fileStem;
 	tmpOut=fileStemTemp;
 	} // ~if
 	else
 	{
 	// cout << "EVEN" << endl;
 	tmpIn=fileStemTemp;
 	tmpOut=fileStem;
 	} // ~else
 	
 	getFileName(fileStem,'B',0,fileName);
 	outDollarBwt=fopen(fileName.c_str(),"w");
 	assert(inSeq!=NULL);
 	
 	
 	for (int j(1);j<alphabetSize;j++)
 	{
 	
 	getFileName(tmpIn,'S',j,fileName);
 	outSeq[j]= fopen(fileName.c_str(), "w");
 	assert(outSeq[j]!=NULL);
 	
 	getFileName(tmpIn,'P',j,fileName);
 	outPtr[j]= fopen(fileName.c_str(), "w");
 	assert(outPtr[j]!=NULL);
 	
 	#ifdef TRACK_SEQUENCE_NUMBER
 	getFileName(tmpIn,'N',j,fileName);
 	outNum[j]= fopen(fileName.c_str(), "w");
 	assert(outNum[j]!=NULL);
 	#endif
 	
 	getFileName(tmpIn,'B',j,fileName);
 	outBwt[j]=new BwtWriter(fileName.c_str());
 	} // ~for
 	
 	LetterCount dollars;
 	// vector<LetterCount> alreadyInPile(alphabetSize);
 	// vector<LetterCount> countedThisIter(alphabetSize);
 	LetterCountEachPile alreadyInPile;
 	LetterCountEachPile countedThisIter;
 	LetterCountEachPile newCharsThisIter;
 	
 	// TBD Rationalize count names wrt first and subsequent iterations
 	LetterCount addedSoFar;
 	LetterCount outputSoFar;
 	
 	LetterCount prevCharsOutputThisIter;
 	LetterCount newCharsAddedThisIter;
 	
 	
 	// LetterCount readSoFar[alphabetSize];
 	
 	// First iteration
 	// - move original sequence into piles based on last character
 	// - work out BWT corresponding to 0-suffixes and 1-suffixes
 	// TBD check for invalid chars, do qual masking
 	
 	ReadBuffer readBuffer(seqSize,-1,-1,-1);
 	assert(readBuffer.blockSize_>seqSize+1);
 	
 	// copy first sequence over
 	strcpy( readBuffer.seqBufBase_, seqBuf);
 	// while ( fgets( readBuffer.seqBufBase_, readBuffer.blockSize_, inSeq)!=NULL)
 	do
 	{
 	thisPile=whichPile[(int)readBuffer.seqBufBase_[seqSize-1]];
 	assert(thisPile>0); // zero is terminator so should not be present
 	assert(thisPile<alphabetSize);
 	#ifdef DEBUG
 	cout << readBuffer.seqBufBase_ << endl;
 	#endif
 	
 	// count characters and output first N chars of BWT
 	// (those preceding terminator characters)
 	dollars.count_[thisPile]++;
 	assert(fwrite( readBuffer.seqBufBase_+seqSize-1, sizeof(char), 1, outDollarBwt )==1);
 	
 	// fprintf( outSeq[thisPile], "%s", readBuffer.seqBufBase_);
 	readBuffer.convertFromASCII();
 	readBuffer.sendTo( outSeq[thisPile] );
 	
 	#ifdef TRACK_SEQUENCE_NUMBER
 	assert(fwrite( &seqNum, sizeof(SequenceNumberType),
 	1, outNum[thisPile] )==1);
 	// seqNum++;
 	#endif
 	
 	// create BWT corresponding to 1-suffixes
 	assert(whichPile[(int)readBuffer.seqBufBase_[seqSize-2]]>=0);
 	assert(whichPile[(int)readBuffer.seqBufBase_[seqSize-2]]<alphabetSize);
 	countedThisIter[thisPile].count_[whichPile[(int)readBuffer.seqBufBase_[seqSize-2]]]++;
 	// assert(fwrite( readBuffer.seqBufBase_+seqSize-2, sizeof(char), 1, outBwt[thisPile] )==1);
 	(*outBwt[thisPile])( readBuffer.seqBufBase_+seqSize-2, 1 );
 	
 	assert(fwrite( addedSoFar.count_+thisPile, sizeof(LetterCountType),
 	1, outPtr[thisPile] )==1);
 	addedSoFar.count_[thisPile]++;
 	seqNum++;
 	} // ~while
 	while ( fgets( readBuffer.seqBufBase_, readBuffer.blockSize_, inSeq)!=NULL);
 	
 	fclose (inSeq);
 	fclose (outDollarBwt);
 	
 	cerr << prefix << "Read " << seqNum << " sequences" << endl;
 	for (int i(1);i<alphabetSize;i++)
 	{
 	fclose(outSeq[i]);
 	fclose(outPtr[i]);
 	#ifdef TRACK_SEQUENCE_NUMBER
 	fclose(outNum[i]);
 	#endif
 	delete outBwt[i];
 	// fclose(outBwt[i]);
 	}
 	// return (0);
 	
 	// ReadBuffer buffer(seqSize);
 	
 	// Main loop
 	for (int i(2);i<=seqSize;i++)
 	{
 	
 	
 	cout << "Starting iteration " << i << ", time now: " << timer.timeNow();
 	cout << "Starting iteration " << i << ", usage: " << timer << endl;
 	
 	// don't do j=0 - this is the $ sign which is done already
 	for (int j(1);j<alphabetSize;j++)
 	{
 	// prep the output files
 	
 	getFileName(tmpOut,'S',j,fileName);
 	outSeq[j]= fopen(fileName.c_str(), "w");
 	assert(outSeq[j]!=NULL);
 	
 	getFileName(tmpOut,'P',j,fileName);
 	outPtr[j]= fopen(fileName.c_str(), "w");
 	assert(outPtr[j]!=NULL);
 	
 	#ifdef TRACK_SEQUENCE_NUMBER
 	getFileName(tmpOut,'N',j,fileName);
 	outNum[j]= fopen(fileName.c_str(), "w");
 	assert(outNum[j]!=NULL);
 	#endif
 	
 	getFileName(tmpOut,'B',j,fileName);
 	
 	if (i==seqSize)
 	{
 	if (compressedOutput==true)
 	outBwt[j] = new BwtWriterRunLength( fileName.c_str() );
 	else
 	outBwt[j] = new BwtWriterASCII( fileName.c_str() );
 	}
 	else
 	outBwt[j] = new BwtWriter( fileName.c_str() );
 	
 	#ifdef DEBUG
 	cout << "Prepping output file " << tmpOut << endl;
 	#endif
 	
 	setvbuf( outSeq[j], NULL, _IOFBF, 262144);
 	// setvbuf( outPtr[j], NULL, _IOFBF, 65536);
 	// setvbuf( outNum[j], NULL, _IOFBF, 65536);
 	// setvbuf( outBwt[j], NULL, _IOFBF, 65536);
 	
 	// prep the input files
 	getFileName(tmpIn,'B',j,fileName);
 	inBwt[j]=new BwtReader( fileName.c_str() );
 	inBwt2[j]=new BwtReader( fileName.c_str() );
 	
 	#ifdef DEBUG
 	cout << "Prepping input file " << tmpIn << endl;
 	#endif
 	
 	} // ~for j
 	
 	addedSoFar.clear();
 	outputSoFar.clear();
 	
 	prevCharsOutputThisIter.clear();
 	newCharsAddedThisIter.clear();
 	
 	#ifdef DEBUG
 	cout << "already in pile" << endl;
 	alreadyInPile.print();
 	cout << "counted this iter" << endl;
 	countedThisIter.print();
 	#endif
 	countedThisIter.clear();
 	newCharsThisIter.clear();
 	
 	#ifdef DEBUG
 	cout << "Count in dollars pile: ";
 	dollars.print();
 	#endif
 	
 	int fdSeq, fdNum, fdPtr;
 	
 	// don't do j=0; $ sign done already
 	for (int j(1);j<alphabetSize;j++)
 	{ // read each input file in turn
 	
 	getFileName(tmpIn,'S',j,fileName);
 	fdSeq=open(fileName.c_str(),O_RDONLY,0);
 	assert(fdSeq!=-1);
 	
 	getFileName(tmpIn,'P',j,fileName);
 	fdPtr=open(fileName.c_str(),O_RDONLY,0);
 	assert(fdPtr!=-1);
 	
 	#ifdef TRACK_SEQUENCE_NUMBER
 	getFileName(tmpIn,'N',j,fileName);
 	fdNum=open(fileName.c_str(),O_RDONLY,0);
 	assert(fdNum!=-1);
 	#ifdef USE_POSIX_FILE_OPTIMIZATIONS
 	assert(posix_fadvise(fdNum,0,0,POSIX_FADV_SEQUENTIAL|POSIX_FADV_NOREUSE|POSIX_FADV_WILLNEED)!=-1);
 	#endif
 	#endif
 	
 	#ifdef USE_POSIX_FILE_OPTIMIZATIONS
 	assert(posix_fadvise(fdSeq,0,0,POSIX_FADV_SEQUENTIAL|POSIX_FADV_NOREUSE|POSIX_FADV_WILLNEED)!=-1);
 	assert(posix_fadvise(fdPtr,0,0,POSIX_FADV_SEQUENTIAL|POSIX_FADV_NOREUSE|POSIX_FADV_WILLNEED)!=-1);
 	#endif
 	
 	#ifdef USE_PREFIX_ONLY
 	ReadBufferPrefix buffer(seqSize, i, fdSeq, fdNum, fdPtr);
 	#else
 	ReadBuffer buffer(seqSize, fdSeq, fdNum, fdPtr);
 	#endif
 	
 	while ( buffer.getNext( seqNum, seqPtr))
 	{
 	
 	thisPile=buffer[seqSize-i];
 	// thisPile=whichPile[seqBuff[seqSize-i]];
 	assert(thisPile>=0);
 	lastPile=buffer[seqSize-i+1];
 	// lastPile=whichPile[seqBuff[seqSize-i+1]];
 	assert(lastPile>=0);
 	
 	
 	
 	
 	#ifdef DEBUG
 	cout << "Read in " << seqPtr << " " << seqNum << " " << thisPile << " " << lastPile << endl;
 	for (int ZZ(0);ZZ<seqSize;ZZ++)
 	{ cout << "ZZ " <<ZZ <<endl; cout << alphabet[buffer[ZZ]] <<endl;}
 	cout << endl;
 	
 	#endif
 	
 	// *** work out position in new pile ***
 	
 	// sum contents of lexicographically smaller piles
 	// ... probably possible to speed this up by storing cumulative sums
 	#ifdef DEBUG
 	cout << "already in pile" << endl;
 	alreadyInPile.print();
 	#endif
 	
 	// posInPile=0;
 	posInPile=dollars.count_[thisPile];
 	// cout << posInPile << " " << thisPile << " " << lastPile << endl;
 	for (int k(1);k<lastPile; k++)
 	{
 	posInPile+=alreadyInPile[k].count_[thisPile];
 	// cout << posInPile << endl;
 	} // ~for k
 	
 	#ifdef DEBUG
 	cout << "posInPile starts at " << posInPile << endl;
 	cout << "counting in pile " << alphabet[lastPile] << endl;
 	#endif
 	
 	// count all chars prior to seqPtr in lastPile
 	// read seqPtr too, but don't add to countedThisIter
 	charsToGrab=seqPtr-addedSoFar.count_[lastPile];//+1;
 	#ifdef DEBUG
 	cout << "charsToGrab " << charsToGrab << endl;
 	#endif
 	
 	// Should now always read at least 1 byte
 	
 	assert(charsToGrab>=0);
 	
 	assert(inBwt[lastPile]->readAndCount
 	(countedThisIter[lastPile],charsToGrab)==charsToGrab);
 	inBwt[lastPile]->readAndCount(newCharsThisIter[lastPile],1);
 	
 	
 	addedSoFar.count_[lastPile]=seqPtr+1;
 	
 	
 	posInPile+=countedThisIter[lastPile].count_[thisPile];
 	
 	
 	
 	#ifdef DEBUG
 	cout << "counted this iter" << endl;
 	countedThisIter.print();
 	#endif
 	
 	// *** add char into new pile ***
 	
 	// read and output bytes up to insertion point
 	
 	charsToGrab=posInPile-prevCharsOutputThisIter.count_[thisPile];
 	
 	assert(inBwt2[thisPile]->readAndSend( *outBwt[thisPile], charsToGrab )
 	==charsToGrab);
 	// bwtBuf[0]=(seqSize-i-1>=0)?baseNames[buffer[seqSize-i-1]]:'$';
 	bwtBuf[0]=(seqSize-i-1>=0)?alphabet[buffer[seqSize-i-1]]:alphabet[0];
 	(*outBwt[thisPile])( &bwtBuf[0], 1 );
 	
 	prevCharsOutputThisIter.count_[thisPile]=posInPile;
 	
 	
 	// pointer into new pile must be offset by number of new entries added
 	// seqPtr+=newCharsAddedThisIter.count_[thisPile];
 	seqPtr=posInPile+newCharsAddedThisIter.count_[thisPile];
 	assert(fwrite( &seqPtr, sizeof(LetterCountType),
 	1, outPtr[thisPile] )==1);
 	#ifdef DEBUG
 	cout << "adding pointer " << seqPtr << " to pile "
 	<< alphabet[thisPile] << endl;
 	#endif
 	// then the offset itself is updated
 	newCharsAddedThisIter.count_[thisPile]++;
 	
 	// do radix sort
 	// fprintf( outSeq[thisPile], "%s\n", seqBuff);
 	// assert(fwrite( seqBuff, sizeof(char),
 	// 1+seqSize, outSeq[thisPile] )==1+seqSize);
 	buffer.sendTo( outSeq[thisPile] );
 	
 	
 	#ifdef TRACK_SEQUENCE_NUMBER
 	assert(fwrite( &seqNum, sizeof(SequenceNumberType),
 	1, outNum[thisPile] )==1);
 	#endif
 	
 	} // ~while
 	
 	
 	// fclose(inSeq);
 	close(fdSeq);
 	close(fdPtr);
 	#ifdef TRACK_SEQUENCE_NUMBER
 	close(fdNum);
 	#endif
 	
 	} // ~for j
 	
 	cout << "All new characters inserted, usage: " << timer << endl;
 	for (int j(1);j<alphabetSize;j++)
 	{
 	
 	while(inBwt[j]->readAndCount
 	(countedThisIter[j],ReadBufferSize)==ReadBufferSize);
 	
 	} // ~for j
 	
 	cout << "finishing off BWT strings" << endl;
 	for (int j(1);j<alphabetSize;j++)
 	{
 	
 	while(inBwt2[j]->readAndSend( *outBwt[j], ReadBufferSize )
 	==ReadBufferSize);
 	
 	} // ~for
 	
 	#ifdef DEBUG
 	cout << "final value of counted this iter" << endl;
 	countedThisIter.print();
 	cout << "final value of new chars this iter" << endl;
 	newCharsThisIter.print();
 	#endif
 	
 	alreadyInPile+=newCharsThisIter;
 	
 	for (int j(1);j<alphabetSize;j++)
 	{
 	fclose(outSeq[j]);
 	fclose(outPtr[j]);
 	#ifdef TRACK_SEQUENCE_NUMBER
 	fclose(outNum[j]);
 	#endif
 	delete (outBwt[j]);
 	delete (inBwt[j]);
 	delete (inBwt2[j]);
 	} // ~for j
 	
 	tmpSwap=tmpIn;
 	tmpIn=tmpOut;
 	tmpOut=tmpSwap;
 	
 	cout << "finished iteration " << i << ", usage: " << timer << endl;
 	
 	// assert(i<2);
 	} // ~for i (main iteration)
 	
 	#ifdef REMOVE_TEMPORARY_FILES
 	string fileTypes("SPB");
 	for (int j(1);j<alphabetSize;j++)
 	{
 	for (int i(0);i<fileTypes.size();i++)
 	{
 	getFileName(tmpOut,fileTypes[i],j,fileName);
 	if (remove(fileName.c_str())!=0)
 	{
 	cerr << "Warning: failed to clean up temporary file " << fileName
 	<< endl;
 	} // ~if
 	else
 	{
 	cerr << "Removed temporary file " << fileName
 	<< endl;
 	} // ~if
 	} // ~for i
 	} // ~for j
 	#endif
 	
 	cerr << "Final output files are named "<< tmpIn << " and similar" << endl;
 	return (0);
 	}
