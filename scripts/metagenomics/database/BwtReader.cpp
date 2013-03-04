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
 	
 	
 	#include "BwtReader.hh"
 	#include "LetterCount.hh"
 	#include "BwtWriter.hh"
 	#include <cstring>
 	#include <cstdlib>
 	
 	//#define DEBUG 1
 	using namespace std;
 	
 	//
 	// BwtReaderBase member function definitions
 	//
 	
 	int BwtReaderBase::readAndCount( LetterCount& c )
 	{ // will call the relevant virtual function
 	return readAndCount( c, maxLetterCountType );
 	} // ~readAndCount
 	
 	
 	//
 	// BwtReaderASCII member function definitions
 	//
 	
 	void BwtReaderASCII::rewindFile(void)
 	{
 	rewind(pFile_);
 	currentPos_=0;
 	} // ~rewindFile
 	
 	LetterCountType BwtReaderASCII::tellg( void ) const
 	{
 	return currentPos_;
 	} // ~tellg
 	
 	
 	
 	int BwtReaderASCII::readAndCount( LetterCount& c, const LetterCountType numChars )
 	{
 	#ifdef DEBUG
 	std::cout << "readAndCount " << numChars << " chars " << endl;
 	#endif
 	LetterCountType charsLeft(numChars), charsToRead, charsRead;
 	while (charsLeft>0)
 	{
 	charsToRead=((charsLeft>ReadBufferSize)?ReadBufferSize:charsLeft);
 	charsRead=fread( buf_, sizeof(char), charsToRead, pFile_ );
 	#ifdef DEBUG
 	std::cout << "Reading " << charsRead << " chars ";
 	#endif
 	for (LetterCountType i(0);i<charsRead;i++)
 	{
 	#ifdef DEBUG
 	std::cout << buf_[i];
 	#endif
 	c.count_[whichPile[(int)buf_[i]]]++;
 	}
 	#ifdef DEBUG
 	std::cout << std::endl;
 	#endif
 	charsLeft-=charsRead;
 	if (charsRead<charsToRead)
 	{ // did not get everything asked for! return num of chars actually found
 	currentPos_+=(numChars-charsLeft);
 	return (numChars-charsLeft);
 	} // ~if
 	} // ~while
 	currentPos_+=numChars;
 	return numChars;
 	} // ~int BwtReaderASCII::readAndCount( LetterCount& c, const int numChars )
 	
 	
 	//
 	// BwtReaderRunLength member function definitions
 	//
 	
 	BwtReaderRunLength::BwtReaderRunLength( const std::string& fileName ):
 	BwtReaderBase(fileName), runLength_(0),
 	pBuf_(buf_+ReadBufferSize),pBufMax_(buf_+ReadBufferSize),
 	lastChar_(notInAlphabet),
 	finished_(false)
 	{
 	uint j;
 	for (uint i(0);i<256;i++)
 	{
 	lengths_[i]=1+(i>>4);
 	j=(i&0xF);
 	codes_[i]=(j<alphabetSize)?alphabet[j]:notInAlphabet;
 	} // ~for i
 	} // ~ctor
 	
 	void BwtReaderRunLength::rewindFile(void)
 	{ // rewind file and set all vars as per constructor
 	rewind(pFile_);
 	runLength_=0;
 	pBuf_=buf_+ReadBufferSize;
 	pBufMax_=buf_+ReadBufferSize;
 	lastChar_=notInAlphabet;
 	currentPos_=0;
 	finished_=false;
 	} // ~rewindFile
 	
 	LetterCountType BwtReaderRunLength::tellg( void ) const
 	{
 	return currentPos_;
 	} // ~tellg
 	
 	
 	int BwtReaderRunLength::readAndCount( LetterCount& c, const LetterCountType numChars )
 	{
 	#ifdef DEBUG
 	std::cout << "readAndCount " << numChars << " chars " << endl;
 	#endif
 	LetterCountType charsLeft(numChars);
 	while(charsLeft>runLength_)
 	{
 	// Below is not great design, at first call of this function it accesses an
 	// out-of-range array element. Fortunately it always adds zero to it! :)
 	c.count_[whichPile[lastChar_]]+=runLength_;
 	charsLeft-=runLength_;
 	if (getRun()==false)
 	{
 	currentPos_+=(numChars-charsLeft);
 	return (numChars-charsLeft);
 	// assert(1==0);
 	} // ~if
 	} // ~while
 	
 	c.count_[whichPile[lastChar_]]+=charsLeft;
 	runLength_-=charsLeft;
 	currentPos_+=numChars;
 	return numChars;
 	} // ~BwtReaderRunLength::readAndCount( LetterCount& c, const int numChars )
 	
 	int BwtReaderRunLength::readAndSend(BwtWriterBase& writer, const int numChars )
 	{
 	int charsLeft(numChars);
 	while(charsLeft>runLength_)
 	{
 	// int fred(whichPile[lastChar_]);
 	writer.sendRun( lastChar_, runLength_ );
 	// c.count_[whichPile[lastChar_]]+=runLength_;
 	charsLeft-=runLength_;
 	if (getRun()==false)
 	{
 	currentPos_+=(numChars-charsLeft);
 	return (numChars-charsLeft);
 	// assert(1==0);
 	} // ~if
 	} // ~while
 	
 	writer.sendRun(lastChar_, charsLeft);
 	// c.count_[whichPile[lastChar_]]+=charsLeft;
 	runLength_-=charsLeft;
 	currentPos_+=numChars;
 	return numChars;
 	} //~BwtReaderRunLength::readAndSend(BwtWriterBase& writer, const int numChars)
 	
 	int BwtReaderRunLength::operator()( char* p, int numChars )
 	{
 	#ifdef DEBUG
 	std::cout << "B asked for " << numChars << " " << lastChar_ << " "
 	<< runLength_ << " " << pFile_ << std::endl;
 	#endif
 	int charsLeft(numChars);
 	// return fread( p, sizeof(char), numChars, pFile_ );
 	while(charsLeft>runLength_)
 	{
 	
 	memset(p, lastChar_, runLength_);
 	p+=runLength_;
 	
 	charsLeft-=runLength_;
 	if (getRun()==false)
 	{
 	// runLength_=0;
 	#ifdef DEBUG
 	std::cout << "B read " << numChars-charsLeft << " out of "
 	<< numChars << std::endl;
 	#endif
 	currentPos_+=(numChars-charsLeft);
 	return (numChars-charsLeft);
 	} // ~if
 	} // ~while
 	runLength_=lengths_[lastChar_];
 	lastChar_=codes_[lastChar_];
 	
 	memset(p, lastChar_, charsLeft);
 	
 	runLength_-=charsLeft;
 	#ifdef DEBUG
 	std::cout << "B delivered " << numChars << " " << charsLeft << " "
 	<< pFile_ << std::endl;
 	#endif
 	currentPos_+=numChars;
 	return numChars;
 	} // ~operator()
 	
 	bool BwtReaderRunLength::getRun(void)
 	{
 	if (pBuf_==pBufMax_)
 	{
 	if (finished_) return false;
 	else
 	{
 	int numRead(fread( buf_, sizeof(uchar),
 	ReadBufferSize, pFile_ ));
 	if (numRead==0) return false;
 	else if (numRead<ReadBufferSize)
 	{
 	finished_=true;
 	pBufMax_=buf_+numRead;
 	}
 	pBuf_=buf_;
 	} // ~else
 	} // ~if
 	runLength_=lengths_[*pBuf_];
 	lastChar_=codes_[*pBuf_];
 	#ifdef DEBUG
 	cout << "Got run: " << runLength_ << " of " << lastChar_ << endl;
 	#endif
 	pBuf_++;
 	
 	return true;
 	
 	} // ~getRun
