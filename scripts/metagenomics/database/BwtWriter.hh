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
 	
 	#ifndef INCLUDED_BWTWRITER_HH
 	#define INCLUDED_BWTWRITER_HH
 	
 	#include <string>
 	#include <cstdio>
 	
 	
 	
 	struct BwtWriterBase
 	{
 	BwtWriterBase( const std::string& fileName ) : pFile_(fopen(fileName.c_str(),"w"))
 	{
 	assert(pFile_!=NULL);
 	// setvbuf( pFile_, NULL, _IOFBF, 262144);
 	}
 	virtual ~BwtWriterBase() { fclose(pFile_); }
 	
 	virtual void operator()( const char* p, int numChars ) =0;
 	
 	
 	virtual void sendRun( char c, int runLength ) =0;
 	
 	
 	FILE* pFile_;
 	}; // ~BwtWriterBase
 	
 	struct BwtWriterASCII : public BwtWriterBase
 	{
 	BwtWriterASCII( const std::string& fileName ) : BwtWriterBase(fileName)
 	{
 	}
 	virtual ~BwtWriterASCII() {}
 	
 	virtual void operator()( const char* p, int numChars )
 	{
 	assert(fwrite( p, sizeof(char), numChars, pFile_ )==(size_t)numChars);
 	} // ~operator()
 	
 	virtual void sendRun( char c, int runLength )
 	{
 	for (int i(0);i<runLength;i++) fprintf( pFile_, "%c", c );
 	}
 	
 	
 	// FILE* pFile_;
 	};
 	
 	
 	struct BwtWriterRunLength : public BwtWriterBase
 	{
 	BwtWriterRunLength( const std::string& fileName ):
 	BwtWriterBase(fileName),
 	runLength_(0),pBuf_(buf_),pBufMax_(buf_+ReadBufferSize), lastChar_('Z')
 	#ifdef REPORT_COMPRESSION_RATIO
 	, charsReceived_(0), bytesWritten_(0)
 	#endif
 	{}
 	
 	virtual ~BwtWriterRunLength()
 	{
 	if (runLength_!=0) encodeRun( lastChar_, runLength_);
 	
 	
 	if (pBuf_!=buf_)
 	{
 	assert(fwrite( buf_, sizeof(char), (pBuf_-buf_),
 	pFile_ )==(size_t)(pBuf_-buf_));
 	#ifdef REPORT_COMPRESSION_RATIO
 	bytesWritten_+=(LetterCountType)(pBuf_-buf_);
 	#endif
 	}
 	
 	#ifdef REPORT_COMPRESSION_RATIO
 	std::cout << "BwtWriterRunLength: received "
 	<< charsReceived_ << " chars, sent "
 	<< bytesWritten_ << " bytes, compression "
 	<< ((double)8*bytesWritten_)/(charsReceived_)
 	<< " bits per char " << std::endl;
 	#endif
 	
 	
 	}
 	
 	virtual void operator()( const char* p, int numChars )
 	{
 	
 	for (int i(0);i<numChars;i++)
 	{
 	if ((*p)==lastChar_) { runLength_++; }
 	else
 	{
 	if (runLength_>0)
 	{
 	encodeRun( lastChar_, runLength_ );
 	} // ~if
 	runLength_=1; lastChar_=*p;
 	} // ~else
 	p++;
 	} // ~for
 	// assert(fwrite( p, sizeof(char), numChars, pFile_ )==numChars);
 	} // ~operator()
 	
 	void sendChar( char c)
 	{
 	*pBuf_=c;
 	if (++pBuf_==pBufMax_)
 	{
 	#ifdef REPORT_COMPRESSION_RATIO
 	bytesWritten_+=ReadBufferSize;
 	#endif
 	assert(fwrite( buf_, sizeof(char), ReadBufferSize, pFile_ )
 	==(size_t)ReadBufferSize);
 	pBuf_=buf_;
 	}
 	}
 	
 	void encodeRun( char c, uint runLength )
 	{
 	#ifdef DEBUG
 	std::cout << "B sending run " << c << " " << runLength << " " << pFile_
 	<< std::endl;
 	#endif
 	#ifdef REPORT_COMPRESSION_RATIO
 	charsReceived_+=runLength;
 	#endif
 	int charIndex(whichPile[(int)c]);
 	assert(charIndex!=nv);
 	uchar outCode(0xF0|((uchar)charIndex));
 	runLength--;
 	const uint numMaxChars(runLength>>4);
 	for (uint i(0);i<numMaxChars;i++)
 	{
 	sendChar(outCode);
 	}
 	runLength&=(uint)0xF;
 	
 	outCode=(((uchar)runLength)<<4);
 	outCode|=charIndex;
 	// assert(((uint)outCode)<256);
 	// assert(fwrite( &outCode, sizeof(char), 1, pFile_ )==1);
 	sendChar(outCode);
 	#ifdef DEBUG
 	std::cout << "B sending " << (uint)outCode << " " << pFile_ << std::endl;
 	#endif
 	} // ~encodeRun
 	
 	virtual void sendRun( char c, int runLength )
 	{
 	if (c==lastChar_) { runLength_+=runLength; }
 	else
 	{
 	if (runLength_!=0) encodeRun(lastChar_,runLength_);
 	lastChar_=c;
 	runLength_=runLength;
 	}
 	} // ~sendRun
 	
 	uint runLength_;
 	uchar buf_[ReadBufferSize];
 	uchar* pBuf_;
 	const uchar* pBufMax_;
 	uchar lastChar_;
 	#ifdef REPORT_COMPRESSION_RATIO
 	LetterCountType charsReceived_;
 	LetterCountType bytesWritten_;
 	#endif
 	}; // ~BwtWriterRunLength
 	
 	
 	
 	#endif
