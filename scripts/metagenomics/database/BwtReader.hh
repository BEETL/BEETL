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
 	
 	
 	
 	#ifndef INCLUDED_BWTREADER_HH
 	#define INCLUDED_BWTREADER_HH
 	
 	#include <string>
 	#include <cassert>
 	#include "Config.hh"
 	#include "Types.hh"
 	//#include "LetterCount.hh"
 	
 	class BwtWriterBase;
 	class LetterCount;
 	
 	
 	
 	class BwtReaderBase
 	{
 	public:
 	BwtReaderBase( const std::string& fileName ) :
 	pFile_(fopen(fileName.c_str(),"r"))
 	{
 	assert(pFile_!=NULL);
 	}
 	virtual ~BwtReaderBase() { fclose(pFile_); }
 	
 	virtual int readAndCount( LetterCount& c, const LetterCountType numChars )=0;
 	int readAndCount( LetterCount& c );
 	virtual int readAndSend( BwtWriterBase& writer, const int numChars )=0;
 	virtual int operator()( char* p, int numChars )=0;
 	virtual void rewindFile(void)=0;
 	virtual LetterCountType tellg( void ) const=0;
 	
 	protected:
 	FILE* pFile_;
 	
 	char buf_[ReadBufferSize];
 	}; // ~class BwtReaderBase
 	
 	class BwtReaderASCII : public BwtReaderBase
 	{
 	public:
 	BwtReaderASCII( const std::string& fileName ) :
 	BwtReaderBase(fileName),
 	currentPos_(0)
 	{
 	}
 	
 	virtual ~BwtReaderASCII() {}
 	
 	virtual int readAndCount( LetterCount& c, const LetterCountType numChars );
 	
 	virtual int readAndSend( BwtWriterBase& writer, const int numChars )
 	{ assert(1==0); // %%% writeme!
 	}
 	
 	virtual int operator()( char* p, int numChars )
 	{
 	return fread( p, sizeof(char), numChars, pFile_ );
 	} // ~operator()
 	
 	virtual void rewindFile(void);
 	
 	virtual LetterCountType tellg( void ) const;
 	
 	protected:
 	LetterCountType currentPos_;
 	}; // ~class BwtReaderASCII
 	
 	class BwtReaderRunLength : public BwtReaderBase
 	{
 	public:
 	BwtReaderRunLength( const std::string& fileName );
 	
 	virtual int readAndCount( LetterCount& c, const LetterCountType numChars );
 	
 	virtual int readAndSend( BwtWriterBase& writer, const int numChars );
 	
 	virtual int operator()( char* p, int numChars );
 	
 	virtual void rewindFile(void);
 	
 	virtual LetterCountType tellg( void ) const;
 	
 	
 	bool getRun(void);
 	
 	protected:
 	uint lengths_[256];
 	uchar codes_[256];
 	uchar buf_[ReadBufferSize];
 	uint runLength_;
 	uchar* pBuf_;
 	uchar* pBufMax_;
 	uchar lastChar_;
 	bool finished_;
 	LetterCountType currentPos_;
 	
 	}; // class ~BwtReaderRunLength
 	
 	
 	#endif
