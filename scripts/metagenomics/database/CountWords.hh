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
 	
 	#ifndef INCLUDED_COUNTWORDS_HH
 	#define INCLUDED_COUNTWORDS_HH
 	
 	
 	static const char* idhh="@(#) $Id: CountWords.hh,v 1.11 2011/11/28 09:47:38 rshaw Exp $";
 	
 	using namespace std;
 	
 	
 	#include <string>
 	#include <LetterCount.hh>
#include <Alphabet.hh>
 	
 	//#define DEBUG 1
 	
#define COMPRESS_SEQ 1
	
struct Range
{
  Range(const string& word, LetterCountType pos, LetterCountType num ) :
    word_(word), pos_(pos), num_(num) {}
  Range(void) {}
  string word_;
  LetterCountType pos_;
  LetterCountType num_;
};

struct AllRanges: public vector<vector<vector<Range> > >
{
  AllRanges( void ) : vector<vector<vector<Range > > >
  ( alphabetSize, vector<vector<Range> >(alphabetSize)) {}
};

//
// abstract RangeStore - manage sets of intervals in BWT
//

struct RangeStore
{
  virtual ~RangeStore() {}
  virtual void swap( void )=0;
  virtual void setPortion( int pileNum, int portionNum )=0;
  virtual bool getRange( Range& thisRange )=0;
  virtual void addRange( int pileNum, int portionNum, const string& seq,
			 LetterCountType pos, LetterCountType num )=0;
  
  
}; // ~struct RangeStore

 	//
 	// RangeStoreRAM - hold BWT intervals in guess where
 	//

struct RangeStoreRAM : public RangeStore
{
  RangeStoreRAM() : pThis(&r1), pNext(&r2) {}
  virtual ~RangeStoreRAM() {}
  
  virtual void swap( void )
  {
    pTemp=pNext;pNext=pThis;pThis=pTemp;
  } // ~swap
 	
 	
  virtual void setPortion( int pileNum, int portionNum )
  {
    i_=(*pThis)[pileNum][portionNum].begin();
    end_=(*pThis)[pileNum][portionNum].end();
  } // ~setPortion
  
  virtual bool getRange( Range& thisRange )
  {
    if (i_!=end_)
      {
 	thisRange=*(i_++);
 	return true;
      }
    else return false;
  } // ~getRange
  
  // for (vector<Range>::iterator k((*r.pThis)[i][j].begin());
  // k!=(*r.pThis)[i][j].end();++k)
  // while (r.getNextRange(thisRange))
  
  
  
  virtual void addRange( int pileNum, int portionNum, const string& seq,
			 LetterCountType pos, LetterCountType num )
  {
    (*pNext)[pileNum][portionNum].push_back(Range(seq,pos,num));
#ifdef DEBUG
    cout << "add range: " << alphabet[pileNum] << " " << alphabet[portionNum]
	 << " " << seq << " " << pos << " " << num << endl;
#endif
  } // ~addRange
 	
 	// clear range store ready for next iter
 	virtual void clear( void)
 	{
 	for (int i(1);i<alphabetSize;++i)
 	{
 	for (int j(1);j<alphabetSize;++j)
 	{
 	(*pThis)[i][j].clear();
 	} // ~for j
 	} // ~for i
 	} // ~clear
 	
 	AllRanges r1, r2;
 	AllRanges *pThis,*pNext,*pTemp;
 	vector<Range>::iterator i_;
 	vector<Range>::iterator end_;
 	
 	}; // ~struct RangeStoreRAM
 	
 	/*****************************************************************************/
 	
 	void writeSeq(FILE* outFileHndl, const char* pcSeq);
 	void readSeq(FILE* inFileHndl, char* pcSeq);
 	
 	/*****************************************************************************/
 	
 	//
 	// RangeStoreExternal
 	//
 	
 	typedef unsigned short NumberFrag;
 	const int bitsPerFrag((8*sizeof(NumberFrag))-1);
 	const NumberFrag needAnotherFrag(1<<bitsPerFrag);
 	const NumberFrag fragMask(~needAnotherFrag);
 	const LetterCountType letterCountMask((LetterCountType)fragMask);
 	const int fragBufSize(8);
 	
 	struct RangeState
 	{
 	RangeState() : pFile_(NULL) { clear(); }
 	void clear(void)
 	{
 	// wordLast_.clear();
 	// next line makes sure string comparison fails at first char
 	// wordLast_+=notInAlphabet;
 	wordLast_[0]=notInAlphabet;
 	wordLast_[1]='\0';
 	posLast_=0;
 	if (pFile_!=NULL) fclose(pFile_);
 	pFile_=NULL;
 	} // ~clear
 	void addSeq( const string& seq )
 	{
 	if (wordLast_[0]!=notInAlphabet)
 	{
 	const char *pSeq(seq.c_str()), *pLast(wordLast_);
 	while (*pSeq==*pLast) { ++pSeq; ++pLast; }
 	#ifdef DEBUG
 	cout << "FF " << wordLast_ << " " << seq << " " << pSeq << endl;
 	#endif
 	
 	#ifdef COMPRESS_SEQ
 	writeSeq(pFile_, pSeq);
 	#else
 	fprintf(pFile_,"%s\n",pSeq);
 	// fprintf(pFile_,"%s\n",pSeq);
 	#endif
 	
 	strcpy(wordLast_,seq.c_str());
 	// strcpy(wordLast_,pSeq);
 	}
 	else
 	#ifdef COMPRESS_SEQ
 	writeSeq(pFile_, seq.c_str());
 	#else
 	fprintf(pFile_,"%s\n",seq.c_str());
 	#endif
 	}
 	
 	void getSeq( string& word )
 	{
 	// assert(fgets(buf_,256,stateIn_.pFile_)!=NULL);
 	// buf_[strlen(buf_)-1]='\0';
 	// thisRange.word_=buf_;
 	
 	
 	if (wordLast_[0]!=notInAlphabet)
 	{
 	word=wordLast_;
 	#ifdef DEBUG
 	cout << "FF " << word;
 	#endif
 	
 	#ifdef COMPRESS_SEQ
 	readSeq(pFile_, wordLast_); // readSeq adds \0.
 	#else
 	assert(fgets(wordLast_,256,pFile_)!=NULL);
 	wordLast_[strlen(wordLast_)-1]='\0'; // Replace \n with \0.
 	#endif
 	
 	for (int i(0);i<strlen(wordLast_);i++)
 	word[word.size()-strlen(wordLast_)+i]=wordLast_[i];
 	#ifdef DEBUG
 	cout << " -> " << word << endl;
 	#endif
 	}
 	else
 	{
 	#ifdef COMPRESS_SEQ
 	readSeq(pFile_, wordLast_); // readSeq adds \0.
 	#else
 	assert(fgets(wordLast_,256,pFile_)!=NULL);
 	wordLast_[strlen(wordLast_)-1]='\0'; // Replace \n with \0
 	#endif
 	
 	word=wordLast_;
 	}
 	
 	
 	// assert(fgets(wordLast_,256,pFile_)!=NULL);
 	// wordLast_[strlen(wordLast_)-1]='\0';
 	// word=wordLast_;
 	}
 	
 	
 	
 	void addNum( LetterCountType num )
 	{
 	#ifdef DEBUG
 	cout << "AN: send " << num << endl;
 	#endif
 	NumberFrag* pFrag(fragBuf_+fragBufSize);
 	do
 	{
 	pFrag--;
 	*pFrag=(NumberFrag)(num&letterCountMask);
 	num>>=bitsPerFrag;
 	#ifdef DEBUG
 	cout << "AN: " << *pFrag << endl;
 	#endif
 	} while (num!=0);
 	fragBuf_[fragBufSize-1]|=needAnotherFrag;
 	const int toWrite((fragBuf_+fragBufSize)-pFrag);
 	#ifdef DEBUG
 	cout << "AN: sent " << toWrite << endl;
 	#endif
 	assert(fwrite(pFrag,sizeof(NumberFrag),toWrite,
 	pFile_)==toWrite);
 	}
 	
 	bool getNum( LetterCountType& num )
 	{
 	
 	num=0;
 	if (fread(fragBuf_,sizeof(NumberFrag),1,pFile_)!=1)
 	return false;
 	#ifdef DEBUG
 	cout << "AN: got " << fragBuf_[0] << endl;
 	#endif
 	while ((fragBuf_[0]&needAnotherFrag)==0)
 	{
 	num|=(LetterCountType)fragBuf_[0];
 	num<<=bitsPerFrag;
 	assert (fread(fragBuf_,sizeof(NumberFrag),1,pFile_)==1);
 	#ifdef DEBUG
 	cout << "AN: get " << fragBuf_[0] << endl;
 	#endif
 	} // ~while
 	
 	fragBuf_[0]&=fragMask;
 	num|=(LetterCountType)fragBuf_[0];
 	#ifdef DEBUG
 	cout << "AN: decoded " << fragBuf_[0] << endl;
 	#endif
 	return true;
 	}
 	
 	
 	
 	NumberFrag fragBuf_[fragBufSize];
 	char wordLast_[256];
 	LetterCountType posLast_;
 	FILE* pFile_;
 	}; // ~struct RangeState
 	
 	
 	//
 	// RangeStoreExternal
 	//
 	
 	struct RangeStoreExternal : public RangeStore
 	{
 	RangeStoreExternal( const string fileStemIn="countA",
 	const string fileStemOut="countB") :
 	fileStemIn_(fileStemIn),fileStemOut_(fileStemOut),stateIn_()
 	{
 	string fileName;
 	for (int i(1);i<alphabetSize;++i)
 	{
 	for (int j(1);j<alphabetSize;++j)
 	{
 	stateOut_[i][j].clear();
 	// stateOut_[i][j].pFile_=NULL;
 	// (*pThis)[i][j].clear();
 	getFileName(fileStemIn_,i,j,fileName);
 	if (remove(fileName.c_str())==0)
 	{
 	#ifdef DEBUG
 	cerr << "Removed " << fileName << endl;
 	#endif
 	}
 	getFileName(fileStemOut_,i,j,fileName);
 	if (remove(fileName.c_str())==0)
 	{
 	#ifdef DEBUG
 	cerr << "Removed " << fileName << endl;
 	#endif
 	}
 	} // ~for j
 	} // ~for i
 	} // ~ctor
 	
 	virtual ~RangeStoreExternal() {}
 	
 	virtual void swap( void )
 	{
 	#ifdef DEBUG
 	cout << "swap" << endl;
 	cout << fileStemIn_ << " " << fileStemOut_ << endl;
 	#endif
 	string temp=fileStemIn_;
 	fileStemIn_=fileStemOut_;
 	fileStemOut_=temp;
 	#ifdef DEBUG
 	cout << fileStemIn_ << " " << fileStemOut_ << endl;
 	#endif
 	} // ~swap
 	
 	virtual void setPortion( int pileNum, int portionNum )
 	{
 	#ifdef DEBUG
 	cout << "set portion " << alphabet[pileNum] << alphabet[portionNum] << endl;
 	#endif
 	stateIn_.clear();
 	// if (stateIn_.pFile_!=NULL) fclose(stateIn_.pFile_);
 	string fileName;
 	getFileName(fileStemIn_,pileNum,portionNum,fileName);
 	#ifdef DEBUG
 	cout << "Made input file name " << fileName << endl;
 	#endif
 	stateIn_.pFile_=fopen(fileName.c_str(),"r");
 	#ifdef DEBUG
 	if (stateIn_.pFile_==NULL)
 	{
 	cerr << "Warning: no file " << fileName
 	<< " found, presuming no ranges of interest in this region"
 	<< endl;
 	} // ~if
 	#endif
 	}
 	
 	virtual bool getRange( Range& thisRange )
 	{
 	#ifdef DEBUG
 	cout << "get range: " << fileStemIn_ << endl;
 	#endif
 	if (stateIn_.pFile_==NULL)
 	return false;
 	// else if (fread(&thisRange.pos_,sizeof(LetterCountType),1,
 	// stateIn_.pFile_)!=1)
 	else if (stateIn_.getNum(thisRange.pos_)==false)
 	{
 	fclose(stateIn_.pFile_);stateIn_.pFile_=NULL;
 	return false;
 	
 	}
 	else
 	{
 	// assert(fread(&thisRange.num_,sizeof(LetterCountType),1,
 	// stateIn_.pFile_)==1);
 	assert(stateIn_.getNum(thisRange.num_));
 	
 	stateIn_.getSeq( thisRange.word_ );
 	// assert(fgets(buf_,256,stateIn_.pFile_)!=NULL);
 	// buf_[strlen(buf_)-1]='\0';
 	// thisRange.word_=buf_;
 	
 	
 	#ifdef DEBUG
 	cout << "got range: " << fileStemIn_ << " " << thisRange.word_ << " " << thisRange.pos_
 	<< " " << thisRange.num_ << endl;
 	#endif
 	return true;
 	}
 	} // ~getRange
 	
 	virtual void addRange( int pileNum, int portionNum, const string& seq,
 	LetterCountType pos, LetterCountType num )
 	{
 	#ifdef DEBUG
 	cout << "set range: " << fileStemIn_ << " " << alphabet[pileNum] << " " << alphabet[portionNum]
 	<< " " << seq << " " << pos << " " << num << endl;
 	#endif
 	if (stateOut_[pileNum][portionNum].pFile_==NULL)
 	{
 	string fileName;
 	getFileName(fileStemOut_,pileNum,portionNum,fileName);
 	#ifdef DEBUG
 	cout << "Made output file name " << fileName << endl;
 	#endif
 	stateOut_[pileNum][portionNum].pFile_=fopen(fileName.c_str(),"w");
 	assert(stateOut_[pileNum][portionNum].pFile_!=NULL);
 	} // ~if
 	// assert(fwrite(&pos,sizeof(LetterCountType),1,
 	// stateOut_[pileNum][portionNum].pFile_)==1);
 	// assert(fwrite(&num,sizeof(LetterCountType),1,
 	// stateOut_[pileNum][portionNum].pFile_)==1);
 	stateOut_[pileNum][portionNum].addNum(pos);
 	stateOut_[pileNum][portionNum].addNum(num);
 	
 	
 	#ifdef DEBUG
 	cout << "add seq " << seq.c_str() << " " << strlen(seq.c_str()) << endl;
 	#endif
 	// fputs(seq.c_str(),fileOut_[pileNum][portionNum]);
 	// fprintf(stateOut_[pileNum][portionNum].pFile_,"%s\n",seq.c_str());
 	stateOut_[pileNum][portionNum].addSeq(seq);
 	}
 	
 	virtual void clear( void)
 	{
 	
 	string fileName;
 	
 	for (int i(1);i<alphabetSize;++i)
 	{
 	for (int j(1);j<alphabetSize;++j)
 	{
 	// if (stateOut_[i][j].pFile_!=NULL)
 	// fclose( stateOut_[i][j].pFile_ );
 	// stateOut_[i][j].pFile_=NULL;
 	stateOut_[i][j].clear();
 	
 	getFileName(fileStemIn_,i,j,fileName);
 	if (remove(fileName.c_str())!=0)
 	{
 	#ifdef DEBUG
 	cerr << "Could not remove " << fileName << endl;
 	#endif
 	}
 	
 	// (*pThis)[i][j].clear();
 	} // ~for j
 	} // ~for i
 	} // ~clear
 	
 	
 	
 	
 	void getFileName( const string& stem, const int pile, const int portion,
 	string& fileName )
 	{
 	fileName=stem;
 	fileName+='-';
 	fileName+='0';
 	assert(pile<=9);
 	fileName+=(char)(48+pile);
 	fileName+='-';
 	fileName+='0';
 	assert(portion<=9);
 	fileName+=(char)(48+portion);
 	
 	}
 	string fileStemIn_;
 	string fileStemOut_;
 	
 	RangeState stateIn_;
 	RangeState stateOut_[alphabetSize][alphabetSize];
 	char buf_[256]; // get rid, not nice
 	
 	}; // ~struct RangeStoreExternal
 	
 	//
 	// IntervalHandler
 	//
 	// Idea here is that different algorithms can be implemented by defining
 	// new subclasses of IntervalHandler
 	// BackTracker class implements the backward search and takes its
 	// cue from IntervalHandler as to whether to continue down a particular
 	// branch of the search tree
 	
 	typedef bool AlphabetFlag[alphabetSize];
 	
 	struct IntervalHandlerBase
 	{
 	virtual ~IntervalHandlerBase() {}
 	virtual void foundInBoth
 	( const LetterCount& countsThisRangeA, const LetterCount& countsThisRangeB,
 	const Range& thisRangeA, const Range& thisRangeB,
 	AlphabetFlag& propagateIntervalA, AlphabetFlag& propagateIntervalB)=0;
 	virtual void foundInAOnly
 	( const LetterCount& countsThisRangeA,
 	const Range& thisRangeA,
 	AlphabetFlag& propagateIntervalA )=0;
 	};
 	
 	struct IntervalHandlerSplice : public IntervalHandlerBase
 	{
 	IntervalHandlerSplice(unsigned int minOcc) : minOcc_(minOcc) {}
 	virtual ~IntervalHandlerSplice() {}
 	virtual void foundInBoth
 	( const LetterCount& countsThisRangeA, const LetterCount& countsThisRangeB,
 	const Range& thisRangeA, const Range& thisRangeB,
 	AlphabetFlag& propagateIntervalA, AlphabetFlag& propagateIntervalB);
 	const LetterCountType minOcc_;
 	virtual void foundInAOnly
 	( const LetterCount& countsThisRangeA,
 	const Range& thisRangeA,
 	AlphabetFlag& propagateIntervalA );
 	};
 	
 	struct IntervalHandlerReference : public IntervalHandlerBase
 	{
 	IntervalHandlerReference(unsigned int minOcc) : minOcc_(minOcc) {}
 	virtual ~IntervalHandlerReference() {}
 	virtual void foundInBoth
 	( const LetterCount& countsThisRangeA, const LetterCount& countsThisRangeB,
 	const Range& thisRangeA, const Range& thisRangeB,
 	AlphabetFlag& propagateIntervalA, AlphabetFlag& propagateIntervalB);
 	const LetterCountType minOcc_;
 	virtual void foundInAOnly
 	( const LetterCount& countsThisRangeA,
 	const Range& thisRangeA,
 	AlphabetFlag& propagateIntervalA );
 	};
 	
 	
 	struct BackTracker
 	{
 	BackTracker( BwtReaderBase* _inBwtA, BwtReaderBase* _inBwtB,
 	LetterCountType& _currentPosA, LetterCountType& _currentPosB,
 	RangeStoreExternal& _rA, RangeStoreExternal& _rB,
 	LetterCount& _countsSoFarA, LetterCount& _countsSoFarB,
 	int _minOcc ) :
 	inBwtA(_inBwtA),inBwtB(_inBwtB),
 	currentPosA(_currentPosA), currentPosB(_currentPosB),
 	rA(_rA),rB(_rB),countsSoFarA(_countsSoFarA),countsSoFarB(_countsSoFarB),
 	minOcc(_minOcc),
 	numRanges(0),numSingletonRanges(0)
 	// intervalHandler_(_minOcc)
 	{}
 	
 	void skipIfNecessary( const Range& thisRange,
 	LetterCountType& currentPos,
 	BwtReaderBase& inBwt,
 	LetterCount& countsSoFar );
 	
 	template<class IntervalHandlerParam> void operator()
 	( int i, string& thisWord, IntervalHandlerParam& intervalHandler_ );
 	
 	BwtReaderBase* inBwtA;
 	BwtReaderBase* inBwtB;
 	
 	LetterCountType& currentPosA;
 	LetterCountType& currentPosB;
 	
 	RangeStoreExternal& rA;
 	RangeStoreExternal& rB;
 	LetterCount& countsSoFarA;
 	LetterCount& countsSoFarB;
 	
 	int minOcc;
 	
 	LetterCountType numRanges;
 	LetterCountType numSingletonRanges;
 	
 	AlphabetFlag propagateIntervalA_;
 	AlphabetFlag propagateIntervalB_;
 	
 	// IntervalHandlerSplice intervalHandler_;
 	
 	};
 	
 	
 	
 	#endif
