//
 	// (c) 2011 Illumina
 	//
 	// Encapsulate the functionality for rank and inverse-rank queries;
 	// operates on a single BWT file (agnostic to whether it's a partial
 	// or complete BWT)
 	
 	#include <cstdio>
 	#include <iostream>
 	#include <fstream>
 	#include <string>
 	#include <vector>
 	#include <map>
 	#include <assert.h>
 	#include <algorithm>
 	
 	// data type definitions
 	#include "Tools.h"
 	
 	
 	using namespace std;
 	
 	
 	/* class OccComparison */
 	/* { */
 	/* public: */
 	/* static int idx; */
 	
 	/* static struct _CompareFloatField */
 	/* { */
 	/* bool operator() (const vector<dataTypeNSeq>& left, int right) */
 	/* { */
 	/* return left[1]<right; */
 	/* } */
 	
 	/* } CompareOcc; */
 	/* }; */
 	
 	
 	class StaticOccComparison
 	{
 	public:
 	int idx;
 	
 	StaticOccComparison(const int& i ) { idx=i; }
 	
 	bool operator() (const vector<dataTypeNSeq>& left, dataTypeNSeq right)
 	{
 	return left[idx]<right;
 	}
 	};
 	
 	
 	
 	class Rank{
 	public:
 	Rank( void );
 	~Rank();
 	
 	// set methods
 	bool setAlpha( dataTypedimAlpha *alpha,const int& s );
 	bool setAlphaInverse( dataTypedimAlpha *alphaInverse,const int& s );
 	void setCharCount( const int& cnt ) { charCount_=cnt; }
 	
 	// read BWT and initialize all the counts
 	bool initialize( const string& fn, const int& blockSize );
 	
 	// Given character x and a number n, get the blockID that contains
 	// the n-th occurrence of x
 	dataTypeNChar getInverseRank( const char& c,const dataTypeNChar& n );
 	
 	private:
 	FILE* ifsBWT; // input file streamn of BWT
 	
 	vector<dataTypedimAlpha> alpha_;
 	vector<dataTypedimAlpha> alphaInverse_;
 	bool alphaInitialized_;
 	bool alphaInverseInitialized_;
 	int charCount_;
 	int blockSize_;
 	
 	// holding the actual counts
 	vector< vector<dataTypeNSeq> > occ_; // addressing: occ_[blockID][char]
 	vector< map<dataTypedimAlpha,dataTypedimAlpha> > mapCountToBlock_; // adressing mapCountToBlock_[char][countToFind]==blockID
 	vector< StaticOccComparison > compareObjects_;
 	
 	};
