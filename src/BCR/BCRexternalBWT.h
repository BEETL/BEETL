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

#ifndef _BCRexternalBWT_H_
#define _BCRexternalBWT_H_

#include <map>
#include "BWTCollection.h"
#include <fstream>
#include <iostream>

//#include "LetterCount.hh" // TEMP

class BCRexternalBWT : public SXSI::BWTCollection {
public:
    /**
     * Constructor
     */
  explicit BCRexternalBWT
    (char* file1, char* fileOut, int mode,
     CompressionFormatType outputCompression );
    ~BCRexternalBWT();
	
	int buildBCR(char const *, char const *);
	int unbuildBCR(char const*, char const*, char const *, char const *);
	int backwardSearchBCR(char const* , char const* , char const * , char const *);
	int decodeBCRnaiveForward(char const *, char const *, char const *);	//Inverse BWT by Forward direction of nText sequences, one sequence at a time, in lexicographic order.
	int decodeBCRmultipleReverse(char const*, char const*, char const *);   //Inverse BWT by Backward direction of nText sequences at the same time by lengthRead iterations.
	int Recover1symbolReverse(char const* , char const* , uchar *, sortElement *);
	int RecoverNsymbolsReverse(char  const*, char const*, uchar *);
	int RecoverNsymbolsReverseByVector(char const* file1, char const* fileOutBwt, uchar * newSymb);
	dataTypeNSeq recover1SequenceForward(char const* , char const * , sortElement , uchar *, dataTypelenSeq *) ;
	vector <int> recoverNSequenceForward(char const* , char const *, dataTypeNSeq);
	int recoverNSequenceForwardSequentially(char const* , char const *, dataTypeNSeq);
	void storeBWT(uchar const *);
	void storeEntireBWT(const char*);
	void storeSA(dataTypelenSeq);
	void storeEntirePairSA(const char*);
	void storeEntireSAfromPairSA(const char*);
	dataTypeNChar rankManySymbols(FILE &, dataTypeNChar *, dataTypeNChar, uchar *);
#ifdef XXX
	dataTypeNChar rankManySymbols(FILE &, LetterCount&, dataTypeNChar, uchar *); // TEMP
#endif
	dataTypeNChar rankManySymbolsByVector(FILE & , dataTypeNChar *, dataTypeNChar, uchar *);
	dataTypeNChar findRankInBWT (char const*, char const*, dataTypedimAlpha, dataTypeNChar, uchar);
	dataTypeNChar findRankInBWTbyVector (char const*, char const*, dataTypedimAlpha, dataTypeNChar, uchar);
	int rankInverseManyByVector (char const* , char const* , dataTypeNSeq , uchar *);
	int backwardSearchManyBCR(char const* , char const*, char const *, vector<string>, dataTypelenSeq);
	int SearchAndLocateKmer (char const* , char const* , char const * , vector<string> , dataTypelenSeq , vector <int>*);
private:  
	void InsertNsymbols(uchar const *, dataTypelenSeq);
	void InsertFirstsymbols(uchar const *);
	int initializeUnbuildBCR(char const*, char const*, dataTypeNChar []);
	int computeNewPositonForBackSearch (char const*, char const*, uchar );
	int computeNewPositonForBackSearchByVector (char const*, char const*, uchar );
	int computeManyNewPositonForBackSearchByVector(char const* , char const* , uchar*, dataTypeNSeq);
	int computeVectorUnbuildBCR(char const*, char const*, dataTypeNChar []);
	int update_Pos_Pile(sortElement *);
	int update_Pos_Pile_Blocks(dataTypeNChar *, dataTypeNChar *, dataTypedimAlpha, uchar);
	int findBlockToRead(dataTypeNChar *, dataTypedimAlpha , dataTypeNChar *, dataTypeNChar *);
};

#endif
