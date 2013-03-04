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

#include "BWTCollection.hh"
#include "BwtReader.hh"
#include "parameters/BwtParameters.hh"

#include <fstream>
#include <iostream>
#include <map>


class BCRexternalBWT : public SXSI::BWTCollection
{
public:
    /**
     * Constructor
     */
    explicit BCRexternalBWT ( char *file1, char *fileOut, int mode, CompressionFormatType outputCompression, const BwtParameters *bwtParams = NULL );
    ~BCRexternalBWT();

    int buildBCR( char const *, char const *, const BwtParameters *bwtParams );
    int unbuildBCR( char const *, char const *, char const *, char const * );
    int backwardSearchBCR( char const * , char const * , char const * , char const * );
    int decodeBCRnaiveForward( char const *, char const *, char const * ); //Inverse BWT by Forward direction of nText sequences, one sequence at a time, in lexicographic order.
    int decodeBCRmultipleReverse( char const *, char const *, char const *, bool processQualities = false ); //Inverse BWT by Backward direction of nText sequences at the same time by lengthRead iterations.
    int Recover1symbolReverse( char const * , char const * , uchar *, sortElement * );
    int RecoverNsymbolsReverse( char  const *, char const *, uchar *, uchar *newQual = 0 );
    int RecoverNsymbolsReverseByVector( char const *file1, char const *fileOutBwt, uchar *newSymb, uchar *newQual = 0 );
    dataTypeNSeq recover1SequenceForward( char const * , char const * , sortElement , uchar *, dataTypelenSeq * ) ;
    vector <int> recoverNSequenceForward( char const * , char const *, dataTypeNSeq );
    int recoverNSequenceForwardSequentially( char const * , char const *, dataTypeNSeq );
    void storeBWT( uchar const *, uchar const *qual = NULL );
    void storeBWT_parallelPile( uchar const *newSymb, uchar const *newQual, unsigned int parallelPile, dataTypeNSeq startIndex, dataTypeNSeq endIndex );
    void storeEntireBWT( const char * );
    void storeSA( dataTypelenSeq );
    void storeEntirePairSA( const char * );
    void storeEntireSAfromPairSA( const char * );
    dataTypeNChar rankManySymbols( FILE &, dataTypeNChar *, dataTypeNChar, uchar * );
#ifdef XXX
    dataTypeNChar rankManySymbols( FILE &, LetterCount &, dataTypeNChar, uchar * ); // TEMP
#endif
    dataTypeNChar rankManySymbolsByVector( FILE & , dataTypeNChar *, dataTypeNChar, uchar *, uchar *foundQual = 0, FILE *InFileBWTQual = 0 );
    dataTypeNChar findRankInBWT ( char const *, char const *, dataTypedimAlpha, dataTypeNChar, uchar );
    dataTypeNChar findRankInBWTbyVector ( char const *, char const *, dataTypedimAlpha, dataTypeNChar, uchar );
    int rankInverseManyByVector ( char const * , char const * , dataTypeNSeq , uchar * );
    int backwardSearchManyBCR( char const * , char const *, char const *, vector<string>, dataTypelenSeq );
    int SearchAndLocateKmer ( char const * , char const * , char const * , vector<string> , dataTypelenSeq , vector <int> * );
private:
    void InsertNsymbols( uchar const *, dataTypelenSeq, uchar const *qual = NULL );
    void InsertNsymbols_parallelPile( uchar const *newSymb, dataTypelenSeq posSymb, uchar const *newQual, unsigned int parallelPile, dataTypeNSeq startIndex, dataTypeNSeq endIndex, vector< vector< sortElement > > &newVectTriplePerNewPile );
    void InsertFirstsymbols( uchar const *,uchar const *qual = NULL );
    int initializeUnbuildBCR( char const *, char const *, dataTypeNChar [] );
    int computeNewPositonForBackSearch ( char const *, char const *, uchar );
    int computeNewPositonForBackSearchByVector ( char const *, char const *, uchar );
    int computeManyNewPositonForBackSearchByVector( char const * , char const * , uchar *, dataTypeNSeq );
    int computeVectorUnbuildBCR( char const *, char const *, dataTypeNChar [] );
    int update_Pos_Pile( sortElement * );
    int update_Pos_Pile_Blocks( dataTypeNChar *, dataTypeNChar *, dataTypedimAlpha, uchar );
    int findBlockToRead( dataTypeNChar *, dataTypedimAlpha , dataTypeNChar *, dataTypeNChar * );

private:
    void dumpRamFileToFile( const char *filenameIn, const char *filenameOut );
    void ReadFilesForCycle( const char *prefix, const dataTypelenSeq cycle, const dataTypeNSeq nText, uchar *newSymb, const bool processQualities, uchar *newQual );
    BwtReaderBase *instantiateBwtReaderForFirstCycle( const char *filenameIn );
    BwtWriterBase *instantiateBwtWriterForFirstCycle( const char *filenameOut );
    BwtReaderBase *instantiateBwtReaderForIntermediateCycle( const char *filenameIn, bool allowDefrag = false );
    BwtWriterBase *instantiateBwtWriterForIntermediateCycle( const char *filenameOut );
    BwtWriterBase *instantiateBwtWriterForLastCycle( const char *filenameOut );
    const BwtParameters *bwtParams_;
};

#endif
