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
#include "libzoo/cli/ToolParameters.hh"
#include "parameters/SearchParameters.hh"
#include "shared/FragmentedVector.hh"

#include <fstream>
#include <iostream>
#include <map>



class BCRexternalBWT : public SXSI::BWTCollection
{
public:
    /**
     * Constructor
     */
    explicit BCRexternalBWT ( char *file1, char *fileOut, int mode, CompressionFormatType outputCompression, ToolParameters *toolParams = NULL );
    ~BCRexternalBWT();

    int buildBCR( char const *, char const *, const BwtParameters *bwtParams );
    int unbuildBCR( char const *, char const *, char const *, char const * );
    int backwardSearchBCR( char const * , char const * , char const * , char const * );
    int decodeBCRnaiveForward( char const *, char const *, char const * ); //Inverse BWT by Forward direction of nText sequences, one sequence at a time, in lexicographic order.
    int decodeBCRmultipleReverse( char const *, char const *, char const *, bool processQualities = false ); //Inverse BWT by Backward direction of nText sequences at the same time by lengthRead iterations.
    int Recover1symbolReverse( char const * , char const * , uchar *, sortElement * );
    int RecoverNsymbolsReverse( char  const *, char const *, uchar *, uchar *newQual = 0 );
    int RecoverNsymbolsReverseByVector( char const *file1, char const *fileOutBwt, uchar *newSymb, uchar *newQual = 0 );
    SequenceNumber recover1SequenceForward( char const * , char const * , sortElement , uchar *, SequenceLength * ) ;
    vector <int> recoverNSequenceForward( char const * , char const *, SequenceNumber );
    int recoverNSequenceForwardSequentially( char const * , char const *, SequenceNumber );
    void storeBWT( uchar const *, uchar const *qual = NULL );
    void storeBWT_parallelPile( uchar const *newSymb, uchar const *newQual, unsigned int parallelPile, SequenceNumber startIndex, SequenceNumber endIndex );
    void storeEntireBWT( const char * );
    void storeSA( SequenceLength );
    void storeEntirePairSA( const char * );
    void storeEntireSAfromPairSA( const char * );
    virtual void storeBWTandLCP( uchar const * );
    virtual void storeEntireLCP( const char * );
    LetterNumber rankManySymbols( FILE &, LetterNumber *, LetterNumber, uchar * );
#ifdef XXX
    LetterNumber rankManySymbols( FILE &, LetterCount &, LetterNumber, uchar * ); // TEMP
#endif
    LetterNumber rankManySymbolsByVector( FILE & , LetterNumber *, LetterNumber, uchar *, uchar *foundQual = 0, FILE *InFileBWTQual = 0 );
    LetterNumber findRankInBWT ( char const *, char const *, AlphabetSymbol, LetterNumber, uchar );
    LetterNumber findRankInBWTbyVector ( char const *, char const *, AlphabetSymbol, LetterNumber, uchar );
    int rankInverseManyByVector ( char const * , char const * , SequenceNumber , uchar * );
    int backwardSearchManyBCR( char const * , char const *, char const *, vector<string>, SequenceLength );
    int SearchAndLocateKmer ( char const * , char const * , char const * , vector<string> , SequenceLength , vector <int> & );
private:
    void InsertNsymbols( uchar const *, SequenceLength, uchar const *qual = NULL );
    void InsertNsymbols_parallelPile( uchar const *newSymb, SequenceLength posSymb, uchar const *newQual, unsigned int parallelPile, SequenceNumber startIndex, SequenceNumber endIndex, vector< FragmentedVector< sortElement > > &newVectTriplePerNewPile );
    void InitialiseTmpFiles();
    void InsertFirstsymbols( uchar const *, uchar const *qual = NULL, const int cycle = 0 );
    int initializeUnbuildBCR( char const *, char const *, LetterNumber [] );
    int computeNewPositionForBackSearch ( char const *, char const *, uchar );
    int computeNewPositionForBackSearchByVector ( char const *, char const *, uchar );
    int computeManyNewPositionForBackSearchByVector( char const * , char const * , uchar *, SequenceNumber );
    int computeVectorUnbuildBCR( char const *, char const *, LetterNumber [] );
    int update_Pos_Pile( sortElement * );
    int update_Pos_Pile_Blocks( LetterNumber *, LetterNumber *, AlphabetSymbol, uchar );
    int findBlockToRead( LetterNumber *, AlphabetSymbol , LetterNumber *, LetterNumber * );

private:
    void pauseBetweenCyclesIfNeeded();
    void convertFileFromIntermediateToFinalFormat( const char *filenameIn, const char *filenameOut );
    void ReadFilesForCycle( const char *prefix, const SequenceLength cycle, const SequenceLength readLength, const SequenceNumber nText, uchar *newSymb, const bool processQualities, uchar *newQual );
    BwtReaderBase *instantiateBwtReaderForIntermediateCycle( const char *filenameIn, bool allowDefrag = false );
    BwtWriterBase *instantiateBwtWriterForIntermediateCycle( const char *filenameOut );
    BwtWriterBase *instantiateBwtWriterForLastCycle( const char *filenameOut );
    BwtWriterBase *pWriterBwt0_; // persistent file, as we only ever need to append (never insert) characters to it
    ToolParameters *toolParams_;
    BwtParameters *bwtParams_;
    UnbwtParameters *unbwtParams_;
    SearchParameters *searchParams_;
};

#endif
