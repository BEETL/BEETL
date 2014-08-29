/**
 ** Copyright (c) 2011-2014 Illumina, Inc.
 **
 ** This file is part of the BEETL software package,
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
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
#include "shared/FragmentedVector.hh"

#include <fstream>
#include <iostream>
#include <map>


class SearchParameters;


class BCRexternalBWT : public SXSI::BWTCollection
{
public:
    /**
     * Constructor
     */
    explicit BCRexternalBWT ( const string &file1, const string &fileOut, const int mode, const CompressionFormatType outputCompression, ToolParameters *toolParams = NULL );
    ~BCRexternalBWT();

    int buildBCR( const string &, const string &, const BwtParameters *bwtParams );
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
    void storeEntireBWT( const string & );
    void storeSA( SequenceLength );
    void storeEntirePairSA( const char * );
    void storeEntireSAfromPairSA( const char * );
    virtual void storeBWTandLCP( uchar const * );
    virtual void storeEntireLCP( const string & );
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
    void InsertFirstsymbols( uchar const *, uchar const *qual = NULL, const int subSequenceNum = 0 );
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
    BwtReaderBase *instantiateBwtReaderForLastCycle( const char *filenameOut );
    void writeEndPosFile( const uint8_t subSequenceNum, const bool lastFile );

    BwtWriterBase *pWriterBwt0_; // persistent file, as we only ever need to append (never insert) characters to it
    shared_ptr< ToolParameters > toolParams_;
    shared_ptr< BwtParameters > bwtParams_;
    shared_ptr< UnbwtParameters > unbwtParams_;
    shared_ptr< SearchParameters > searchParams_;
};

#endif
