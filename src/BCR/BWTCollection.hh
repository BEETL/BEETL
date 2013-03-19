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

#ifndef _SXSI_BWTCollection_h_
#define _SXSI_BWTCollection_h_

#include "Algorithm.hh"
#include "LetterCount.hh"
#include "Sorting.hh"
#include "Tools.hh" // Defines ulong and uchar.
#include "parameters/BwtParameters.hh"
#include "parameters/UnbwtParameters.hh"

#include <iostream>


// by Tobias, small class interface to call from beetl executable

class BCR : public Algorithm
{
    int mode_;
    char *inFile_;
    char *outFile_;
    CompressionFormatType outputCompression_;
public:
    BCR( int mode, string in, string out,
         CompressionFormatType outputCompression );
    ~BCR()
    {
        delete[] inFile_;
        delete outFile_;
    }
    void run( void );
};
// end


namespace SXSI
{
/**
 * General interface for a bwt collection
 *
 * Class is virtual, make objects by calling
 * the static method InitBWTCollection().
 */
class BWTCollection
{
public:

    vector <sortElement> vectTriple;  //Is is used both encoding, decoding, searching.
    //ulong seqN;  //contains a number of a sequence
    //ulong posN;  //contains the position of the last inserted symbol of the sequence seqN[i]
    //uchar pileN; //contains the number of the pile of the last inserted symbol of the sequence seqN[i]

    dataTypeNSeq nText;  //number total of texts in filename1
    //dataTypeNSeq middle; // number of sequence in filename1
    //dataTypeNChar middleLength; //text[middleLength] = the first letter of the second database (filename2)
    dataTypelenSeq lengthRead; //number of char in each text + $
    dataTypeNChar lengthTot;   //length of the all texts without $
    dataTypeNChar lengthTot_plus_eof; //length of the BWT
    dataTypeNChar **tableOcc; //contains the number of occurrences of each symbol
    LetterCountEachPile tableOcc_; // replace tableOcc
    dataTypedimAlpha alpha[256]; //Corresponding between the alphabet, the piles and tableOcc
    dataTypedimAlpha sizeAlpha;  //number of the different symbols in the input texts
    dataTypedimAlpha *alphaInverse;  //Corresponding between alpha[i] and the symbol as char

    vector< vector< vector<dataTypeNChar> > > vectorOcc;
    vector <dataTypeNChar> numBlocksInPartialBWT;

    //dataTypeNChar*** vectorOcc;
    vector<sortElement> FirstVector, LastVector;

    CompressionFormatType outputCompression_;


    /**
    * Init an instance of a text collection object
    *
    * Returns a pointer to an object implementing this interface.
    */
    static BWTCollection *InitBWTCollection
    ( char *file1, char *fileOut, int mode,
      CompressionFormatType outputCompression );

    /**
     * Virtual destructor
     */
    virtual ~BWTCollection() { };
    /**
     *
     * The i'th text insertion gets an identifier value i-1.
     * In other words, document identifiers start from 0.
     */
    virtual int buildBCR( char const *, char const *, const BwtParameters *bwtParams ) = 0;
    virtual int unbuildBCR( char const *, char const *, char const *, char const * ) = 0;
    virtual int backwardSearchBCR( char const * , char const * , char const * , char const * ) = 0;
    virtual int decodeBCRnaiveForward( char const *, char const *, char const * ) = 0; //Inverse BWT by Forward direction of nText sequences, one sequence at a time, in lexicographic order.
    virtual int decodeBCRmultipleReverse( char const *, char const *, char const *, bool processQualities = false ) = 0; //Inverse BWT by Backward direction of nText sequences at the same time by lengthRead iterations.
    virtual int RecoverNsymbolsReverse( char const *, char const *, uchar *, uchar *newQual = 0 ) = 0;
    virtual int RecoverNsymbolsReverseByVector( char const *file1, char const *fileOutBwt, uchar *newSymb, uchar *newQual = 0 ) = 0;
    virtual int Recover1symbolReverse( char const * , char const * , uchar *, sortElement * ) = 0;
    virtual dataTypeNSeq recover1SequenceForward( char const * , char const * , sortElement , uchar *, dataTypelenSeq * ) = 0 ;
    virtual vector <int> recoverNSequenceForward( char const * , char const *, dataTypeNSeq ) = 0;
    virtual int recoverNSequenceForwardSequentially( char const * , char const *, dataTypeNSeq ) = 0;
    virtual void storeBWT( uchar const *, uchar const *qual = NULL ) = 0;
    virtual void storeEntireBWT( const char * ) = 0;
    virtual void storeSA( dataTypelenSeq ) = 0;
    virtual void storeEntirePairSA( const char * ) = 0;
    virtual void storeEntireSAfromPairSA( const char * ) = 0;
    virtual void storeBWTandLCP( uchar const * ) = 0;
    virtual void storeEntireLCP( const char * ) = 0;
    virtual dataTypeNChar rankManySymbols( FILE &, dataTypeNChar *, dataTypeNChar, uchar * ) = 0;
    virtual dataTypeNChar rankManySymbolsByVector( FILE & , dataTypeNChar *, dataTypeNChar, uchar *, uchar *foundQual = 0, FILE *InFileBWTQual = 0 ) = 0;
    virtual dataTypeNChar findRankInBWT ( char const *, char const *, dataTypedimAlpha, dataTypeNChar, uchar ) = 0;
    virtual dataTypeNChar findRankInBWTbyVector ( char const *, char const *, dataTypedimAlpha, dataTypeNChar, uchar ) = 0;
    virtual int rankInverseManyByVector ( char const * , char const * , dataTypeNSeq , uchar * ) = 0;
    virtual int backwardSearchManyBCR( char const * , char const *, char const *, vector<string>, dataTypelenSeq ) = 0;
    virtual int SearchAndLocateKmer ( char const * , char const * , char const * , vector<string> , dataTypelenSeq, vector <int> * ) = 0;
private:

    virtual void InsertNsymbols( uchar const *, dataTypelenSeq, uchar const *qual = NULL ) = 0;
    virtual void InsertFirstsymbols( uchar const *, uchar const *qual = NULL ) = 0;
    virtual int initializeUnbuildBCR( char const *, char const *, dataTypeNChar [] ) = 0;
    virtual int computeNewPositonForBackSearch ( char const *, char const *, uchar ) = 0;
    virtual int computeNewPositonForBackSearchByVector ( char const *, char const *, uchar ) = 0;
    virtual int computeVectorUnbuildBCR( char const *, char const *, dataTypeNChar [] ) = 0;
    virtual int computeManyNewPositonForBackSearchByVector( char const * , char const * , uchar *, dataTypeNSeq ) = 0;
    virtual int update_Pos_Pile( sortElement * ) = 0;
    virtual int update_Pos_Pile_Blocks( dataTypeNChar *, dataTypeNChar *, dataTypedimAlpha, uchar ) = 0;
    virtual int findBlockToRead( dataTypeNChar *, dataTypedimAlpha , dataTypeNChar *, dataTypeNChar * ) = 0;
protected:
    // Protected constructor; call the static function InitBWTCollection().
    BWTCollection() { };

    // No copy constructor or assignment
    BWTCollection( BWTCollection const & );
    BWTCollection &operator = ( BWTCollection const & );
};
}
#endif
