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

#ifndef _SXSI_BWTCollection_h_
#define _SXSI_BWTCollection_h_

#include "Algorithm.hh"
#include "LetterCount.hh"
#include "Sorting.hh"
#include "Tools.hh" // Defines ulong and uchar.
#include "parameters/BwtParameters.hh"
#include "parameters/UnbwtParameters.hh"

#include <iostream>
#include <string>

using std::string;


// by Tobias, small class interface to call from beetl executable

class BCR : public Algorithm
{
    int mode_;
    string inFile_;
    string outFile_;
    CompressionFormatType outputCompression_;
public:
    BCR( const int mode, const string &in, const string &out,
         const CompressionFormatType outputCompression );
    ~BCR()
    {
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

    SequenceNumber nText;  //number total of texts in filename1
    //SequenceNumber middle; // number of sequence in filename1
    //LetterNumber middleLength; //text[middleLength] = the first letter of the second database (filename2)
    SequenceLength lengthRead; //number of char in each text + $
    LetterNumber lengthTot;   //length of the all texts without $
    LetterNumber lengthTot_plus_eof; //length of the BWT
    LetterNumber **tableOcc; //contains the number of occurrences of each symbol
    LetterCountEachPile tableOcc_; // replace tableOcc
    vector<AlphabetSymbol> alpha; //Corresponding between the alphabet, the piles and tableOcc
    AlphabetSymbol sizeAlpha;  //number of the different symbols in the input texts
    AlphabetSymbol *alphaInverse;  //Corresponding between alpha[i] and the symbol as char

    vector< vector< vector<LetterNumber> > > vectorOcc;
    vector <LetterNumber> numBlocksInPartialBWT;

    //LetterNumber*** vectorOcc;
    vector<sortElement> FirstVector, LastVector;

    CompressionFormatType outputCompression_;


    /**
    * Init an instance of a text collection object
    *
    * Returns a pointer to an object implementing this interface.
    */
    static BWTCollection *InitBWTCollection
    ( const string &file1, const string &fileOut, const int mode,
      const CompressionFormatType outputCompression );

    /**
     * Virtual destructor
     */
    virtual ~BWTCollection() {}
    /**
     *
     * The i'th text insertion gets an identifier value i-1.
     * In other words, document identifiers start from 0.
     */
    virtual int buildBCR( const string &, const string &, const BwtParameters *bwtParams ) = 0;
    virtual int unbuildBCR( char const *, char const *, char const *, char const * ) = 0;
    virtual int backwardSearchBCR( char const * , char const * , char const * , char const * ) = 0;
    virtual int decodeBCRnaiveForward( char const *, char const *, char const * ) = 0; //Inverse BWT by Forward direction of nText sequences, one sequence at a time, in lexicographic order.
    virtual int decodeBCRmultipleReverse( char const *, char const *, char const *, bool processQualities = false ) = 0; //Inverse BWT by Backward direction of nText sequences at the same time by lengthRead iterations.
    virtual int RecoverNsymbolsReverse( char const *, char const *, uchar *, uchar *newQual = 0 ) = 0;
    virtual int RecoverNsymbolsReverseByVector( char const *file1, char const *fileOutBwt, uchar *newSymb, uchar *newQual = 0 ) = 0;
    virtual int Recover1symbolReverse( char const * , char const * , uchar *, sortElement * ) = 0;
    virtual SequenceNumber recover1SequenceForward( char const * , char const * , sortElement , uchar *, SequenceLength * ) = 0 ;
    virtual vector <int> recoverNSequenceForward( char const * , char const *, SequenceNumber ) = 0;
    virtual int recoverNSequenceForwardSequentially( char const * , char const *, SequenceNumber ) = 0;
    virtual void storeBWT( uchar const *, uchar const *qual = NULL ) = 0;
    virtual void storeEntireBWT( const string & ) = 0;
    virtual void storeSA( SequenceLength ) = 0;
    virtual void storeEntirePairSA( const char * ) = 0;
    virtual void storeEntireSAfromPairSA( const char * ) = 0;
    virtual void storeBWTandLCP( uchar const * ) = 0;
    virtual void storeEntireLCP( const string & ) = 0;
    virtual LetterNumber rankManySymbols( FILE &, LetterNumber *, LetterNumber, uchar * ) = 0;
    virtual LetterNumber rankManySymbolsByVector( FILE & , LetterNumber *, LetterNumber, uchar *, uchar *foundQual = 0, FILE *InFileBWTQual = 0 ) = 0;
    virtual LetterNumber findRankInBWT ( char const *, char const *, AlphabetSymbol, LetterNumber, uchar ) = 0;
    virtual LetterNumber findRankInBWTbyVector ( char const *, char const *, AlphabetSymbol, LetterNumber, uchar ) = 0;
    virtual int rankInverseManyByVector ( char const * , char const * , SequenceNumber , uchar * ) = 0;
    virtual int backwardSearchManyBCR( char const * , char const *, char const *, vector<string>, SequenceLength ) = 0;
    virtual int SearchAndLocateKmer ( char const * , char const * , char const * , vector<string> , SequenceLength, vector <int> & ) = 0;
private:

    virtual void InsertNsymbols( uchar const *, SequenceLength, uchar const *qual = NULL ) = 0;
    virtual void InsertFirstsymbols( uchar const *, uchar const *qual = NULL, const int subSequenceNum = 0 ) = 0;
    virtual int initializeUnbuildBCR( char const *, char const *, LetterNumber [] ) = 0;
    virtual int computeNewPositionForBackSearch ( char const *, char const *, uchar ) = 0;
    virtual int computeNewPositionForBackSearchByVector ( char const *, char const *, uchar ) = 0;
    virtual int computeVectorUnbuildBCR( char const *, char const *, LetterNumber [] ) = 0;
    virtual int computeManyNewPositionForBackSearchByVector( char const * , char const * , uchar *, SequenceNumber ) = 0;
    virtual int update_Pos_Pile( sortElement * ) = 0;
    virtual int update_Pos_Pile_Blocks( LetterNumber *, LetterNumber *, AlphabetSymbol, uchar ) = 0;
    virtual int findBlockToRead( LetterNumber *, AlphabetSymbol , LetterNumber *, LetterNumber * ) = 0;
protected:
    // Protected constructor; call the static function InitBWTCollection().
    BWTCollection() : tableOcc( NULL ), alpha( 256 ), sizeAlpha( 0 ), alphaInverse( NULL ), outputCompression_( compressionASCII ) { }

    // No copy constructor or assignment
    BWTCollection( BWTCollection const & );
    BWTCollection &operator = ( BWTCollection const & );
};
}
#endif
