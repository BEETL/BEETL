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

#ifndef BWT_INDEX_HH
#define BWT_INDEX_HH

#include "BwtReader.hh"
#include "BwtWriter.hh"

#include <string>
#include <vector>

using std::string;
using std::vector;

//#define USE_COMPACT_STRUCTURES
#ifdef USE_COMPACT_STRUCTURES
# define LETTER_COUNT_CLASS LetterCountCompact
#else
# define LETTER_COUNT_CLASS LetterCount
#endif


const vector<char> indexV1Header = { 'B', 'W', 'I', 13, 10, 26, 1, 0 };
const vector<char> indexV2Header = { 'B', 'W', 'I', 13, 10, 26, 2, 0 };

template< class T >
class BwtReaderIndex : public T
{
public:
    BwtReaderIndex( const string &filename, const string &optionalSharedMemoryPath );
    //    BwtReaderIndex( const BwtReaderIndex & );

    BwtReaderIndex( const BwtReaderIndex &obj ) :
        T( obj ),
        indexFilename_( obj.indexFilename_ ),
        pIndexFile_( obj.pIndexFile_ ),
        indexPosBwt_( obj.indexPosBwt_ ),
        indexPosFile_( obj.indexPosFile_ ),
        indexCount_( obj.indexCount_ ),
        indexNext_( obj.indexNext_ )
    {
        assert( pIndexFile_ == NULL ); // If it's not NULL, we may try to fclose it multiple times
    }


    virtual ~BwtReaderIndex() {}
    virtual BwtReaderIndex *clone() const
    {
        return new BwtReaderIndex( *this );
    };

    virtual LetterNumber readAndCount( LetterCount &c, const LetterNumber numChars );

    virtual LetterNumber readAndSend( BwtWriterBase &writer, const LetterNumber numChars )
    {
        assert( 1 == 0 );
    }

    virtual LetterNumber operator()( char *p, LetterNumber numChars )
    {
        return T::operator()( p, numChars );
        //        assert( 1 == 0 );
    }

    virtual void rewindFile( void );

    //  virtual LetterNumber tellg( void ) const;

    void initIndex( const string &optionalSharedMemoryPath );

    //  bool getRun(void);
protected:


    string indexFilename_;


    FILE *pIndexFile_;

    vector<LetterNumber> indexPosBwt0_;
    vector<LetterNumber> indexPosFile0_;
    vector<LETTER_COUNT_CLASS> indexCount0_;

    // Pointers to the same structure, used in case of mmapped files
    LetterNumber *indexPosBwt_;
    LetterNumber *indexPosFile_;
    LETTER_COUNT_CLASS *indexCount_;
    uint32_t indexSize_;

    uint32_t indexNext_;
};


void buildIndex( BwtReaderBase *reader, FILE *pFile, const int indexBinSize );


#endif //ifdef BWT_INDEX_HH
