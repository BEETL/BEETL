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

#ifndef INCLUDED_RANGE_HH
#define INCLUDED_RANGE_HH

#include "Config.hh"
#include "Alphabet.hh"
#include "LetterCount.hh"
#include "Tools.hh"
#include "libzoo/util/TemporaryFilesManager.hh"

#include <string>


#define ENCODE_POSITIONS_AS_OFFSETS


// The following helper methods should be moved to a stream class
void writeCompressedNum( TemporaryFile *pFile, LetterNumber num );
void writeBytes( TemporaryFile *pFile, bool val, int byteCount );
bool readCompressedNum( TemporaryFile *pFile, LetterNumber &num );
void readBytes( TemporaryFile *pFile, bool &val, const int byteCount );


class RangeState;

class Range
{
public:
    Range(
        const string &word,
        const LetterNumber pos,
        const LetterNumber num,
        const bool isBkptExtension = false
    ) :
#ifdef PROPAGATE_PREFIX
        word_( word ),
#endif
        pos_( pos ),
        num_( num ),
        isBkptExtension_( isBkptExtension )
    {}

    Range( void ) :
        pos_( 0 ),
        num_( 0 ),
        isBkptExtension_( false )
    {}

    virtual ~Range() {}

    virtual void clear()
    {
        *this = Range();
    }

    virtual bool writeTo( TemporaryFile *pFile, RangeState &currentState ) const;
    virtual bool readFrom( TemporaryFile *pFile, RangeState &currentState );

#ifdef PROPAGATE_PREFIX
    string word_;
#endif
    LetterNumber pos_;
    LetterNumber num_;
    bool isBkptExtension_;
};


//
// RangeState: this is a proxy class that sits in RangeStoreExternal
// and handles the optional storage of compressed sequence
// when PROPAGATE_PREFIX is set
//
struct RangeState
{
    RangeState();
    void clear( void );

    TemporaryFile *pFile_;
#ifdef ENCODE_POSITIONS_AS_OFFSETS
    LetterNumber lastProcessedPos_;
#endif

    RangeState &operator<<( const Range & );
    RangeState &operator>>( Range & );
    bool good();

#ifdef PROPAGATE_PREFIX
private:
    void addSeq( const string &seq );
    void getSeq( string &word );

    char wordLast_[256];
#endif //ifdef PROPAGATE_PREFIX
};


#endif
