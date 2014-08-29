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

#ifndef DEFINED_LETTERCOUNT_HH
#define DEFINED_LETTERCOUNT_HH

#include "Alphabet.hh"
#include "Types.hh"

#include <cassert>
#include <iostream>
#include <vector>

using std::vector;


template<typename T> struct LetterCountTemplate
{
    LetterCountTemplate()
    {
        clear();
    } // ~ctor
    void clear( void )
    {
        for ( int i( 0 ); i < alphabetSize; i++ ) count_[i] = 0;
    } // ~clear

    void operator+=( const char c )
    {
        assert( whichPile[( int )c] < alphabetSize );
        count_[whichPile[( int )c]]++;
        assert( count_[whichPile[( int )c]] != 0 );
    }

    template<typename TT>void operator+=( const LetterCountTemplate<TT> &rhs )
    {
        for ( int i( 0 ); i < alphabetSize; i++ )
        {
            T newCount = count_[i] + rhs.count_[i];
            assert( newCount >= count_[i] && "Overflow error in LetterCountTemplate" ); 
            count_[i] = newCount;
        }
    } // ~clear

    template<typename TT>void operator-=( const LetterCountTemplate<TT> &rhs )
    {
        // on your own head be it if you make an unsigned quantity negative...
        for ( int i( 0 ); i < alphabetSize; i++ ) count_[i] -= rhs.count_[i];
    } // ~clear

    void countString( const char *const s, const T length )
    {
        for ( T i = 0; i < length; ++i )
            operator+=( s[i] );
    }

    friend std::ostream &operator<<( std::ostream &os, const LetterCountTemplate &obj )
    {
        for ( int i( 0 ); i < alphabetSize; i++ )
            os << " " << alphabet[i] << ":" << obj.count_[i];
        return os;
    }

    friend std::istream &operator>>( std::istream &is, LetterCountTemplate &obj )
    {
        char letter, colon;
        for ( int i( 0 ); i < alphabetSize; i++ )
        {
            is >> letter >> colon >> obj.count_[i];
            if ( letter != alphabet[i] || colon != ':' )
            {
                std::cerr << "Error reading LetterCount from file: letter=" << letter << ", colon=" << colon << std::endl;
                exit( 1 );
            }
        }
        return is;
    }

    //  LetterCountData count_;
    T count_[alphabetSize];
}; // ~LetterCountTemplate

typedef LetterCountTemplate<LetterNumber> LetterCount;
typedef LetterCountTemplate<LetterNumberCompact> LetterCountCompact;

struct LetterCountEachPile : public vector<LetterCount>
{
    LetterCountEachPile()
    {
        resize( alphabetSize );
        clear();
    }  // ~ctor
    void clear( const int startIndex = 0 )
    {
        for ( int i( startIndex ); i < alphabetSize; i++ ) ( *this )[i].clear();
    } // ~clear

    void print( void )
    {
        for ( int i( 0 ); i < alphabetSize; i++ )
        {
            std::cout << alphabet[i] << " pile" << ( *this )[i] << std::endl;
        } // ~for
    } // ~clear

    void operator+=( const LetterCountEachPile &rhs )
    {
        for ( int i( 0 ); i < alphabetSize; i++ ) ( *this )[i] += rhs[i];
    } // ~clear

};


#endif
