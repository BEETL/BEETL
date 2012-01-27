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

#ifndef DEFINED_LETTERCOUNT_HH
#define DEFINED_LETTERCOUNT_HH

#include <vector>
#include <iostream>
#include <assert.h>
#include "Types.hh"
#include "Alphabet.hh"

struct LetterCount
{
  LetterCount() 
  { 
    clear();
  } // ~ctor
  void clear( void )
  { 
    for (int i(0);i<alphabetSize;i++) count_[i]=0; 
  } // ~clear

  void print( void )
  {
    for (int i(0);i<alphabetSize;i++) 
      std::cout << " " << alphabet[i] << ":" << count_[i];
    std::cout << std::endl;
  } // ~print

  void operator+=( char c)
  {
    assert(whichPile[(int)c]<alphabetSize);
    count_[whichPile[(int)c]]++;
  }

  void operator+=( const LetterCount& rhs )
  { 
    for (int i(0);i<alphabetSize;i++) count_[i]+=rhs.count_[i]; 
  } // ~clear

  //  LetterCountData count_;
  LetterCountType count_[alphabetSize];
}; // ~LetterCount

struct LetterCountEachPile : public std::vector<LetterCount>
{
  LetterCountEachPile()
  {
    resize(alphabetSize); clear();
  }  // ~ctor
  void clear( void )
  { 
    for (int i(0);i<alphabetSize;i++) (*this)[i].clear(); 
  } // ~clear

  void print( void )
  { 
    for (int i(0);i<alphabetSize;i++) 
    {
      std::cout << alphabet[i] << " pile";
      (*this)[i].print(); 
    } // ~for
  } // ~clear

  void operator+=( const LetterCountEachPile& rhs )
  { 
    for (int i(0);i<alphabetSize;i++) (*this)[i]+=rhs[i]; 
  } // ~clear


};



#endif
