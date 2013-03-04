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

#ifndef INCLUDED_SEQREADER_HH
#define INCLUDED_SEQREADER_HH

#include "Alphabet.hh"
#include "Config.hh"
#include "Types.hh"

#include <cassert>
#include <cstdio>
#include <string>


class SeqReaderBase
{
public:
    SeqReaderBase();
    virtual ~SeqReaderBase();
    virtual void readNext( char *seqBuf = NULL )=0;
    virtual const char *thisSeq( void )=0;
    virtual const char *thisQual( void )=0;
    virtual const char *thisName( void )=0;
    virtual bool allRead( void ) const=0;
    virtual int length( void ) const=0;
}; // ~class BwtReaderBase

class SeqReaderFile : public SeqReaderBase
{
public:
    static SeqReaderFile *getReader( FILE *pFile );


    SeqReaderFile( FILE *pFile );
    virtual ~SeqReaderFile();

    virtual void readNext( char *seqBuf = NULL )=0;
    virtual const char *thisSeq( void );
    virtual const char *thisQual( void );
    virtual const char *thisName( void );
    virtual bool allRead( void ) const;
    virtual int length( void ) const;
protected:
    FILE *pFile_;
    char bufSeq_[1+maxSeqSize];
    char bufQual_[1+maxSeqSize];
    char bufName_[1+maxSeqSize];
    bool allRead_;
    int length_;
}; // ~class BwtReaderBase

class SeqReaderRaw: public SeqReaderFile
{
public:
    SeqReaderRaw( FILE *pFile );
    virtual ~SeqReaderRaw();

    virtual void readNext( char *seqBuf = NULL );
    //  virtual const char* thisSeq( void );
    //  virtual const char* thisQual( void );
    // virtual const char* thisName( void );
    //  virtual bool allRead( void ) const=0;

};

class SeqReaderFasta: public SeqReaderFile
{
public:
    SeqReaderFasta( FILE *pFile );
    virtual ~SeqReaderFasta();

    virtual void readNext( char *seqBuf = NULL );
    //  virtual const char* thisSeq( void );
    //  virtual const char* thisQual( void );
    // virtual const char* thisName( void );
    //  virtual bool allRead( void ) const=0;
};

class SeqReaderFastq: public SeqReaderFile
{
public:
    SeqReaderFastq( FILE *pFile );
    virtual ~SeqReaderFastq();

    virtual void readNext( char *seqBuf = NULL );
    //  virtual const char* thisSeq( void );
    //  virtual const char* thisQual( void );
    //  virtual const char* thisName( void );
    // virtual bool allRead( void ) const=0;
};




#endif
