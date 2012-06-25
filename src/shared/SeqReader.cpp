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


#include <cstring>
#include <cstdlib>
#include <math.h>
#include <iostream>

//#include "BwtReader.hh"
//#include "LetterCount.hh"
//#include "BwtWriter.hh"
//#include "CountWords.hh"
#include "SeqReader.hh"

//#define DEBUG 1



using namespace std;





//
// SeqReaderBase member function definitions
//



SeqReaderBase::SeqReaderBase() {}

SeqReaderBase::~SeqReaderBase() {}


//
// SeqReaderFile member function definitions
//

SeqReaderFile::SeqReaderFile( FILE* pFile ) : 
  pFile_(pFile), allRead_(false), length_(-1) {}

SeqReaderFile::~SeqReaderFile() {}


SeqReaderFile* SeqReaderFile::getReader( FILE* pFile )
{
  int i(fgetc(pFile));
  char c((char)i); // TBD check for error condition
  ungetc(i,pFile);
  if (c=='>')
  {
    return new SeqReaderFasta(pFile);
  }
  else if (c=='@')
  {
    return new SeqReaderFastq(pFile);
  }
  else if (whichPile[i]!=nv)
  {
    return new SeqReaderRaw(pFile);
  }
  else
  {
    cerr << "Unable to deduce file type from first char (char code = "
	 << i << " )" << endl;
    exit (-1);
  } // ~else
} // ~getReader

const char* SeqReaderFile::thisSeq( void )
{ return bufSeq_; }
const char* SeqReaderFile::thisQual( void )
{ return bufQual_; }
const char* SeqReaderFile::thisName( void )
{ return bufName_; }
bool SeqReaderFile::allRead( void ) const
{ return allRead_; }
int SeqReaderFile::length( void ) const
{ return length_; }



//
// SeqReaderRaw member function definitions
//

SeqReaderRaw::SeqReaderRaw( FILE* pFile ) : SeqReaderFile(pFile) 
{
  cout << "Creating SeqReaderRaw" << endl;
  readNext();
  if (allRead()==true)
  {
    cerr << "No sequences in file!" << endl;
    exit(-1);
  }
  else
  {
    length_=strlen(bufSeq_)-1;
    cerr << "Deducing read length of " << length_ << endl;
  }
}

SeqReaderRaw::~SeqReaderRaw() {}


void SeqReaderRaw::readNext(void)
{
  //  cout << "readNext" << endl;
  if (allRead_==true)
  {
    cerr << "Tried to read an empty sequence stream" << endl;
    exit (-1);
  }
  else if( fgets( bufSeq_, maxSeqSize, pFile_)==NULL)
  {
    allRead_=true;
  }
  else if ((length_!=-1)&&(((int)strlen(bufSeq_))!=length_+1))
  {
    cerr << "Length of current sequence does not match length of first" << endl;    exit (-1);

  }
}


//
// SeqReaderFasta member function definitions
//

SeqReaderFasta::SeqReaderFasta( FILE* pFile ) : SeqReaderFile(pFile) 
{
  cout << "Creating SeqReaderFasta" << endl;
  readNext();
  if (allRead()==true)
  {
    cerr << "No sequences in file!" << endl;
    exit(-1);
  }
  else
  {
    length_=strlen(bufSeq_)-1;
    cerr << "Deducing read length of " << length_ << endl;
  }
}

SeqReaderFasta::~SeqReaderFasta() {}


void SeqReaderFasta::readNext(void)
{
  if (allRead_==true)
  {
    cerr << "Tried to read an empty sequence stream" << endl;
    exit (-1);
  }
  else if( fgets( bufName_, maxSeqSize, pFile_)==NULL)
  {
    allRead_=true;
  }
  else
  {
    if (bufName_[0]!='>')
    {
      cerr << "Expected FASTA header, got " << bufName_ << endl;
      exit (-1);
    }
    if( fgets( bufSeq_, maxSeqSize, pFile_)==NULL)
    {
      cerr << "read FASTA header with no entry, incomplete file?" << endl;
      exit (-1);
    }
    else if ((length_!=-1)&&(((int)strlen(bufSeq_))!=length_+1))
      //else if (strlen(bufSeq_)!=length_)
    {
      cerr << "Length of current sequence does not match length of first" << endl;    exit (-1);      
    }

  }

}




//
// SeqReaderFastq member function definitions
//

SeqReaderFastq::SeqReaderFastq( FILE* pFile ) : SeqReaderFile(pFile) 
{
  cout << "Creating SeqReaderFastq" << endl;
  readNext();
  if (allRead()==true)
  {
    cerr << "No sequences in file!" << endl;
    exit(-1);
  }
  else
  {
    length_=strlen(bufSeq_)-1;
    cerr << "Deducing read length of " << length_ << endl;
  }
}

SeqReaderFastq::~SeqReaderFastq() {}


void SeqReaderFastq::readNext(void)
{  
  if (allRead_==true)
  {
    cerr << "Tried to read an empty sequence stream" << endl;
    exit (-1);
  }
  else if( fgets( bufName_, maxSeqSize, pFile_)==NULL)
  {
    allRead_=true;
  }
  else
  {
    if (bufName_[0]!='@')
    {
      cerr << "Expected FASTQ header, got " << bufName_ << endl;
      exit (-1);
    }
    if( fgets( bufSeq_, maxSeqSize, pFile_)==NULL)
    {
      cerr << "read FASTA header with no entry, incomplete file?" << endl;
      exit (-1);
    }
    else 
    {
      if ((length_!=-1)&&(((int)strlen(bufSeq_))!=length_+1))
	//else if (strlen(bufSeq_)!=length_)
      {
	cerr << "Length of current sequence does not match length of first" << endl;    exit (-1);      
      }
      else if ( fgets( bufQual_, maxSeqSize, pFile_)==NULL)
      {
	cerr << "Could not read FASTQ quality spacer, incomplete file?" << endl;
	exit (-1);
      }
      else if (bufQual_[0]!='+')
      {
	cerr << "Expected FASTQ quality spacer, got " << bufQual_ << endl;
	exit (-1);
      }
      else if ( fgets( bufQual_, maxSeqSize, pFile_)==NULL)
      {
	cerr << "Could not read FASTQ quality string, incomplete file?" << endl;
	exit (-1);
      }
    } // ~else
  } // ~else
} // ~Fastq::readNext



