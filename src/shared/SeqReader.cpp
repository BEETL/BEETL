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

#include "SeqReader.hh"

#include "Logger.hh"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace std;


//#define DEBUG 1

//
// SeqReaderBase member function definitions
//

SeqReaderBase::SeqReaderBase() {}

SeqReaderBase::~SeqReaderBase() {}


//
// SeqReaderFile member function definitions
//

SeqReaderFile::SeqReaderFile( FILE *pFile ) :
    pFile_( pFile ), allRead_( false ), length_( -1 )
{
    bufSeq_[0]  = 0;
    bufQual_[0] = 0;
    bufName_[0] = 0;
}

SeqReaderFile::~SeqReaderFile() {}


SeqReaderFile *SeqReaderFile::getReader( FILE *pFile )
{
    int i( fgetc( pFile ) );
    char c( ( char )i ); // TBD check for error condition
    ungetc( i,pFile );
    if ( c=='>' )
    {
        return new SeqReaderFasta( pFile );
    }
    else if ( c=='@' )
    {
        return new SeqReaderFastq( pFile );
    }
    else if ( whichPile[i]!=nv )
    {
        return new SeqReaderRaw( pFile );
    }
    else
    {
        Logger::error() << "Unable to deduce file type from first char (char code = "
                        << i << " )" << endl;
        exit( EXIT_FAILURE );
    } // ~else
} // ~getReader

const char *SeqReaderFile::thisSeq( void )
{
    return bufSeq_;
}
const char *SeqReaderFile::thisQual( void )
{
    return bufQual_;
}
const char *SeqReaderFile::thisName( void )
{
    return bufName_;
}
bool SeqReaderFile::allRead( void ) const
{
    return allRead_;
}
int SeqReaderFile::length( void ) const
{
    return length_;
}



//
// SeqReaderRaw member function definitions
//

SeqReaderRaw::SeqReaderRaw( FILE *pFile ) : SeqReaderFile( pFile )
{
    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "Creating SeqReaderRaw" << endl;
    readNext();
    if ( allRead()==true )
    {
        Logger::error() << "No sequences in file!" << endl;
        exit( EXIT_FAILURE );
    }
    else
    {
        length_=strlen( bufSeq_ )-1;
        Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "Deducing read length of " << length_ << endl;
    }
}

SeqReaderRaw::~SeqReaderRaw() {}


void SeqReaderRaw::readNext( char *seqBuf )
{
    //  cout << "readNext" << endl;
    if ( allRead_==true )
    {
        Logger::error() << "Tried to read an empty sequence stream" << endl;
        exit( EXIT_FAILURE );
    }
    else if ( fgets( seqBuf?:bufSeq_, maxSeqSize, pFile_ )==NULL )
    {
        allRead_=true;
    }
    else if ( ( length_!=-1 )&&( ( ( int )strlen( seqBuf?:bufSeq_ ) )!=length_+1 ) )
    {
        Logger::error() << "Length of current sequence does not match length of first" << endl;
        exit( EXIT_FAILURE );

    }
}


//
// SeqReaderFasta member function definitions
//

SeqReaderFasta::SeqReaderFasta( FILE *pFile ) : SeqReaderFile( pFile )
{
    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "Creating SeqReaderFasta" << endl;
    readNext();
    if ( allRead()==true )
    {
        Logger::error() << "No sequences in file!" << endl;
        exit( EXIT_FAILURE );
    }
    else
    {
        length_=strlen( bufSeq_ )-1;
        Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "Deducing read length of " << length_ << endl;
    }
}

SeqReaderFasta::~SeqReaderFasta() {}


void SeqReaderFasta::readNext( char *seqBuf )
{
    if ( allRead_==true )
    {
        Logger::error() << "Tried to read an empty sequence stream" << endl;
        exit( EXIT_FAILURE );
    }
    else if ( fgets( bufName_, maxSeqSize, pFile_ )==NULL )
    {
        allRead_=true;
    }
    else
    {
        if ( bufName_[0]!='>' )
        {
            Logger::error() << "Expected FASTA header, got " << bufName_ << endl;
            exit( EXIT_FAILURE );
        }
        if ( fgets( seqBuf?:bufSeq_, maxSeqSize, pFile_ )==NULL )
        {
            Logger::error() << "read FASTA header with no entry, incomplete file?" << endl;
            exit( EXIT_FAILURE );
        }
        else if ( ( length_!=-1 )&&( ( ( int )strlen( seqBuf?:bufSeq_ ) )!=length_+1 ) )
            //else if (strlen(seqBuf?:bufSeq_)!=length_)
        {
            Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "Length of current sequence does not match length of first" << endl;
            exit( EXIT_FAILURE );
        }

    }

}




//
// SeqReaderFastq member function definitions
//

SeqReaderFastq::SeqReaderFastq( FILE *pFile ) : SeqReaderFile( pFile )
{
    Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "Creating SeqReaderFastq" << endl;
    readNext();
    if ( allRead()==true )
    {
        Logger::error() << "No sequences in file!" << endl;
        exit( EXIT_FAILURE );
    }
    else
    {
        length_=strlen( bufSeq_ )-1;
        Logger::out( LOG_SHOW_IF_VERY_VERBOSE ) << "Deducing read length of " << length_ << endl;
    }
}

SeqReaderFastq::~SeqReaderFastq() {}


void SeqReaderFastq::readNext( char *seqBuf )
{
    if ( allRead_==true )
    {
        Logger::error() << "Tried to read an empty sequence stream" << endl;
        exit( EXIT_FAILURE );
    }
    else if ( fgets( bufName_, maxSeqSize, pFile_ )==NULL )
    {
        allRead_=true;
    }
    else
    {
        if ( bufName_[0]!='@' )
        {
            Logger::error() << "Expected FASTQ header, got " << bufName_ << endl;
            exit( EXIT_FAILURE );
        }
        if ( fgets( seqBuf?:bufSeq_, maxSeqSize, pFile_ )==NULL )
        {
            Logger::error() << "read FASTA header with no entry, incomplete file?" << endl;
            exit( EXIT_FAILURE );
        }
        else
        {
            if ( ( length_!=-1 )&&( ( ( int )strlen( seqBuf?:bufSeq_ ) )!=length_+1 ) )
                //else if (strlen(seqBuf?:bufSeq_)!=length_)
            {
                Logger::error() << "Length of current sequence does not match length of first" << endl;
                exit( EXIT_FAILURE );
            }
            else if ( fgets( bufQual_, maxSeqSize, pFile_ )==NULL )
            {
                Logger::error() << "Could not read FASTQ quality spacer, incomplete file?" << endl;
                exit( EXIT_FAILURE );
            }
            else if ( bufQual_[0]!='+' )
            {
                Logger::error() << "Expected FASTQ quality spacer, got " << bufQual_ << endl;
                exit( EXIT_FAILURE );
            }
            else if ( fgets( bufQual_, maxSeqSize, pFile_ )==NULL )
            {
                Logger::error() << "Could not read FASTQ quality string, incomplete file?" << endl;
                exit( EXIT_FAILURE );
            }
        } // ~else
    } // ~else
} // ~Fastq::readNext



