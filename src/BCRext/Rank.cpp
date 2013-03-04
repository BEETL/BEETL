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

// Encapsulate the functionality for rank and inverse-rank queries;
// operates on a single BWT file (agnostic to whether it's a partial
// or complete BWT)

#include "Rank.hh"


Rank::Rank( void )
{
    charCount_ = 0;

    for ( int i=0; i<255; i++ )
    {
        StaticOccComparison tmp( i );
        compareObjects_.push_back( tmp );
    }
}


bool Rank::initialize( const string &fn, const int &blockSize )
{
    assert( ( alphaInitialized_==true ) && ( alphaInverseInitialized_==true ) );
    assert( charCount_ != 0 );

    // open the BWT
    ifsBWT = fopen( fn.c_str(),"r" );

    // get length of file:
    fseek( ifsBWT,0,SEEK_END ); // TODO, move to the end!!
    off_t bwtLength = ftell( ifsBWT );
    fseek( ifsBWT,0,SEEK_SET ); // TODO, move to the end!!


    cerr << "Reading " << bwtLength << " characters from " << fn << "." << endl;
    cerr << "Blocksize is set to " << blockSize << endl;

    blockSize_ = blockSize;

    dataTypeNChar noOfBlocks = ( dataTypeNChar )bwtLength/blockSize + 1; // + 1 to hold the fraction of bases
    // resize occ_
    cerr << "Resizing occ_ to hold " << noOfBlocks << " blocks of size " << charCount_ << "." << endl;
    cerr << "character count = " << charCount_ << endl;
    occ_.resize( noOfBlocks );
    for ( dataTypeNChar i=0; i<noOfBlocks; i++ )
    {
        occ_[i].resize( charCount_ );
    }



    char *buffer = new char[blockSize];
    for ( int i=0; i<blockSize; i++ )
    {
        buffer[i]='\0';
    }

    int totalCharRead = 0;
    int blockIdx = 0;

    // parse the current block
    vector<dataTypeNSeq> curBlock( charCount_ );

    while ( !feof( ifsBWT ) )
    {
        int charsRead = fread( buffer,sizeof( char ),blockSize,ifsBWT );

        assert( ( charsRead>=0 ) && ( charsRead<=blockSize ) );

        for ( int i=0; i<charsRead; i++ )
        {

#ifdef DEBUG
            cerr << "buffer[" << i << "] = " << buffer[i] << "\t"
                 << "(int)alpha_[buffer[i]] = " << ( int )alpha_[buffer[i]] << "\t";

            cerr << "(int)buffer[i] = " << ( int )buffer[i] << "\t" << "(int)alpha_[(int)buffer[i]] = "
                 << ( int )alpha_[( int )buffer[i]] << "\t" << "totalRead = " << totalCharRead << endl;
#endif

            int curMapping = ( int )alpha_[( int )buffer[i]];
            assert( ( curMapping>=0 )&&( curMapping<charCount_ ) );

            curBlock[ curMapping ]++;
            totalCharRead++;
        }



        occ_[blockIdx]=curBlock;

        // increment current block idx
        blockIdx++;

    }


    delete[] buffer;
    return true;
}


Rank::~Rank()
{


}


bool Rank::setAlpha( dataTypedimAlpha *alpha,const int &s )
{
    // resize internally
    alpha_.resize( s,0 );

    for ( int i=0; i<s; i++ )
    {
        alpha_[i]=alpha[i];
    }


    alphaInitialized_=true;
    return true;
}


bool Rank::setAlphaInverse( dataTypedimAlpha *alphaInverse,const int &s )
{
    // resize internally
    alphaInverse_.resize( s,0 );

    for ( int i=0; i<s; i++ )
    {
        alphaInverse_[i]=alphaInverse[i];
    }


    alphaInverseInitialized_=true;
    return true;
}



dataTypeNChar Rank::getInverseRank( const char &c,const dataTypeNChar &n )
{

    // use STL lower_bound to get the first element that does not
    // violate the ordering
    int mappedIdx = ( int )alpha_[( int )c];

    //  vector< vector<dataTypeNSeq> >::iterator it = lower_bound( occ_.begin(),occ_.end(),n,OccComparison::CompareOcc );
    vector< vector<dataTypeNSeq> >::iterator it = lower_bound( occ_.begin(),occ_.end(),n,compareObjects_[mappedIdx] );

    // calculate position within the block
    dataTypeNChar blockIdx = it-occ_.begin();

    // seek to position blockIdx*blockSize_
    dataTypeNChar bwtPos = blockIdx*blockSize_;
    dataTypeNChar characterCount = ( *it )[mappedIdx];

    fseek( ifsBWT,bwtPos,SEEK_SET ); // TODO, move to the end!!

    char *buffer = new char[blockSize_];
    int charsRead = fread( buffer,sizeof( char ),blockSize_,ifsBWT );

    for ( int i=0; i<charsRead; i++ )
    {
        if ( buffer[i]==c )
        {
            characterCount++;
        }
        bwtPos++; // moving on, increment current position

        if ( characterCount==n )
        {
            break; // stop
        }
    }



    return bwtPos;
}


