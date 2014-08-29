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

#include "EndPosFile.hh"

#include <cassert>
#include <fstream>

using namespace std;


EndPosFile::EndPosFile( const string &bwtFilenamePrefix )
    : file_( bwtFilenamePrefix + "-end-pos" )
    , sequenceGroupCount_( 0 )
    , sequenceCountInGroup_( 0 )
    , hasRevComp_( 0 )
{
    file_.read( reinterpret_cast< char * >( &sequenceGroupCount_ ), sizeof( SequenceNumber ) );
    file_.read( reinterpret_cast< char * >( &sequenceCountInGroup_ ), sizeof( uint8_t ) );
    file_.read( reinterpret_cast< char * >( &hasRevComp_ ), sizeof( uint8_t ) );
    //    assert( file_.good() );
    dollarSignCount_ = sequenceGroupCount_ * sequenceCountInGroup_ * ( hasRevComp_ ? 2 : 1 );


}

SequenceNumber EndPosFile::convertDollarNumToSequenceNum( const SequenceNumber dollarNum )
//SequenceNumber EndPosFile_convertDollarNumToSequenceNum( const SequenceNumber dollarNum )
{
    assert( file_.good() && "Error: -end-pos file not readable" );

    assert( dollarNum < dollarSignCount_ );
    /*
      if ( dollarPos >= numDollarEntries )
      {
      cout << "Warning: dollarPos " << dollarPos << " >= numDollarEntries " << numDollarEntries << endl;
      //                    continue;
      dollarPos %= numDollarEntries;
      }
    */
    file_.seekg( sizeof( SequenceNumber ) + 2 * sizeof( uint8_t ) + ( dollarNum ) * ( sizeof( SequenceNumber ) + sizeof( uint8_t ) ) );

    SequenceNumber sequenceGroupNum;
    file_.read( reinterpret_cast< char * >( &sequenceGroupNum ), sizeof( SequenceNumber ) );
    uint8_t positionInGroup;
    file_.read( reinterpret_cast< char * >( &positionInGroup ), sizeof( uint8_t ) );

    SequenceNumber sequenceNum = sequenceGroupNum + positionInGroup * sequenceGroupCount_;

    return sequenceNum;
}
