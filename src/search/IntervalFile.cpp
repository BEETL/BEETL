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

#include "IntervalFile.hh"

#include <cassert>
#include <iostream>

using namespace std;


ostream &operator<<( std::ostream &os, const IntervalRecord &obj )
{
    return os << "{ " << obj.kmer << " interval, position: " << obj.position << ", count: " << obj.count << " }";
}

IntervalWriter::IntervalWriter( std::ostream &file ): file_( file )
{
    assert( file_.good() );
}

void IntervalWriter::write( const IntervalRecord &ir ) const
{
    file_ << ir.kmer;
    file_ << ' ' << ir.position;
    file_ << ' ' << ir.count << endl;
}

void IntervalWriter::writeV2( const IntervalRecord &ir ) const
{
    file_ << ir.kmer;
    file_ << ' ' << ir.position;
    file_ << ' ' << ir.count << ':';
    for ( auto pos : ir.dollarSignPositions )
        file_ << ' ' << pos;
    file_ << endl;
}

bool IntervalReader::read( IntervalRecord &ir )
{
    file_ >> ir.kmer;
    file_ >> ir.position;
    file_ >> ir.count;
    if ( file_.good() )
        return true;
    else
        return false;
}

vector<IntervalRecord> IntervalReader::readFullFileAsVector()
{
    vector<IntervalRecord> result;
    IntervalRecord rec;

    while ( read( rec ) )
        result.push_back( rec );

    return result;
}
