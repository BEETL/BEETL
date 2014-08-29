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

#ifndef FRAGMENTED_VECTOR_HH
#define FRAGMENTED_VECTOR_HH

#include <vector>

using std::vector;


template <class T>
class FragmentedVector
{
public:
    FragmentedVector()
    {
        chunks_.resize( 1 );
        chunks_[0].reserve( 1024 );
        assert( chunks_[0].capacity() == 1024 );
    }

    FragmentedVector( const FragmentedVector<T> &obj )
    {
        // Warning: this copy constructor only copies empty objects
        assert( obj.empty() );
        chunks_.resize( 1 );
        chunks_[0].reserve( 1024 );
        assert( chunks_[0].capacity() == 1024 );
    }

    bool empty() const
    {
        //        assert( valid.empty() == chunks_[0].empty() );
        return chunks_[0].empty();
    }

    void push_back( const T &obj )
    {
        size_t lastSize = chunks_.back().size();
        if ( lastSize == chunks_.back().capacity() )
        {
            chunks_.resize( chunks_.size() + 1 );
            chunks_.back().reserve( min<size_t>( lastSize * 2, 1024 * 1024 * 1024 / sizeof( T ) ) ); // max growing chunk of 1GB
        }
        chunks_.back().push_back( obj );
    }

    void appendTo( vector<T> &dest )
    {
        for ( unsigned int i = 0; i < chunks_.size(); ++i )
            dest.insert( dest.end(), chunks_[i].begin(), chunks_[i].end() );
    }

    /*
    // If this method is to be used often, implement it using a size_t totalSize member;
        size_t size() const
        {
            size_t result = 0;
            for (unsigned int i=0; i<chunks_.size(); ++i)
                result += chunks_[i].size();
            return result;
        }
    */

private:
    vector< vector< T > > chunks_;
};

#endif // FRAGMENTED_VECTOR_HH
