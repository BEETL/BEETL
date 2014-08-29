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

#ifndef SEARCH_USING_BACKTRACKER_HH
#define SEARCH_USING_BACKTRACKER_HH

#include "Algorithm.hh"

#include <string>

using std::string;
class SearchParameters;


class SearchUsingBacktracker : public Algorithm
{

public:
    SearchUsingBacktracker(
        const SearchParameters &searchParams
    );

    virtual ~SearchUsingBacktracker() {}
    virtual void run( void );

private:
    const SearchParameters &searchParams_;
};

#endif // SEARCH_USING_BACKTRACKER_HH
