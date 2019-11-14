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

#ifndef FILENAME_HH
#define FILENAME_HH

#include <iostream>
#include <string>
#include <sstream>

using std::string;


typedef const char *constCharStar;


class Filename
{
public:
    Filename( const std::string &str )
        : str_( str )
    {
    }

    Filename( const int i )
    {
        std::ostringstream fn;
        fn << i;
        str_ = fn.str();
    }

    Filename( const std::string &part1, const std::string &part2 )
    {
        str_ = part1 + part2;
    }

    Filename( const char *part1, const std::string &part2 )
    {
        std::ostringstream fn;
        fn << part1 << part2;
        str_ = fn.str();
    }

    Filename( const std::string &part1, const int part2, const std::string &part3 = "" )
    {
        std::ostringstream fn;
        fn << part1 << part2 << part3;
        str_ = fn.str();
    }

    Filename( const char *part1, const std::string &part2, const int part3, const std::string &part4 = "" )
    {
        std::ostringstream fn;
        fn << part1 << part2 << part3 << part4;
        str_ = fn.str();
    }

    Filename( const std::string &part1, const std::string &part2, const int part3, const std::string &part4 = "" )
    {
        std::ostringstream fn;
        fn << part1 << part2 << part3 << part4;
        str_ = fn.str();
    }

    Filename( const std::string &part1, const int part2, const std::string &part3, const int part4, const std::string &part5 = "" )
    {
        std::ostringstream fn;
        fn << part1 << part2 << part3 << part4 << part5;
        str_ = fn.str();
    }

    // automatic cast operator to string
    operator std::string( void ) const
    {
        return str();
    }

    // automatic cast operator to const char*
    operator constCharStar( void )
    {
        str_ = str();
        return str_.c_str();
    }

    virtual ~Filename() {}
    virtual std::string str() const;

protected:
    string str_;
};


class TmpFilename : public Filename
{
public:
    TmpFilename( const std::string &str ) : Filename( str ) {}
    TmpFilename( const int i ) : Filename( i ) {}
    TmpFilename( const std::string &part1, const std::string &part2 ) : Filename( part1, part2 ) {}
    TmpFilename( const char *part1, const std::string &part2 ) : Filename( part1, part2 ) {}
    TmpFilename( const std::string &part1, const std::string &part2, const int part3 ) : Filename( part1, part2, part3 ) {}
    TmpFilename( const std::string &part1, const int part2, const std::string &part3 = "" ) : Filename( part1, part2, part3 ) {}
    TmpFilename( const char *part1, const std::string &part2, const int part3, const std::string &part4 = "" ) : Filename( part1, part2, part3, part4 ) {}
    TmpFilename( const std::string &part1, const int part2, const std::string &part3, const int part4, const std::string &part5 = "" ) : Filename( part1, part2, part3, part4, part5 ) {}

    virtual ~TmpFilename() {}
    virtual std::string str() const;

private:
    mutable string fullPathStr_;
};


#endif // FILENAME_HH
