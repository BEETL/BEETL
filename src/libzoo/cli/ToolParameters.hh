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

#ifndef BEETL_TOOL_PARAMETERS_HH
#define BEETL_TOOL_PARAMETERS_HH

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using std::string;
using std::vector;
using namespace std;


class ToolParameters;

enum { TYPE_STRING = 1, TYPE_CHOICE = 2, TYPE_INT = 4, TYPE_SWITCH = 8 };
enum { REQUIRED = 16, OPTIONAL = 0, ENVIRONMENT = 32 };
enum { AUTOMATED = 64, NOT_AUTOMATED = 0 };
struct ParameterEntry
{
    string stringId;
    int numId;
    string longCmdLineName;
    string shortCmdLineName;
    string description;
    string defaultValue;
    int flags;
    const string *possibleValues;
    ToolParameters *parent_;

    string userValue;
    int parsedValue;

    ParameterEntry( string stringId, int numId, string longCmdLineName, string shortCmdLineName, string description, string defaultValue, int flags, const string *possibleValues, ToolParameters *parent )
        : stringId( stringId )
        , numId( numId )
        , longCmdLineName( longCmdLineName )
        , shortCmdLineName( shortCmdLineName )
        , description( description )
        , defaultValue( defaultValue )
        , flags( flags )
        , possibleValues( possibleValues )
        , parent_( parent )
        , userValue( "" )
        , parsedValue( -1 )
    {}

    bool operator==( const int rhs );
    bool operator==( const string &rhs );
    operator int () const;
    operator string () const;
    ParameterEntry &operator=( const int rhs );
    ParameterEntry &operator=( const string &rhs );
    bool isSet() const;
    void silentSet( const int val )
    {
        set( val, true );
    }

private:
    void set( const int val, const bool isSilent = false );
    void set( const string &valString );
};

ostream &operator<<( std::ostream &os, const ParameterEntry &obj );



static const int MULTIPLE_OPTIONS = 99;



// options: --color

enum ColorFormat
{
    COLOR_FORMAT_NEVER = 0,
    COLOR_FORMAT_ALWAYS = 1,
    COLOR_FORMAT_AUTO = 2
};

static const string colorLabels[] =
{
    "never",
    "always",
    "auto",
    "" // end marker
};


class ToolParameters
{
public:
    virtual ~ToolParameters() {}

    int getValue( const int key ) const;
    int getValue( const string &key ) const;

    string getStringValue( const int key ) const;
    string getStringValue( const string &key ) const;

    ParameterEntry &operator[]( const string &key )
    {
        return getEntry( key );
    }
    ParameterEntry operator[]( const string &key ) const
    {
        return getEntry( key );
    }
    ParameterEntry &operator[]( const int key )
    {
        return getEntry( key );
    }
    ParameterEntry operator[]( const int key ) const
    {
        return getEntry( key );
    }
    void print( std::ostream &os, const bool singleLine, const int flagMask = 0 ) const;

private:
    ParameterEntry &getEntry( const string &key );
    ParameterEntry getEntry( const string &key ) const;
    ParameterEntry &getEntry( const int key );
    ParameterEntry getEntry( const int key ) const;


protected:
    void addDefaultVerbosityAndHelpEntries();

public:
    void addEntry( int numId, string stringId, string longCmdLineName, string shortCmdLineName, string description, string defaultValue, int flags, const string *inputFormatLabels = NULL );
    void printUsage() const;
    bool parseArgv( const int argc, const char **argv );
    bool chechRequiredParameters();
    void commitDefaultValues();
    void mergeWith( const ToolParameters &other );

    vector<ParameterEntry> entries_;
private:
    void setLoggerVerbosityAndTempDir();
    void printUsageForCategory( int trueMask, int falseMask ) const;
    map<string, int> stringIdMap; // Map stringId -> internalId
    map<int, int> numIdMap; // Map numId -> internalId
};


#endif //ifndef BEETL_TOOL_PARAMETERS_HH
