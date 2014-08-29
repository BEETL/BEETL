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

#include "ToolParameters.hh"

#include "libzoo/cli/Common.hh"
#include "libzoo/util/Logger.hh"
#include "libzoo/util/TemporaryFilesManager.hh"

#include <cassert>
#include <cstring>
#include <sstream>

using namespace std;


bool ParameterEntry::operator==( const int rhs )
{
    return ( parsedValue == rhs );
}

bool ParameterEntry::operator==( const string &rhs )
{
    return ( userValue == rhs );
}

ParameterEntry::operator int () const
{
    return parsedValue;
}

ParameterEntry::operator string () const
{
    if ( flags & TYPE_INT )
    {
        ostringstream ss;
        ss << parsedValue;
        return ss.str();
    }
    else if ( flags & TYPE_SWITCH )
    {
        switch ( parsedValue )
        {
            case 0:
                return "off";
            case 1:
                return "on";
            default:
                assert( false );
        }
    }
    else if ( flags & TYPE_CHOICE )
    {
        assert( parsedValue >= 0 );
        return possibleValues[parsedValue];
    }
    return userValue;
}

ParameterEntry &ParameterEntry::operator=( const int rhs )
{
    set( rhs );
    return *this;
}

ParameterEntry &ParameterEntry::operator=( const string &rhs )
{
    set( rhs );
    return *this;
}

bool ParameterEntry::isSet() const
{
    if ( flags & TYPE_SWITCH )
        return ( parsedValue == 1 ); // SWITCHes always get a default value
    else
        return !userValue.empty();
}

void ParameterEntry::set( const string &valString )
{
    if ( !userValue.empty() )
    {
        cerr << "Warning: Overwriting value for parameter \"" << stringId << "\" from " << userValue << " to \"" << valString << "\"" << endl;
    }
    userValue = valString;

    // If this parameter had multiple-choice values, make sure one of them matches the user's
    if ( flags & TYPE_CHOICE )
    {
        assert ( possibleValues );

        bool valueFound = false;
        for ( int j = 0; !possibleValues[j].empty(); ++j )
        {
            if ( strcasecmp ( valString.c_str(), possibleValues[j].c_str() ) == 0 )
            {
                parsedValue = j;
                valueFound = true;
                break;
            }
        }
        if ( !valueFound )
        {
            cerr << "Error: Invalid parameter value \"" << valString << "\" for parameter \"" << stringId << "\"\n" << endl;
            exit( -1 );
        }
    }

    // Integer types store their main value in parsedValue
    if ( flags & TYPE_INT )
    {
        istringstream ss( userValue );
        ss >> parsedValue;
    }
}

void ParameterEntry::set( const int val, const bool isSilent )
{
    if ( parsedValue != -1 && !isSilent )
    {
        cerr << "Warning: Overwriting value for parameter \"" << stringId << "\" from " << parsedValue << " to " << val << endl;
    }
    parsedValue = val;

    // userValue also acts as a flag to know if the value has been set
    if ( userValue.empty() )
    {
        userValue = "set";
    }
}


ostream &operator<<( std::ostream &os, const ParameterEntry &obj )
{
    int spaces = 25;
    os << "    " << obj.longCmdLineName;
    spaces -= obj.longCmdLineName.length();
    if ( obj.shortCmdLineName != "" )
    {
        os << " (" << obj.shortCmdLineName << ")";
        spaces -= 3 + obj.shortCmdLineName.length();
    }
    if ( spaces > 0 )
    {
        os << string( spaces, ' ' );
    }

    spaces = 85;
    if ( !obj.description.empty() )
    {
        os << obj.description << ' ';
        spaces -= obj.description.length() + 1;
    }

    if ( obj.flags & TYPE_CHOICE )
    {
        os << '[';
        if ( !obj.possibleValues[0].empty() )
        {
            for ( int i = 0; ; ++i )
            {
                const string &val = obj.possibleValues[i];
                if ( val.empty() ) break;

                if ( i )
                    os << '|';
                os << val;
                spaces -= val.length() + 1;
            }
            ++spaces;
        }
        os << ']';
        spaces -= 2;
    }

    if ( spaces > 0 )
    {
        os << string( spaces, ' ' );
    }
    if ( obj.defaultValue != "" )
        os << " (Default: " << obj.defaultValue << ")";

    return os;
}


//////////////////////////////////////////


ParameterEntry &ToolParameters::getEntry( const string &key )
{
    map<string, int>::iterator it = stringIdMap.find( key );
    if ( it != stringIdMap.end() )
    {
        int internalId = it->second;
        return entries_[ internalId ];
    }
    std::cerr << "Error: Invalid parameter key \"" << key << "\"" << std::endl;
    assert( false );
}

ParameterEntry &ToolParameters::getEntry( const int key )
{
    map<int, int>::iterator it = numIdMap.find( key );
    if ( it != numIdMap.end() )
    {
        int internalId = it->second;
        return entries_[ internalId ];
    }
    std::cerr << "Error: Invalid parameter key \"" << key << "\"" << std::endl;
    assert( false );
}

ParameterEntry ToolParameters::getEntry( const string &key ) const
{
    map<string, int>::const_iterator it = stringIdMap.find( key );
    if ( it != stringIdMap.end() )
    {
        int internalId = it->second;
        return entries_[ internalId ];
    }
    std::cerr << "Error: Invalid parameter key \"" << key << "\"" << std::endl;
    assert( false );
}

ParameterEntry ToolParameters::getEntry( const int key ) const
{
    map<int, int>::const_iterator it = numIdMap.find( key );
    if ( it != numIdMap.end() )
    {
        int internalId = it->second;
        return entries_[ internalId ];
    }
    std::cerr << "Error: Invalid parameter key \"" << key << "\"" << std::endl;
    assert( false );
}

int ToolParameters::getValue( const int key ) const
{
    const ParameterEntry &entry = getEntry( key );
    return entry.parsedValue;
}

int ToolParameters::getValue( const string &key ) const
{
    const ParameterEntry &entry = getEntry( key );
    return entry.parsedValue;
}

string ToolParameters::getStringValue( const int key ) const
{
    const ParameterEntry &entry = getEntry( key );

    if ( entry.parsedValue == -1 )
        return entry.userValue;

    return entry.possibleValues[entry.parsedValue];
}

string ToolParameters::getStringValue( const string &key ) const
{
    const ParameterEntry &entry = getEntry( key );

    if ( entry.parsedValue == -1 )
        return entry.userValue;

    return entry.possibleValues[entry.parsedValue];
}

void ToolParameters::addEntry( int numId, string stringId, string longCmdLineName, string shortCmdLineName, string description, string defaultValue, int flags, const string *possibleValues )
{
    // If a value is "automated", there shouldn't be a default value
    if ( flags & AUTOMATED )
    {
        assert( defaultValue.empty() );
        defaultValue = "auto";
    }

    int internalId = entries_.size();
    entries_.push_back( ParameterEntry( stringId,  numId,  longCmdLineName,  shortCmdLineName,  description,  defaultValue, flags, possibleValues, this ) );

    // Map stringId -> internalId
    stringIdMap[stringId] = internalId;

    // Map numId -> internalId
    numIdMap[numId] = internalId;
}

void ToolParameters::printUsage() const
{
    cout << "Parameters:" << endl;
    printUsageForCategory( REQUIRED, ENVIRONMENT );
    cout << endl;
    cout << "Options:" << endl;
    printUsageForCategory( 0, REQUIRED | ENVIRONMENT );
    cout << endl;
    cout << "Environment:" << endl;
    printUsageForCategory( ENVIRONMENT, REQUIRED );
    cout << endl;
}


void ToolParameters::addDefaultVerbosityAndHelpEntries()
{
    addEntry( -1, "temp directory", "--temp-directory", "-T", "Path for temporary files (hint: choose a fast drive)", ".", TYPE_STRING | ENVIRONMENT );
    addEntry( -1, "no temp subdir", "--no-temp-subdir", "", "Prevent creation of a uniquely named temporary sub-directory", "", TYPE_SWITCH | ENVIRONMENT );
    addEntry( -1, "use shm", "--use-shm", "", "Use shared memory across processes (faster initialisation of beetl-search and beetl-extend) e.g. --use-shm=/dev/shm", "", TYPE_STRING | ENVIRONMENT );
    addEntry( -1, "use color", "--color", "", "For beetl-extend to highlight the matching k-mers", "auto", TYPE_CHOICE | ENVIRONMENT, colorLabels );
    addEntry( -1, "verbosity", "--verbosity", "", "[quiet|normal|verbose|very-verbose|debug] or [0|1|2|3|4]", "normal", TYPE_STRING | ENVIRONMENT );
    addEntry( -1, "verbose", "", "-v", "Shortcut to --verbosity = verbose", "", TYPE_SWITCH | ENVIRONMENT );
    addEntry( -1, "very-verbose", "", "-vv", "Shortcut to --verbosity = very-verbose", "", TYPE_SWITCH | ENVIRONMENT );
    addEntry( -1, "help", "--help", "-h", "Help", "", TYPE_SWITCH | ENVIRONMENT );
}

void ToolParameters::printUsageForCategory( int trueMask, int falseMask ) const
{
    for ( vector<ParameterEntry>::const_iterator it = entries_.begin(); it != entries_.end(); ++it )
    {
        if ( ( it->flags & trueMask ) == trueMask && ( ~it->flags & falseMask ) == falseMask )
        {
            cout << *it << endl;
        }
    }
}


bool ToolParameters::parseArgv( const int argc, const char **argv )
{
    for ( int i = 1; i < argc; ++i )
    {
        bool found = false;
        for ( vector<ParameterEntry>::iterator it = entries_.begin(); it != entries_.end(); ++it )
        {
            bool isSwitch = ( it->flags & TYPE_SWITCH );
            string val;
            if ( isNextArgument( it->shortCmdLineName, it->longCmdLineName, argc, argv, i, isSwitch ? NULL : &val ) )
            {
                if ( isSwitch )
                    *it = 1;
                else
                    *it = val;
                found = true;
                break;
            }
        }

        if ( !found )
        {
            cerr << "Error: Invalid parameter: " << argv[i] << "\n" << endl;
            return false;
        }
    }

    setLoggerVerbosityAndTempDir();
    return true;
}

void ToolParameters::setLoggerVerbosityAndTempDir()
{
    if ( getEntry( "verbose" ) == 1 )
    {
        Logger::setVerbosity( "verbose" );
    }
    if ( getEntry( "very-verbose" ) == 1 )
    {
        Logger::setVerbosity( "very-verbose" );
    }
    if ( getEntry( "verbosity" ).isSet() )
    {
        Logger::setVerbosity( getEntry( "verbosity" ) );
    }

    // Initialise temporary directory
    TemporaryFilesManager::get().setTempPath( getEntry( "temp directory" ), getEntry( "no temp subdir" ) != true );
}

void ToolParameters::print( std::ostream &os, const bool singleLine, const int flagMask ) const
{
    if ( singleLine )
    {
        os << "{ ";
        for ( vector<ParameterEntry>::const_iterator it = entries_.begin(); it != entries_.end(); ++it )
            if ( it->flags & AUTOMATED )
                os << it->stringId << " = " << static_cast<int>( *it ) << "/\"" << static_cast<string>( *it ) << "\", ";;
        os << "}";
    }
    else
        for ( unsigned int i = 0; i < entries_.size(); ++i )
        {
            Logger::out() << "  " << entries_[i].stringId << " = " << static_cast<string>( entries_[i] ) << endl;;
            /*
                    if ( config.first[i] == MULTIPLE_OPTIONS )
                    {
                        Logger::out() << "* => ";
                        config.first[i] = 0;
                    }
                    Logger::out() << getOptionPossibleValue( i, config.first[i] ) << endl;
            */
        }
}

bool ToolParameters::chechRequiredParameters()
{
    bool checkPassed = true;
    for ( vector<ParameterEntry>::iterator it = entries_.begin(); it != entries_.end(); ++it )
    {
        if ( ( it->flags & REQUIRED ) && !it->isSet() && it->defaultValue.empty() )
        {
            cerr << "Error: Missing argument: " << it->longCmdLineName << " is required" << endl;
            checkPassed = false;
        }
    }
    if ( !checkPassed )
        cerr << endl;
    return checkPassed;
}

void ToolParameters::commitDefaultValues()
{
    for ( vector<ParameterEntry>::iterator it = entries_.begin(); it != entries_.end(); ++it )
    {
        if ( !it->isSet() && !it->defaultValue.empty() && it->defaultValue != "detect" && !( it->flags & AUTOMATED ) )
        {
            *it = it->defaultValue;
        }
        else if ( ( it->flags & TYPE_SWITCH ) && ( it->parsedValue == -1 ) )
        {
            // switches have default value 0
            *it = 0;
        }
    }
}

void ToolParameters::mergeWith( const ToolParameters &other )
{
    assert( false && "todo" );
    /*
        for ( unsigned int j = 0; j < allEstimates[i - 1].first.size(); ++j )
            if ( allEstimates[i - 1].first[j] != allEstimates[i].first[j] )
                allEstimates[i - 1].first[j] = MULTIPLE_OPTIONS;
    */
}
