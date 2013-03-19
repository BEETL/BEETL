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

#ifndef BEETL_TOOL_PARAMETERS_HH
#define BEETL_TOOL_PARAMETERS_HH

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using std::string;
using std::vector;


static const int MULTIPLE_OPTIONS = 99;

typedef unsigned int ParameterId;
//typedef string (*OptionPossibleValues)[];


class ToolParameters : public vector< int >
{
public:
    virtual ~ToolParameters() {}

    int getValue( const ParameterId key ) const
    {
        if ( key >= this->size() )
            return MULTIPLE_OPTIONS;

        int valNum = this->operator[]( key );
        return valNum;
    }

    string getValueAsString( const ParameterId key ) const
    {
        int valNum = getValue( key );
        return getOptionPossibleValues()[key][valNum];
    }

    void set( const ParameterId key, const int val )
    {
        if ( this->size() <= key )
            this->resize( key + 1, MULTIPLE_OPTIONS );
        this->operator[]( key ) = val;
    }

    void set( const ParameterId key, const string valString )
    {
        const string *possibleValues = getOptionPossibleValues()[key];
        for ( unsigned int i = 0; possibleValues[i] != ""; ++i )
        {
            if ( possibleValues[i] == valString )
            {
                set( key, i );
                return;
            }
        }
        std::cerr << "Error: Invalid value " << valString << " for " << getOptionName( key ) << std::endl;
        exit( 1 );
    }

    const string getOptionPossibleValue( const unsigned i, const unsigned j ) const
    {
        return getOptionPossibleValues()[i][j];
    }

    virtual const string getOptionName( const unsigned i ) const = 0;

private:
    virtual const string **getOptionPossibleValues() const = 0;
};

#endif //ifndef BEETL_TOOL_PARAMETERS_HH
