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

#include "ColorText.hh"

#include <cstring>
#include <unistd.h>

using namespace std;


string ColorText::startRed = "";
string ColorText::endRed = "";

void ColorText::init( int activateColor )
{
    static bool firstTime = true;
    if ( !firstTime )
        return;
    firstTime = false;

    if ( activateColor != 0 )
    {
        char const *t = getenv ( "TERM" );
        if ( activateColor == 1
             || ( isatty ( STDOUT_FILENO ) && ( t == 0 || strcmp ( t, "dumb" ) != 0 ) )
           )
        {
            startRed = "\033[7;31m"; // Note: 7=background to 1=text colour
            endRed = "\033[0m";
        }
    }
}
