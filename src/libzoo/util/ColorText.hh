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

#ifndef COLOR_TEXT_HH
#define COLOR_TEXT_HH

#include <iostream>
#include <string>

using std::string;


class ColorText
{
public:
    static string startRed;
    static string endRed;

    static void init( int activateColor = -1 ); // -1 = auto: on for tty supporting colors, off for files
};

#endif // COLOR_TEXT_HH
