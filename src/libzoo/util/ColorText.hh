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
