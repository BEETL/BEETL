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

#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <math.h>
#include <limits>

#include "AlignmentParameters.hh"
#include "ErrorInfo.hh"

#ifndef INCLUDED_CORRECTIONALIGNER_HH
#define INCLUDED_CORRECTIONALIGNER_HH

using namespace std;
using namespace BeetlAlignParameters;

template <class T>
void zapMatrix( T **matrix, int rows );
template <class T>
T **makeMatrix( int rows, int cols );

class CorrectionAligner
{
public:
    virtual ~CorrectionAligner() {}
    static string MakeFastaRecord( int currentRead, string name, string sequence, string quality );
    static string MakeFastqRecord( int currentRead, string name, string sequence, string quality );
    static bool SortByLastCycle( ErrorInfo *a, ErrorInfo *b );

    void ApplyCorrections(
        SeqReaderFile *readsFile,
        vector<ErrorInfo> &corrections,
        ostream &correctedReadsOut,
        bool correctionsOnly,
        ReadsFormat fileType
    );

    void ApplyCorrections(
        SeqReaderFile *readsFile,
        vector<ErrorInfo> &corrections,
        const string &outFile,
        bool correctionsOnly,
        ReadsFormat fileType
    );

    //override this method to make a new aligner capable of correcting reads without quality scores
    virtual string Correct( const string &errorContainingRead, vector<ErrorInfo *> &corrections );

    //override this method to make a new aligner capable of correcting reads with quality scores
    virtual void CorrectRead(
        vector<ErrorInfo *> &corrections,
        const string &errorContainingRead,
        const string &inQstr,
        string &outRead,
        string &outQstr
    );

};

class SmithWatermanCorrectionAligner : public CorrectionAligner
{
public:
    SmithWatermanCorrectionAligner( int m, int mm, int d, int i ): matchScore_( m ), mismatchScore_( mm ), deletionScore_( d ), insertionScore_( i ) {}

    void Align( const string &seq1, const string &seq2, int &lengthOnSeq1, int &lengthOnSeq2 );
    void Align( const string &seq1, const string &seq2, int &lengthOnSeq1, int &lengthOnSeq2, bool correctForwards );
    string Replace( const string &original, const string &correction, int lineUpPosition, bool correctForwards );
    string Replace( const string &original, const string &correction, int lineUpPosition, bool correctForwards, int &lengthOnOriginal );
    string Correct( const string &errorContainingRead, vector<ErrorInfo *> &corrections );

private:
    int matchScore_;
    int mismatchScore_;
    int deletionScore_;
    int insertionScore_;
};

class StitchAligner : public CorrectionAligner
{
public:
    string Correct( const string &errorContainingRead, vector<ErrorInfo *> &corrections );
};

class NoIndelAligner : public CorrectionAligner
{
public:
    NoIndelAligner( char correctionQuality, int minLastCycle, bool trim ):
        correctionQuality_( correctionQuality ),
        minLastCycle_( minLastCycle ),
        trim_( trim )
    {}

    void CorrectRead(
        vector<ErrorInfo *> &corrections,
        const string &errorContainingRead,
        const string &inQstr,
        string &outRead,
        string &outQstr
    );

private:
    char correctionQuality_;
    int minLastCycle_;
    bool trim_;
};

#endif
