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

#include "CorrectionAligner.hh"

using namespace std;

template <class T>
T **makeMatrix( int rows, int cols )
{
    T **result = new T*[rows];
    for ( int row = 0; row < rows; row++ )
        result[row] = new T[cols];
    return result;
}

template <class T>
void zapMatrix( T **matrix, int rows )
{
    for ( int row = 0; row < rows; row++ )
        delete[] matrix[row];
    delete[] matrix;
}

struct CorrectionInterval
{
    CorrectionInterval( int inStart, int inLengthOnRead, int inCorrectionLength ): correctionLength( inCorrectionLength ), lengthOnRead( inLengthOnRead ), start( inStart ) {}
    int correctionLength, lengthOnRead, start;
};

string CorrectionAligner::MakeFastaRecord( int number, string name, string sequence, string quality )
{
    stringstream ss;
    ss << ">" << number << endl << sequence << endl;
    return ss.str();
}

string CorrectionAligner::MakeFastqRecord( int number, string name, string sequence, string quality )
{
    stringstream ss;
    ss << name << sequence << endl << "+" << endl << quality << endl;
    return ss.str();
}

bool CorrectionAligner::SortByLastCycle( ErrorInfo *a, ErrorInfo *b )
{
    return a->lastCycle > b->lastCycle;
}

void CorrectionAligner::ApplyCorrections(
    SeqReaderFile *readsFile,
    vector<ErrorInfo> &corrections,
    const string &outFile,
    bool correctionsOnly,
    ReadsFormat fileType
)
{
    ofstream correctedReadsFile ( outFile.c_str(), fstream::out );
    ApplyCorrections( readsFile, corrections, correctedReadsFile, correctionsOnly, fileType );
    correctedReadsFile.close();
}

void CorrectionAligner::ApplyCorrections(
    SeqReaderFile *readsFile,
    vector<ErrorInfo> &corrections,
    ostream &correctedReadsOut,
    bool correctionsOnly,
    ReadsFormat fileType
)
{

    uint readLength = readsFile->length();
    readsFile->rewindFile();
    uint currentCorrection = 0;
    int currentRead = 0;

    while ( readsFile->readNext(), !readsFile->allRead() )
    {
        string name = string( readsFile->thisName() );
        string readStr = string( readsFile->thisSeq() ).substr( 0, readLength );

        string qStr = string( readsFile->thisQual() );
        if ( qStr.size() > readLength )
            qStr = qStr.substr( 0, readLength );

        if ( currentCorrection < corrections.size() && corrections[currentCorrection].seqNum == currentRead )
        {
            vector<ErrorInfo *> correctionsForCurrentRead;
            while ( currentCorrection < corrections.size() && corrections[currentCorrection].seqNum == currentRead )
            {
                correctionsForCurrentRead.push_back( &corrections[currentCorrection] );
                ++currentCorrection;
            }

            string correctedRead, correctedQstr;

            CorrectRead( correctionsForCurrentRead, readStr, qStr, correctedRead, correctedQstr );
            if ( fileType == READS_FORMAT_FASTQ )
                correctedReadsOut << MakeFastqRecord( currentRead, name, correctedRead, correctedQstr );
            else if ( fileType == READS_FORMAT_FASTA )
                correctedReadsOut << MakeFastaRecord( currentRead, name, correctedRead, correctedQstr );
        }
        else if ( !correctionsOnly )
        {
            if ( fileType == READS_FORMAT_FASTQ )
                correctedReadsOut << MakeFastqRecord( currentRead, name, readStr, qStr );
            else if ( fileType == READS_FORMAT_FASTA )
                correctedReadsOut << MakeFastaRecord( currentRead, name, readStr, qStr );
        }
        currentRead++;
    }

}

string CorrectionAligner::Correct( const string &errorContainingRead, vector<ErrorInfo *> &corrections )
{
    return errorContainingRead;
}

void CorrectionAligner::CorrectRead(
    vector<ErrorInfo *> &corrections,
    const string &errorContainingRead,
    const string &inQstr,
    string &outRead,
    string &outQstr
)
{
    outRead = Correct( errorContainingRead, corrections );
}

enum AlignType
{
    POSITION_MATCH = 0,
    SEQ1_GAP,
    SEQ2_GAP
};

void SmithWatermanCorrectionAligner::Align( const string &seq1, const string &seq2, int &lengthOnSeq1, int &lengthOnSeq2 )
{
    int **matrix = makeMatrix<int>( seq1.size() + 1, seq2.size() + 1 );
    AlignType **pointers = makeMatrix<AlignType>( seq1.size() + 1, seq2.size() + 1 );

    for ( uint seq1pos = 0; seq1pos <= seq1.size(); seq1pos++ )
        matrix[seq1pos][0] = 0;
    for ( uint seq2pos = 0; seq2pos <= seq2.size(); seq2pos++ )
        matrix[0][seq2pos] = 0;

    for ( uint seq1pos = 1; seq1pos <= seq1.size(); seq1pos++ )
        for ( uint seq2pos = 1; seq2pos <= seq2.size(); seq2pos++ )
        {
            AlignType alignType = POSITION_MATCH;
            int score = 0;
            int matchScore;
            if ( seq1[seq1pos - 1] == seq2[seq2pos - 1] )
                matchScore = matchScore_;
            else
                matchScore = mismatchScore_;

            if ( score < matrix[seq1pos - 1][seq2pos] + deletionScore_ )
            {
                score = matrix[seq1pos - 1][seq2pos] + deletionScore_;
                alignType = SEQ2_GAP;
            }
            if ( score < matrix[seq1pos][seq2pos - 1] + insertionScore_ )
            {
                score = matrix[seq1pos][seq2pos - 1] + insertionScore_;
                alignType = SEQ1_GAP;
            }
            if ( score < matrix[seq1pos - 1][seq2pos - 1] + matchScore )
            {
                score = matrix[seq1pos - 1][seq2pos - 1] + matchScore;
                alignType = POSITION_MATCH;
            }
            matrix[seq1pos][seq2pos] = score;
            pointers[seq1pos][seq2pos] = alignType;
        }

    int pos1 = seq1.size() - 1;
    int pos2 = seq2.size() - 1;

    while ( ( pos1 > 0 ) && ( pos2 > 0 ) )
    {
        switch ( pointers[pos1][pos2] )
        {
            case POSITION_MATCH:
                pos1--;
                pos2--;
                break;
            case SEQ1_GAP:
                pos2--;
                break;
            case SEQ2_GAP:
                pos1--;
                break;
            default:
                exit( 1 );
        }
    }

    lengthOnSeq1 = seq1.size() - pos1;
    lengthOnSeq2 = seq2.size() - pos2;

    zapMatrix<AlignType>( pointers, seq1.size() + 1 );
    zapMatrix<int>( matrix, seq1.size() + 1 );
}

void SmithWatermanCorrectionAligner::Align( const string &seq1, const string &seq2, int &lengthOnSeq1, int &lengthOnSeq2, bool correctForwards )
{
    if ( correctForwards )
    {
        string seq1_ = strreverse( seq1 );
        string seq2_ = strreverse( seq2 );

        Align( seq1_, seq2_, lengthOnSeq1, lengthOnSeq2 );
    }
    else
        Align( seq1, seq2, lengthOnSeq1, lengthOnSeq2 );
}

string SmithWatermanCorrectionAligner::Replace( const string &original, const string &correction, int lineUpPosition, bool correctForwards )
{
    int lengthOnOriginal;
    return Replace( original, correction, lineUpPosition, correctForwards, lengthOnOriginal );
}

string SmithWatermanCorrectionAligner::Replace( const string &original, const string &correction, int lineUpPosition, bool correctForwards, int &lengthOnOriginal )
{
    int lengthOnCorrection;
    string originalPartToAlign = correctForwards ? original.substr( lineUpPosition ) : original.substr( 0, lineUpPosition + 1 );
    Align( correction, originalPartToAlign, lengthOnCorrection, lengthOnOriginal, correctForwards );
    if ( correctForwards )
    {
        return
            original.substr( 0, lineUpPosition )
            +
            correction
            +
            original.substr( min<int>( lengthOnOriginal + lineUpPosition, original.size() ) )
            ;
    }
    else
    {
        return
            original.substr( 0, max<int>( lineUpPosition - lengthOnOriginal + 1, 0 ) )
            +
            correction
            +
            original.substr( min<int>( original.size(), lineUpPosition + 1 ) )
            ;
    }
}

string SmithWatermanCorrectionAligner::Correct( const string &errorContainingRead, vector<ErrorInfo *> &corrections )
{
    bool firstCorrection = true;
    string result( errorContainingRead );
    for ( uint currentCorrection = 0; currentCorrection < corrections.size(); currentCorrection++ )
    {
        ErrorInfo *current = corrections[currentCorrection];
        string witness = ( current->reverseStrand ) ?
                         errorContainingRead.substr( current->positionInRead - current->lastCycle, current->lastCycle ) :
                         errorContainingRead.substr( current->positionInRead + 1, current->lastCycle );

        int lineUpAt;
        int canCorrect = true;
        if ( firstCorrection )
        {
            //if its the first correction we're applying to the read then we know the exact position the putative error
            //occurred, so we can line up exactly and avoid calling str.find
            lineUpAt = current->positionInRead;
        }
        else
        {
            int witnessLocation = result.find( witness );
            if ( witnessLocation == -1 )
                canCorrect = false;
            lineUpAt = ( current->reverseStrand ) ? witnessLocation + current->lastCycle : witnessLocation - 1;
        }
        if ( canCorrect )
        {
            result = Replace(
                         result,
                         current->corrector,
                         lineUpAt,
                         current->reverseStrand
                     );
        }
        firstCorrection = false;
    }

    return result;

}

string StitchAligner::Correct( const string &errorContainingRead, vector<ErrorInfo *> &corrections )
{
    return errorContainingRead.substr( 0, corrections[0]->correctorStart ) + corrections[0]->corrector;
}

void NoIndelAligner::CorrectRead(
    vector<ErrorInfo *> &corrections,
    const string &errorContainingRead,
    const string &inQstr,
    string &outRead,
    string &outQstr
)
{
    bool firstCorrection = true;
    outRead = errorContainingRead;
    outQstr = inQstr;

    sort( corrections.begin(), corrections.end(), SortByLastCycle );

    for ( uint currentCorrection = 0; currentCorrection < corrections.size(); currentCorrection++ )
    {
        ErrorInfo *current = corrections[currentCorrection];

        bool canCorrect = true;

        if ( current->lastCycle < minLastCycle_ )
            canCorrect = false;

        string witness = ( current->reverseStrand ) ?
                         errorContainingRead.substr( current->positionInRead - current->lastCycle, current->lastCycle ) :
                         errorContainingRead.substr( current->positionInRead + 1, current->lastCycle );

        int lineUpAt;
        if ( firstCorrection )
            lineUpAt = current->positionInRead;
        else
        {
            int witnessLocation = outRead.find( witness );
            if ( witnessLocation == -1 )
                canCorrect = false;
            lineUpAt = ( current->reverseStrand ) ? witnessLocation + current->lastCycle : witnessLocation - 1;
        }

        if ( canCorrect )
        {
            string corrector = current->corrector;
            //int corrStart = current->correctorStart;
            int corrStart = ( current->reverseStrand ) ? lineUpAt : lineUpAt - corrector.size() + 1;

            outRead =
                outRead.substr( 0, max<int>( corrStart, 0 ) ) +
                corrector +
                outRead.substr( min<int>( corrStart + corrector.size(), outRead.size() ) )
                ;

            outQstr =
                outQstr.substr( 0, max<int>( corrStart, 0 ) ) +
                string( corrector.size(), correctionQuality_ ) +
                outQstr.substr( min<int>( corrStart + corrector.size(), outQstr.size() ) )
                ;

            if ( trim_ )
            {
                outRead = outRead.substr( max<int>( 0, -corrStart ), errorContainingRead.size() );
                outQstr = outQstr.substr( max<int>( 0, -corrStart ), errorContainingRead.size() );
            }

        }
        firstCorrection = false;
    }
}
