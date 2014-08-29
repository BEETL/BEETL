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

#include "ErrorInfo.hh"

using namespace std;

bool ErrorInfo::SortByRead( ErrorInfo const & a, ErrorInfo const & b )
{
    // Multi-level comparison to ensure uniqueness of results when using different implementations of std::sort
    if (a.seqNum < b.seqNum)
        return true;
    else if (a.seqNum == b.seqNum)
    {
        if (a.positionInRead < b.positionInRead)
            return true;
        else if (a.positionInRead == b.positionInRead)
        {
            if (a.correctorStart < b.correctorStart)
                return true;
            else if (a.correctorStart == b.correctorStart)
                return ( a.corrector.compare( b.corrector ) < 0 );
        }
    }
    return false;
}


void ErrorInfo::print() const
{
    cout << "Error info: " << endl;
    cout << "   " << "First cycle:" << firstCycle << endl; //first cycle the error is noticed
    cout << "   " << "Last cycle:" << lastCycle << endl; // last cycle where error is noticed
    cout << "   " << "Read End:" << readEnd << endl; //cycle at which $ for read is reached
    cout << "   " << "Corrector: " << corrector << endl; // what read should be corrected to
    cout << "   " << "Sequence number: " << seqNum << endl; // position of read in original list
    cout << endl;
}

void ErrorInfo::SetReadNumbersToOriginal( char *endPosFileName, vector<ErrorInfo> &errorsInSortedReads )
{
    //loop through all the errors and for each one look up which read it comes from
    LetterNumber numchar;
    FILE *InFileEndPos = fopen( endPosFileName, "rb" );
    if ( InFileEndPos == NULL )
    {
        std::cerr << "could not open file!" << endl;
        exit ( EXIT_FAILURE );
    }

    SequenceNumber numText = 0;
    numchar = fread ( &numText, sizeof( SequenceNumber ), 1 , InFileEndPos );
    checkIfEqual( numchar, 1 );
    uint8_t subSequenceCount = 0;
    numchar = fread ( &subSequenceCount, sizeof( uint8_t ), 1 , InFileEndPos );
    checkIfEqual( numchar, 1 );
    uint8_t hasRevComp = 0;
    numchar = fread ( &hasRevComp, sizeof( uint8_t ), 1 , InFileEndPos );
    checkIfEqual( numchar, 1 );

    sortElement triple;
    uint currentSortedReadIndex = 0;
    SequenceNumber i = 0;

    while ( currentSortedReadIndex < errorsInSortedReads.size() )
    {
        numchar = fread ( &triple.seqN, sizeof( SequenceNumber ), 1 , InFileEndPos );
        checkIfEqual( numchar, 1 );
        //        numchar = fread ( &triple.posN, sizeof( LetterNumber ), 1 , InFileEndPos );
        //        checkIfEqual( numchar, 1 );
        //        numchar = fread ( &triple.pileN, sizeof( AlphabetSymbol ), 1 , InFileEndPos );
        //        checkIfEqual( numchar, 1 );
        uint8_t subSequenceNum;
        numchar = fread ( &subSequenceNum, sizeof( uint8_t ), 1 , InFileEndPos );
        checkIfEqual( numchar, 1 );

        while (
            i == ( SequenceNumber )( errorsInSortedReads[currentSortedReadIndex].seqNum )
            &&
            currentSortedReadIndex < errorsInSortedReads.size()
        )
            errorsInSortedReads[currentSortedReadIndex++].seqNum = triple.seqN;
        i++;
    }

    fclose( InFileEndPos );
}

static const char complementaryAlphabet[] = "$TGCNAZ";

string strreverse( const string &inStr )
{
    string result = "";
    for ( int i = inStr.size() - 1; i >= 0; i-- )
        result += inStr[i];
    return result;
}

void ErrorInfo::ConvertRCCorrectionsToOriginal( vector<ErrorInfo> &corrections, int numberOfReads, int readLength )
{
    for ( uint errNo = 0; errNo < corrections.size(); errNo++ )
        if ( corrections[errNo].seqNum >= numberOfReads )
        {
            corrections[errNo].seqNum -= numberOfReads;
            corrections[errNo].positionInRead = readLength - 1 - corrections[errNo].positionInRead;
            for ( uint i = 0; i < corrections[errNo].corrector.size(); i++ )
                corrections[errNo].corrector[i] = complementaryAlphabet[whichPile[( int )corrections[errNo].corrector[i]]];
            corrections[errNo].correctorStart = corrections[errNo].positionInRead;
            corrections[errNo].reverseStrand = true;
        }
        else
        {
            corrections[errNo].correctorStart = corrections[errNo].positionInRead - ( corrections[errNo].corrector.size() - 1 );
            corrections[errNo].corrector = strreverse( corrections[errNo].corrector );
            corrections[errNo].reverseStrand = false;
        }
}

void ErrorInfo::CorrectionsToCsv( const string &fileName, vector<ErrorInfo> &corrections )
{
    ofstream correctionsFile;
    correctionsFile.open( fileName.c_str() );
    correctionsFile << "read position reverse_strand correction corrector_start shortest_witness longest_witness" << endl; //the columns of the output csv file.

    for ( uint errNo = 0; errNo < corrections.size(); errNo++ )
        correctionsFile
                << corrections[errNo].seqNum << " "
                << corrections[errNo].positionInRead << " "
                << corrections[errNo].reverseStrand << " "
                << corrections[errNo].corrector << " "
                << corrections[errNo].correctorStart << " "
                << corrections[errNo].firstCycle << " "
                << corrections[errNo].lastCycle << endl;

    correctionsFile.close();
}

vector<ErrorInfo> ErrorInfo::ReadCorrectionsFromCsv( const string &fileName )
{
    vector<ErrorInfo> result;
    ifstream in( fileName.c_str() );
    string correctionRecord;
    bool firstLine = true;

    while ( getline( in, correctionRecord ) )
    {
        if ( firstLine )
        {
            firstLine = false;
            continue;
        }

        ErrorInfo correction;
        stringstream ss( correctionRecord );

        ss >> correction.seqNum;
        ss >> correction.positionInRead;
        ss >> correction.reverseStrand;
        ss >> correction.corrector;
        ss >> correction.correctorStart;
        ss >> correction.firstCycle;
        ss >> correction.lastCycle;

        correction.readEnd = correction.positionInRead + correction.firstCycle + 1;

        result.push_back( correction );
    }

    in.close();
    return result;
}

