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

#include "HiTECErrorLocation.hh"

using namespace std;

bool HiTECErrorLocation::ErrorLocationSorter( HiTECErrorLocation a, HiTECErrorLocation b )
{
    if ( a.readNum == b.readNum )
        return a.positionInRead < b.positionInRead;
    else
        return a.readNum < b.readNum;
}

void HiTECErrorLocation::SetReadNumbersToOriginal( char *endPosFileName, vector<HiTECErrorLocation> &errorsInSortedReads )
{
    //loop through all the errors and for each one look up which read it comes from

    LetterNumber numchar = 0;
    FILE *InFileEndPos;                  // input file of the end positions;
    InFileEndPos = fopen( endPosFileName, "rb" );
    if ( InFileEndPos == NULL )
    {
        std::cerr << "could not open file!" << endl;
        exit ( EXIT_FAILURE );
    }

    SequenceNumber numText = 0;
    numchar = fread ( &numText, sizeof( SequenceNumber ), 1 , InFileEndPos );

    numchar = 0;
    sortElement triple;

    int currentSortedReadIndex = 0;
    SequenceNumber i = 0;

    while ( currentSortedReadIndex < errorsInSortedReads.size() )
    {
        numchar = fread ( &triple.seqN, sizeof( SequenceNumber ), 1 , InFileEndPos );
        checkIfEqual( numchar, 1 ); // we should always read the same number of characters
        numchar = fread ( &triple.posN, sizeof( LetterNumber ), 1 , InFileEndPos ); //it is the relative position of the $ in the partial BWT
        checkIfEqual( numchar, 1 ); // we should always read the same number of characters
        numchar = fread ( &triple.pileN, sizeof( AlphabetSymbol ), 1 , InFileEndPos );
        checkIfEqual( numchar, 1 ); // we should always read the same number of characters

        while (
            i == errorsInSortedReads[currentSortedReadIndex].readNum
            &&
            currentSortedReadIndex < errorsInSortedReads.size()
        )
            errorsInSortedReads[currentSortedReadIndex++].readNum = triple.seqN;
        i++;
    }

    fclose( InFileEndPos );
}

static const char complementaryAlphabet[] = "$TGCNAZ";

void HiTECErrorLocation::ConvertRCCorrectionsToOriginal( vector<HiTECErrorLocation> &errors, int numberOfReads, int readLength )
{
    for ( int errNo = 0; errNo < errors.size(); errNo++ )
        if ( errors[errNo].readNum >= numberOfReads )
        {
            errors[errNo].readNum -= numberOfReads;
            errors[errNo].positionInRead = readLength - 1 - errors[errNo].positionInRead;
        }
        else
        {
            string corr = "";
            corr += errors[errNo].corrector;
            errors[errNo].corrector = corr[0];
        }
}

void HiTECErrorLocation::CorrectionsToCsv( const string &fileName, vector<HiTECErrorLocation> &errorLocations )
{
    ofstream correctionsFile;
    correctionsFile.open( fileName.c_str() );
    correctionsFile << "read position correction" << endl; //the columns of the output csv file.

    for ( int errNo = 0; errNo < errorLocations.size(); errNo++ )
        correctionsFile
                << errorLocations[errNo].readNum << " "
                << errorLocations[errNo].positionInRead << " "
                << errorLocations[errNo].corrector << endl;

    correctionsFile.close();
}

vector<HiTECErrorLocation> HiTECErrorLocation::ReadCorrectionsFromCsv( const string &fileName )
{
    vector<HiTECErrorLocation> result;
    ifstream in( fileName.c_str() );
    string correctionRecord;
    bool firstLine = true;

    int readNum, positionInRead;
    char corrector;

    while ( getline( in, correctionRecord ) )
    {
        if ( firstLine )
        {
            firstLine = false;
            continue;
        }

        stringstream ss( correctionRecord );

        ss >> readNum;
        ss >> positionInRead;
        ss >> corrector;

        HiTECErrorLocation rel = HiTECErrorLocation( readNum, positionInRead, corrector );
        result.push_back( rel );
    }

    in.close();
    return result;
}

