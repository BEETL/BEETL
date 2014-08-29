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

#include "HiTECStats.hh"

#include "BCRext.hh"
#include "BCRexternalBWT.hh"

#include "Timer.hh"
#include "config.h"

#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#ifdef HAVE_POSIX_FADVISE
#endif

using namespace std;
using namespace BeetlBwtParameters;

HiTECStats::HiTECStats(
    const double errorRate,
    const double genomeLength,
    const int numberOfReads,
    const int readLength
):
    numberOfReads_( numberOfReads ),
    readLength_( readLength ),
    genomeLength_( genomeLength ),
    errorRate_( errorRate )
{

}

double nCr( int n, int r )
{
    //use logs to avoid overflow
    double result = 0;
    for ( int i = r + 1; i <= n; i++ )
        result += log( i );
    for ( int i = 1; i <= n - r; i++ )
        result -= log( i );
    return exp( result );
}

double HiTECStats::probDestructible( int witnessLength )
{
    double result = 1;
    //probability witness has errors
    result *= (
                  1 - pow( 1 - errorRate_, witnessLength )
              );
    //probability that support is correct
    result *= ( 1 - errorRate_ );

    //probability occurs in genome
    result *= 0.75 * (
                  1 - pow(
                      1 - pow( 0.25, witnessLength ),
                      genomeLength_
                  )
              );
    return result;
}

double HiTECStats::expectedDestructible( int witnessLength )
{
    double result = 1;
    result *= (
                  1 - pow(
                      1 - probDestructible( witnessLength ),
                      readLength_ - witnessLength
                  )
              );
    result *= pow( 1 - errorRate_, readLength_ );
    result *= numberOfReads_;
    return result;
}

double HiTECStats::expectedUncorrectable( int witnessLength )
{
    //create matrix of 'uncorrectable read counts'...
    //this is essentially an implementation of the formula:
    //f_w(k,l) as found in the HiTEC paper
    double **storedErrors = new double*[readLength_ + 1];
    for ( int i = 0; i < readLength_ + 1; i++ )
        storedErrors[i] = new double[readLength_ + 1];

    for ( int length = 0; length < witnessLength; length++ )
        for ( int numErrors = 0; numErrors <= length; numErrors++ )
            storedErrors[numErrors][length] = nCr( length, numErrors );
    for ( int length = witnessLength; length <= readLength_; length++ )
        for ( int numErrors = 0; numErrors <= length; numErrors++ )
            if ( numErrors >= ( int )length / witnessLength )
            {
                int sum = 0;
                for ( int knockOff = 1; knockOff <= witnessLength; knockOff++ )
                    sum += storedErrors[numErrors - 1][length - knockOff];
                storedErrors[numErrors][length] = sum;
            }

    double result = 0;
    for ( int numErrors = 0; numErrors <= readLength_; numErrors++ )
        result += storedErrors[numErrors][readLength_] * pow( errorRate_, numErrors ) * pow( 1 - errorRate_, readLength_ - numErrors );
    for ( int i = 0; i < readLength_ + 1; i++ )
        delete[] storedErrors[i];
    delete[] storedErrors;

    return result * numberOfReads_;
}

double HiTECStats::expectedCorrectWitnessNeighbourPairs( int support, int witnessLength )
{
    double qc = ( readLength_ - witnessLength ) * pow( 1 - errorRate_, witnessLength + 1 ) / genomeLength_;
    return nCr( numberOfReads_, support ) * pow( qc, support ) * pow( 1 - qc, numberOfReads_ - support ) * genomeLength_;
}

double HiTECStats::expectedIncorrectWitnessNeighbourPairs( int support, int witnessLength )
{
    double qe = ( readLength_ - witnessLength ) * pow( 1 - errorRate_, witnessLength ) * errorRate_ / 3 / genomeLength_;
    return nCr( numberOfReads_, support ) * pow( qe, support ) * pow( 1 - qe, numberOfReads_ - support ) * genomeLength_;
}

double HiTECStats::expectedErroneousReads()
{
    return numberOfReads_ * ( 1 - pow( 1 - errorRate_, readLength_ ) );
}

int HiTECStats::Calculate_wm()
{
    return Calculate_wm( readLength_ );
}

int HiTECStats::Calculate_wm( int maxValue )
{
    double lowestScore = ( double )numberOfReads_;
    int best_wm = 1;
    for ( int candidate_wm = 1; candidate_wm <= maxValue; candidate_wm++ )
    {
        double candidateScore = expectedDestructible( candidate_wm ) + expectedUncorrectable( candidate_wm );
        if ( candidateScore < lowestScore )
        {
            lowestScore = candidateScore;
            best_wm = candidate_wm;
        }
    }
    return best_wm;
}

int HiTECStats::Calculate_wM( double maxDestructibleRate )
{
    //return the smallest witness length (at least wm) which gives desired specificity
    for ( int candidate_Wm = 1; candidate_Wm <= readLength_; candidate_Wm++ )
        if ( expectedDestructible( candidate_Wm ) < maxDestructibleRate * expectedErroneousReads() )
            return candidate_Wm;
    return readLength_;
}

int HiTECStats::CalculateSupport( int witnessLength )
{
    int candidateThreshold = 1;
    while (
        expectedCorrectWitnessNeighbourPairs( candidateThreshold, witnessLength )
        <
        expectedIncorrectWitnessNeighbourPairs( candidateThreshold, witnessLength )
        &&
        candidateThreshold < numberOfReads_
    )
        candidateThreshold++;

    return candidateThreshold + 2;
}
