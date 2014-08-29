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

#include "AlignCorrectorStrings.hh"

#include "BCRext.hh"
#include "BCRexternalBWT.hh"
#include "Common.hh"
#include "DatasetMetadata.hh"
#include "parameters/BwtParameters.hh"
#include "config.h"
#include "libzoo/cli/Common.hh"
#include "libzoo/util/Logger.hh"
#include "errors/ErrorInfo.hh"
#include "errors/WitnessReader.hh"
#include "errors/AlignmentParameters.hh"
#include "errors/BwtCorrector.hh"
#include "errors/CorrectionAligner.hh"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

int main( const int argc, const char **argv )
{
    AlignmentParameters params;
    if ( !params.parseArgv( argc, argv ) || params["help"] == 1 || !params.chechRequiredParameters() )
    {
        params.printUsage();
        exit( params["help"] == 0 );
    }

    //read the corrections to be applied out of the corrections file...
    vector<ErrorInfo> corrections = ErrorInfo::ReadCorrectionsFromCsv( params.getStringValue( "input corrections file" ) );
    std::sort( corrections.begin(), corrections.end(), ErrorInfo::SortByRead );

    cout << "Attempting to make " << corrections.size() << " corrections..." << endl;
    unique_ptr<CorrectionAligner> aligner;

    int alignmentType = params.getValue( "alignment type" );

    if ( alignmentType == ALIGNMENT_TYPE_SW )
    {
        if ( !params["mismatch penalty"].isSet() || !params["deletion penalty"].isSet() || !params["insertion penalty"].isSet() )
        {
            cout << "Smith waterman alignment requires you to set the mismatch/deletion/insertion penalties and the witness length used to generate the corrections..." << endl;
            exit( 1 );
        }
        cout << "Using Smith-Water man local alignment to position the corrections..." << endl;
        aligner.reset( new SmithWatermanCorrectionAligner(
                           2,
                           params.getValue( "mismatch penalty" ),
                           params.getValue( "deletion penalty" ),
                           params.getValue( "insertion penalty" )
                       )
                     );
    }
    else if ( alignmentType == ALIGNMENT_TYPE_NO_INDELS )
    {
        if ( !params["correction quality"].isSet() || !params["min witness length"] )
        {
            cout << "Alignment with no indels requires you to set the witness length and min witness length parameter, and correction string quality parameter..." << endl;
            exit( 1 );
        }
        cout << "Superimposing correction strings onto reads (without indels) and trimming to original length" << endl;
        aligner.reset( new NoIndelAligner(
                           params.getStringValue( "correction quality" )[0],
                           params.getValue( "min witness length" ),
                           ( params.getValue( "trim corrected reads" ) == 1 )
                       )
                     );

    }
    else if ( alignmentType == ALIGNMENT_TYPE_STITCH )
    {
        aligner.reset( new StitchAligner() );
    }
    else
    {
        cerr << "Error: unexpected alignment type" << endl;
        assert( false );
    }

    string readsFileName = params.getStringValue( "input reads file" );
    FILE *reads = fopen( readsFileName.c_str(), "r" );

    string outputReadsFile = params.getStringValue( "corrected reads output file" );
    cout << "Writing corrector-aligned reads to " << outputReadsFile << "..." << endl;

    SeqReaderFile *readsFile = NULL;
    switch ( params.getValue( "input reads format" ) )
    {
        case READS_FORMAT_FASTA:
            readsFile = new SeqReaderFasta( reads );
            break;
        case READS_FORMAT_FASTQ:
            readsFile = new SeqReaderFastq( reads );
            break;
        default:
            cout << "Unsupported input reads file type" << endl;
            exit( 1 );
    }

    ReadsFormat outFormat;

    switch ( params.getValue( "output reads format" ) )
    {
        case READS_FORMAT_FASTQ:
            outFormat = READS_FORMAT_FASTQ;
            break;
        case READS_FORMAT_FASTA:
            outFormat = READS_FORMAT_FASTA;
            break;
        default:
            cout << "Unsupported reads format!" << endl;
            exit( 1 );
    }

    aligner->ApplyCorrections( readsFile, corrections, outputReadsFile, false, outFormat );

    fclose( reads );
    delete readsFile;
    cout << "Done" << endl;
    return 0;
}
