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

#include "BeetlBwt.hh"

#include "BCRexternalBWT.hh"
#include "Common.hh"
#include "DatasetMetadata.hh"
#include "Logger.hh"
#include "parameters/BwtParameters.hh"
#include "config.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;


class HardwareConstraints
{
public:
    /*
        float diskSpeed;
        float ramSpeed;
        float ramCongestionLimit;
    */
    float diskMaxSpeedPerProcess; // in MB/s
    float diskMaxTotalSpeed; // in MB/s
    float ramMaxSpeedPerProcess;  // in MB/s
    float ramMaxTotalSpeed;  // in MB/s

    HardwareConstraints()
    // Default values
        : diskMaxSpeedPerProcess( 50 ) // in MB/s
        , diskMaxTotalSpeed( 100 ) // in MB/s
        , ramMaxSpeedPerProcess( 50 ) // in MB/s
        , ramMaxTotalSpeed( 100 ) // in MB/s
    {}

    void init( const string &hardwareConstraintsFilename )
    {
        if ( !hardwareConstraintsFilename.empty() )
        {
            ifstream is( hardwareConstraintsFilename.c_str() );
            string line;
            while ( getline( is, line ) )
            {
                searchForKeywordInString( line, "diskMaxSpeedPerProcess=", diskMaxSpeedPerProcess );
                searchForKeywordInString( line, "diskMaxTotalSpeed=", diskMaxTotalSpeed );
                searchForKeywordInString( line, "ramMaxSpeedPerProcess=", ramMaxSpeedPerProcess );
                searchForKeywordInString( line, "ramMaxTotalSpeed=", ramMaxTotalSpeed );
            }
        }
    }

private:
    void searchForKeywordInString( const string &line, const string &keyword, float &result )
    {
        if ( line.substr( 0, keyword.length() ) == keyword )
        {
            string rhs = line.substr( keyword.length() );
            result = atoi( rhs.c_str() );
        }
    }

} hardwareConstraints;

struct ResourceEstimate
{
    ResourceEstimate()
        : ramRssMBytes( 0xFFFFFFFF )
        , ramVszMBytes( 0xFFFFFFFF )
        , timeSeconds( 0xFFFFFFFFFFFFFFFF )
    {}

    bool operator==( const ResourceEstimate &rhs )
    {
        return ramRssMBytes==rhs.ramRssMBytes
               && ramVszMBytes==rhs.ramVszMBytes
               && timeSeconds == rhs.timeSeconds;
    }

    bool isSet()
    {
        return ramRssMBytes != 0xFFFFFFFF
               || ramVszMBytes != 0xFFFFFFFF
               || timeSeconds != 0xFFFFFFFFFFFFFFFF;
    }

    unsigned int ramRssMBytes;
    unsigned int ramVszMBytes;
    unsigned long timeSeconds;
};

std::ostream &operator<<( std::ostream &os, const BwtParameters &obj )
{
    os << "{";
    for ( int i=0; i<obj.size(); ++i )
        if ( obj[i] != MULTIPLE_OPTIONS )
            os << optionNames[i] << "=" << optionPossibleValues[i][obj[i]] << ", ";
        else
            os << "*, ";

    os << "}";

    return os;
}

std::ostream &operator<<( std::ostream &os, const ResourceEstimate &obj )
{
    os << "RSS: " << obj.ramRssMBytes << "MB";
    os << ",\tVSZ: " << obj.ramVszMBytes << "MB";
    os << ",\ttime: " << obj.timeSeconds << "sec";
    if ( obj.timeSeconds >= 60 )
    {
        unsigned long time = obj.timeSeconds;
        // we don't show the seconds if we show the hours (and we round to the closest minutes)
        if ( time > 3600 )
        {
            time = ( ( time+30 )/60 )*60;
        }

        int hours = ( time / 60 ) / 60;
        int minutes = ( time / 60 ) % 60;
        int seconds = time % 60;
        os << " (=";
        if ( hours )
        {
            os << hours << "h";
        }
        if ( minutes )
            os << minutes << "min";
        if ( seconds )
            os << seconds << "s";
        os << ")";
    }

    return os;
}

typedef pair< BwtParameters, ResourceEstimate > BwtParamsAndEstimate;
std::ostream &operator<<( std::ostream &os, const BwtParamsAndEstimate &obj )
{
    os << obj.second << "\tfor " << obj.first;
    return os;
}


vector< BwtParamsAndEstimate > allEstimates;

void calculateResourceRequirements( const vector< pair< string, string > > & );
ResourceEstimate estimateResourcesFor( BwtParameters &paramValues );


void recursiveFillOptionValuesAndRunEstimate( BwtParameters &paramValues )
{
    int depth = paramValues.size();
    if ( optionNames[depth] != "" )
    {
        for ( int i = 0; optionPossibleValues[depth][i] != ""; ++i )
        {
            paramValues.push_back( i );
            recursiveFillOptionValuesAndRunEstimate( paramValues );
            paramValues.pop_back();
        }
    }
    else
    {
        ResourceEstimate resourceEstimate = estimateResourcesFor( paramValues );
        if ( resourceEstimate.isSet() )
            allEstimates.push_back( make_pair( paramValues, resourceEstimate ) );
    }
}

bool sortEstimatesByTime( const BwtParamsAndEstimate &lhs, const BwtParamsAndEstimate &rhs )
{
    if ( lhs.second.timeSeconds != rhs.second.timeSeconds )
        return ( lhs.second.timeSeconds < rhs.second.timeSeconds );
    else
    {
        if ( lhs.second.ramRssMBytes != rhs.second.ramRssMBytes )
            return ( lhs.second.ramRssMBytes < rhs.second.ramRssMBytes );
        else
            return ( lhs.second.ramVszMBytes < rhs.second.ramVszMBytes );
    }
}

void filterEstimates( const vector< pair< string, string > > &filters )
{
    for ( vector< pair< string, string > >::const_iterator filter = filters.begin(); filter != filters.end(); ++filter )
    {
        Logger::out( LOG_SHOW_IF_VERBOSE ) << "Filtering with " << filter->first << " = " << filter->second << endl;
        vector< BwtParamsAndEstimate > filteredEstimates;

        // Identify resource num
        int resourceNum = 0;
        for ( resourceNum = 0; optionNames[resourceNum] != ""; ++resourceNum )
        {
            if ( strcasecmp( filter->first.c_str(), optionNames[resourceNum].c_str() ) == 0 )
            {
                break;
            }
        }
        if ( optionNames[resourceNum] == "" )
        {
            cerr << "Error: Impossible to identify resource " << filter->first << endl;
            exit( 1 );
        }

        // Identify resource value
        int resourceValue = 0;
        for ( resourceValue = 0; optionPossibleValues[resourceNum][resourceValue] != ""; ++resourceValue )
        {
            if ( strcasecmp( filter->second.c_str(), optionPossibleValues[resourceNum][resourceValue].c_str() ) == 0 )
            {
                break;
            }
        }
        if ( optionPossibleValues[resourceNum][resourceValue] == "" )
        {
            cerr << "Error: Impossible to identify resource value " << filter->second << endl;
            exit( 1 );
        }

        // Do the filtering
        // cout << " identified filter: " << optionNames[resourceNum] << " = " << optionPossibleValues[resourceNum][resourceValue] << endl;
        for ( int i=0; i<allEstimates.size(); ++i )
        {
            bool keep = ( allEstimates[i].first[resourceNum] == resourceValue );
            if ( keep )
            {
                filteredEstimates.push_back( allEstimates[i] );
            }
        }
        allEstimates.swap( filteredEstimates );
    }
}

void calculateResourceRequirements( const vector< pair< string, string > > &filters )
{
    BwtParameters paramValues;
    recursiveFillOptionValuesAndRunEstimate( paramValues );

    // Filter out estimates based on specified inputs
    filterEstimates( filters );

    // Sort by estimated processing time
    sort( allEstimates.begin(), allEstimates.end(), sortEstimatesByTime );

    // Merge identical estimates
    vector< BwtParamsAndEstimate > mergedEstimates;
    for ( int i=1; i<allEstimates.size(); ++i )
    {
        if ( allEstimates[i-1].second == allEstimates[i].second )
        {
            for ( int j=0; j<allEstimates[i-1].first.size(); ++j )
                if ( allEstimates[i-1].first[j] != allEstimates[i].first[j] )
                    allEstimates[i-1].first[j] = MULTIPLE_OPTIONS;
            allEstimates.erase( allEstimates.begin() + i );
            --i;
        }
    }

    // Print out top 10
    for ( int i=0; /*i<10 &&*/ i<allEstimates.size(); ++i )
    {
        Logger::out( LOG_SHOW_IF_VERBOSE ) << allEstimates[i] << endl;
    }
}


ResourceEstimate estimateResourcesFor( BwtParameters &paramValues )
{
    struct ResourceEstimate result;
    result.ramRssMBytes = 0;
    result.ramVszMBytes = 0;
    result.timeSeconds = 0;

    /*
        cout << "Estimating resources for ";
        for (int i=0; i<paramValues.size(); ++i)
            cout << paramValues[i] << ", ";
        cout << endl;
        OPTION_INPUT_FORMAT,
        OPTION_ALGORITHM,
        OPTION_INTERMEDIATE_FORMAT,
        OPTION_INTERMEDIATE_STORAGE_MEDIUM,
        OPTION_PARALLEL_PREFETCH,
        OPTION_PARALLEL_PROCESSING,
        OPTION_GENERATE_QUALITIES,
        OPTION_GENERATE_LCP,
        OPTION_SINGLE_CYCLE,
        OPTION_COUNT // end marker
    */

    //    unsigned long dataSize = static_cast<unsigned long>( datasetMetadata.nBases * datasetMetadata.rleCompressibility / 8 );
    unsigned long dataRead = 0;
    unsigned long dataWritten = 0;
    for ( int i=0; i<datasetMetadata.nCycles; ++i )
    {
        switch ( paramValues[OPTION_INTERMEDIATE_FORMAT] )
        {
            case INTERMEDIATE_FORMAT_ASCII:
            {
                dataRead += static_cast<unsigned long>( datasetMetadata.nReads * i );
                dataWritten += static_cast<unsigned long>( datasetMetadata.nReads * i );
                break;
            }
            case INTERMEDIATE_FORMAT_HUFFMAN:
                // unknown => special value so that it only gets selected if specified by the user
                result.timeSeconds = 10*60*60;
            case INTERMEDIATE_FORMAT_RLE:
            {
                float cycleRleCompressibility = 2 * datasetMetadata.rleCompressibility - datasetMetadata.rleCompressibility*i/datasetMetadata.nCycles;
                dataRead += static_cast<unsigned long>( datasetMetadata.nReads * i * cycleRleCompressibility / 8 );
                dataWritten += static_cast<unsigned long>( datasetMetadata.nReads * i * cycleRleCompressibility / 8 );
                break;
            }
            case INTERMEDIATE_FORMAT_MULTIRLE:
            {
                float cycleRleCompressibility = 2 * datasetMetadata.rleCompressibility - datasetMetadata.rleCompressibility*i/datasetMetadata.nCycles;
                dataRead += static_cast<unsigned long>( datasetMetadata.nReads * i * cycleRleCompressibility / 8 );
                dataWritten += static_cast<unsigned long>( 4 * datasetMetadata.nReads * cycleRleCompressibility / 8 );
                break;
            }
            default:
                assert( 0 && "should never reach here" );
        }
    }

    switch ( paramValues[OPTION_ALGORITHM] )
    {
        case ALGORITHM_BCR:
            /*
                    // Beetl bcr -r
                OPTION_INPUT_FORMAT                     = not taken into consideration
                OPTION_ALGORITHM                        = bcr
                OPTION_INTERMEDIATE_FORMAT         = RLE
                OPTION_INTERMEDIATE_STORAGE_MEDIUM = disk
                OPTION_PARALLEL_PREFETCH                = OFF
                OPTION_PARALLEL_PROCESSING              = OFF
                OPTION_GENERATE_QUALITIES               = n/a
                OPTION_GENERATE_LCP                     = n/a
                OPTION_SINGLE_CYCLE                     = n/a

                    // Beetl bcr -r + parallel (branch027)
                OPTION_INPUT_FORMAT                     = not taken into consideration
                OPTION_ALGORITHM                        = bcr
                OPTION_INTERMEDIATE_FORMAT         = RLE
                OPTION_INTERMEDIATE_STORAGE_MEDIUM = disk
                OPTION_PARALLEL_PREFETCH                = ON
                OPTION_PARALLEL_PROCESSING              = OFF/2cores/4cores
                OPTION_GENERATE_QUALITIES               = n/a
                OPTION_GENERATE_LCP                     = n/a
                OPTION_SINGLE_CYCLE                     = n/a

                    // Beetl bcr -t (parallel)
                OPTION_INPUT_FORMAT                     = not taken into consideration
                OPTION_ALGORITHM                        = bcr
                OPTION_INTERMEDIATE_FORMAT         = MultiRLE
                OPTION_INTERMEDIATE_STORAGE_MEDIUM = RAM
                OPTION_PARALLEL_PREFETCH                = ON
                OPTION_PARALLEL_PROCESSING              = OFF/2cores/4cores
                OPTION_GENERATE_QUALITIES               = n/a
                OPTION_GENERATE_LCP                     = n/a
                OPTION_SINGLE_CYCLE                     = n/a
            */
        {
            float diskSpeed = min( hardwareConstraints.diskMaxSpeedPerProcess * ( 1 << paramValues[OPTION_PARALLEL_PROCESSING] )
                                   , hardwareConstraints.diskMaxTotalSpeed );
            float ramSpeed = min( hardwareConstraints.ramMaxSpeedPerProcess * ( 1 << paramValues[OPTION_PARALLEL_PROCESSING] )
                                  , hardwareConstraints.ramMaxTotalSpeed );

            // base RAM
            result.ramRssMBytes += 14 * datasetMetadata.nReads / 1024 / 1024 * ( paramValues[OPTION_PARALLEL_PROCESSING]?2:1 );

            // TransposeFasta
            switch ( paramValues[OPTION_INPUT_FORMAT] )
            {
                case INPUT_FORMAT_FASTQ:
                    dataRead += static_cast<unsigned long>( 2 * datasetMetadata.nReads * datasetMetadata.nCycles );
                    dataWritten += static_cast<unsigned long>( datasetMetadata.nReads * datasetMetadata.nCycles );
                    break;
                case INPUT_FORMAT_FASTA:
                case INPUT_FORMAT_SEQ:
                    dataRead += static_cast<unsigned long>( datasetMetadata.nReads * datasetMetadata.nCycles );
                    dataWritten += static_cast<unsigned long>( datasetMetadata.nReads * datasetMetadata.nCycles );
                    break;
                case INPUT_FORMAT_CYC:
                case INPUT_FORMAT_BCL:
                    break;
            }

            unsigned long dataRw = dataRead + dataWritten;

            if ( paramValues[OPTION_PARALLEL_PREFETCH] == PARALLEL_PREFETCH_ON )
            {
                unsigned int prefetchRam = datasetMetadata.nReads / 1024 / 1024;
                result.ramRssMBytes += prefetchRam;
            }
            else
            {
                unsigned long prefetchTime = static_cast<unsigned long>( datasetMetadata.nBases /  hardwareConstraints.diskMaxSpeedPerProcess / 1024 / 1024 ) + 1;
                result.timeSeconds += prefetchTime;
            }

            if ( paramValues[OPTION_INTERMEDIATE_STORAGE_MEDIUM] == INTERMEDIATE_STORAGE_MEDIUM_DISK )
            {
                unsigned long dataRwTime = static_cast<unsigned long>( dataRw / diskSpeed / 1024 / 1024 ) + 16;
                result.timeSeconds += dataRwTime;
            }
            else // == INTERMEDIATE_STORAGE_MEDIUM_RAM
            {
                unsigned long dataRwTime = static_cast<unsigned long>( dataRw / ramSpeed / 1024 / 1024 ) + 16;
                result.timeSeconds += dataRwTime;

                result.ramRssMBytes += 14 * datasetMetadata.nReads / 1024 / 1024; //?
                result.ramRssMBytes += static_cast<unsigned long>( datasetMetadata.nReads * datasetMetadata.nCycles * datasetMetadata.rleCompressibility / 8 * 2 / 1024 / 1024 );
            }
            break;
        }

        case ALGORITHM_EXT:
        {
            // Beetl ext
            result.ramRssMBytes = 2;
            result.ramVszMBytes = 19;

            // TransposeFasta
            switch ( paramValues[OPTION_INPUT_FORMAT] )
            {
                case INPUT_FORMAT_CYC:
                case INPUT_FORMAT_BCL:
                    dataRead += static_cast<unsigned long>( datasetMetadata.nReads * datasetMetadata.nCycles );
                    dataWritten += static_cast<unsigned long>( datasetMetadata.nReads * datasetMetadata.nCycles );
                    break;
                case INPUT_FORMAT_FASTQ:
                case INPUT_FORMAT_FASTA:
                case INPUT_FORMAT_SEQ:
                    break;
            }

            // Main EXT formula
            unsigned long dataRw = dataRead + dataWritten;
            unsigned long dataSfiles_RW = static_cast<unsigned long>( datasetMetadata.nReads * datasetMetadata.nCycles );
            unsigned long dataPfiles_RW = static_cast<unsigned long>( datasetMetadata.nReads * 8 * datasetMetadata.nCycles );
            result.timeSeconds = static_cast<unsigned long>( ( dataRw + dataSfiles_RW + dataPfiles_RW ) / hardwareConstraints.diskMaxSpeedPerProcess / 1024 / 1024 ) + 19;
            break;
        }

        default:
            cerr << "Error: should never reach here" << endl;
            exit ( -1 );
    }

    Logger::out( LOG_FOR_DEBUGGING ) << paramValues << " => " << result << endl;

    return result;
}


void launchBeetlBwt( BwtParamsAndEstimate &config, const string &inputFilename, const string &outputFilename )
{
    Logger::out( LOG_ALWAYS_SHOW ) << "\nLaunching the following configuration of Beetl:" << endl;
    for ( int i=0; i<config.first.size(); ++i )
    {
        Logger::out( LOG_ALWAYS_SHOW ) << "  " << optionNames[i] << " = ";
        if ( config.first[i] == MULTIPLE_OPTIONS )
        {
            Logger::out( LOG_ALWAYS_SHOW ) << "* => ";
            config.first[i] = 0;
        }
        Logger::out( LOG_ALWAYS_SHOW ) << optionPossibleValues[i][config.first[i]] << endl;
    }

    if ( config.first[OPTION_ALGORITHM] == ALGORITHM_BCR )
    {
        int bcrMode = 0 ; // 0=build BWT
        CompressionFormatType outputCompression = compressionIncrementalRunLength; // todo
        /*
          compressionASCII,
          compressionRunLength,
          compressionIncrementalRunLength,
          compressionHuffman
        */
        BCRexternalBWT bwt( ( char * )inputFilename.c_str(), ( char * )outputFilename.c_str(), bcrMode, outputCompression, &config.first );
        return;
    }

    stringstream cmdLineParams;
    cmdLineParams << optionPossibleValues[OPTION_ALGORITHM][config.first[OPTION_ALGORITHM]]
                  << " -i " << inputFilename;

    switch ( config.first[OPTION_ALGORITHM] )
    {
        case MULTIPLE_OPTIONS:
        case ALGORITHM_BCR:
            cmdLineParams << " -o " << outputFilename;
            cmdLineParams << " -m 0 ";
            break;
        case ALGORITHM_EXT:
            cmdLineParams << " -p " << outputFilename;
            break;
        default:
            assert( false && "shouldn't reach here" );
    }

    if ( config.first[OPTION_INTERMEDIATE_FORMAT] == INTERMEDIATE_FORMAT_MULTIRLE )
    {
        cmdLineParams << " -t ";
    }
    else switch ( config.first[OPTION_OUTPUT_FORMAT] )
        {
            case MULTIPLE_OPTIONS:
            case OUTPUT_FORMAT_ASCII:
                cmdLineParams << " -a ";
                break;
            case OUTPUT_FORMAT_RLE:
                cmdLineParams << " -r ";
                break;
            case OUTPUT_FORMAT_HUFFMAN:
                cmdLineParams << " -h ";
                break;
            default:
                assert( false && "shouldn't reach here" );
        }

    if ( config.first[OPTION_SAP_ORDERING] == SAP_ORDERING_ON )
    {
        cmdLineParams << " -sap ";
    }

    launchBeetl( cmdLineParams.str() );
}


void printUsage()
{
    cout << "Parameters:" << endl;
    cout << "    --input (-i)             Input file name or prefix" << endl;
    cout << "    --output (-o)            Output file name or prefix" << endl;
    cout << endl;
    cout << "Automatic parameters:" << endl;
    cout << "    --input-format           = autodetect [fasta|fastq|cyc|seq|bcl]" << endl;
    cout << "    --output-format (-f)     = rle [ascii|rle|huffman]" << endl;
    cout << "    --intermediate-format    = auto [ascii|rle|multirle|huffman]" << endl;
    cout << "    --intermediate-medium    = auto [disk|ram] (multirle->ram, others->disk)" << endl;
    cout << "    --algorithm (-a)         = auto [bcr|ext]" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "    --qualities (-q)         = ignore [ignore|permute]: Ignore/Permute qualities (only available with bcr algorithm)" << endl;
    cout << "    --concatenate-output     Concatenate BWT files at the end" << endl;
    //TODO in future release:    cout << "    --generate-lcp           " << endl;
    cout << "    --sap-ordering           Use SAP ordering (see SAP note below)" << endl;
    cout << "    --generate-end-pos-file  Generates outFileEndPos.bwt" << endl;
    //TODO in future release:    cout << "    --single-cycle           " << endl;
#ifdef USE_OPENMP
    cout << "    --no-parallel-prefetch   Disable parallel prefetch of cycle files" << endl;
    cout << "    --no-parallel-processing Disable parallel processing by letter" << endl;
#endif //ifdef USE_OPENMP
    cout << "    --hw-constraints         File describing hardware constraints for speed estimates" << endl;
    cout << "    --verbose                = quiet [quiet|verbose|very-verbose|debug] or [0|1|2|3]" << endl;
    cout << "    --help (-h)              Help" << endl;
    cout << endl;
    cout << "Notes:" << endl;
    cout << "    RLE      : run-length-encoded format" << endl;
    cout << "    multiRLE : run-length-encoded using an incremental strategy with multiple files" << endl;
    cout << "    SAP      : implicit permutation to obtain more compressible BWT (Note: forces algorithm=ext and intermediate-format=ascii)" << endl;
#ifndef USE_OPENMP
    cout << endl;
    cout << "Warning:" << endl;
    cout << "    This version of BEETL hasn't been compiled with OpenMP parallelisation. Performance may be degraded." << endl;
#endif //ifndef USE_OPENMP
    cout << endl;
}

struct BeetlBwtArguments
{
    string argInput;
    string argOutput             ;// = "outBWT"
    string argInputFormat;
    string argOutputFormat       ;// = "rle"
    bool   argConcatenateOutput  ;//= false;
    string argAlgorithm          ;//= "auto";
    string argIntermediateFormat ;//= "auto";
    string argIntermediateMedium ;//= "auto";
    bool   argParallelPrefetch   ;//= true;
    bool   argParallelProcessing ;//= true;
    string argQualities          ;//= "ignore";
    bool   argGenerateLcp        ;//= false;
    bool   argSapOrdering        ;//= false;
    bool   argGenerateEndPosFile ;//= false;
    bool   argSingleCycle        ;//= false;
    string argHardwareConstraints;
    string argVerbosityLevel     ;//= "0"

    BeetlBwtArguments()
    {
        argOutput             = "outBWT";
        argOutputFormat       = "rle";
        argConcatenateOutput  = false;
        argAlgorithm          = "";//"auto";
        argIntermediateFormat = "";//"auto";
        argIntermediateMedium = "";//"auto";
        argParallelPrefetch   = true;
        argParallelProcessing = true;
        argQualities          = "ignore";
        argGenerateLcp        = false;
        argSapOrdering        = false;
        argGenerateEndPosFile = false;
        argSingleCycle        = false;
        argVerbosityLevel     = "0";
    }

    int getMaxRamInMB()
    {
        return 0;
    }
};

int main( const int argc, const char **argv )
{
    // Generated using: http://patorjk.com/software/taag/#p=display&f=Soft&t=BEETL%20BWT
    cout << ",-----.  ,------.,------.,--------.,--.       ,-----.  ,--.   ,--.,--------. " << endl;
    cout << "|  |) /_ |  .---'|  .---''--.  .--'|  |       |  |) /_ |  |   |  |'--.  .--' " << endl;
    cout << "|  .-.  \\|  `--, |  `--,    |  |   |  |       |  .-.  \\|  |.'.|  |   |  |    " << endl;
    cout << "|  '--' /|  `---.|  `---.   |  |   |  '--.    |  '--' /|   ,'.   |   |  |    " << endl;
    cout << "`------' `------'`------'   `--'   `-----'    `------' '--'   '--'   `--'    " << endl;
    cout << "Version " << PACKAGE_VERSION << endl;
    cout << endl;

    cout << "Command called:" << endl << "   ";
    for ( int i=0; i < argc; ++i )
    {
        cout << " " << argv[i];
    }
    cout << "\n" << endl;

    BeetlBwtArguments args;

    for ( int i = 1; i < argc; ++i )
    {
        if ( isNextArgument( "-h", "--help", argc, argv, i ) )
        {
            printUsage();
            exit( 0 );
        }
        else if ( isNextArgument( "-i", "--input"                  , argc, argv, i, &args.argInput               ) ) {}
        else if ( isNextArgument( "-o", "--output"                 , argc, argv, i, &args.argOutput              ) ) {}
        else if ( isNextArgument( "-f", "--output-format"          , argc, argv, i, &args.argOutputFormat        ) ) {}
        else if ( isNextArgument( ""  , "--input-format"           , argc, argv, i, &args.argInputFormat         ) ) {}
        else if ( isNextArgument( ""  , "--intermediate-format"    , argc, argv, i, &args.argIntermediateFormat  ) ) {}
        else if ( isNextArgument( ""  , "--intermediate-medium"    , argc, argv, i, &args.argIntermediateMedium  ) ) {}
        else if ( isNextArgument( "-a", "--algorithm"              , argc, argv, i, &args.argAlgorithm           ) ) {}
        else if ( isNextArgument( "-q", "--qualities"              , argc, argv, i, &args.argQualities           ) ) {}
        else if ( isNextArgument( ""  , "--concatenate-output"     , argc, argv, i                               ) )
        {
            args.argConcatenateOutput = true;
        }
        else if ( isNextArgument( ""  , "--sap-ordering"           , argc, argv, i                               ) )
        {
            args.argSapOrdering = true;
        }
        else if ( isNextArgument( ""  , "--generate-end-pos-file"  , argc, argv, i                               ) )
        {
            args.argGenerateEndPosFile = true;
        }
        else if ( isNextArgument( ""  , "--no-parallel-prefetch"   , argc, argv, i                               ) )
        {
            args.argParallelPrefetch = false;
        }
        else if ( isNextArgument( ""  , "--no-parallel-processing" , argc, argv, i                               ) )
        {
            args.argParallelProcessing = false;
        }
        //        else if (isNextArgument( ""  , "--generate-lcp"           , argc, argv, i                               )) { args.argGenerateLcp = true; }
        //        else if (isNextArgument( ""  , "--single-cycle"           , argc, argv, i                               )) { args.argSingleCycle = true; }
        else if ( isNextArgument( ""  , "--hw-constraints"         , argc, argv, i, &args.argHardwareConstraints ) ) {}
        else if ( isNextArgument( ""  , "--verbose"                , argc, argv, i, &args.argVerbosityLevel      ) )
        {
            Logger::setVerbosity( args.argVerbosityLevel );
        }
        else
        {
            cerr << "Error: Invalid parameter: " << argv[i] << "\n" << endl;
            printUsage();
            exit( 1 );
        }
    }

    // Checking for required parameters
    if ( args.argInput.empty() || args.argOutput.empty() )
    {
        cerr << "Error: Missing arguments: --input is required\n" << endl;
        printUsage();
        exit( 1 );
    }

    // Auto-detection of missing arguments
    if ( args.argInputFormat.empty() )
    {
        args.argInputFormat = detectFileFormat( args.argInput );
    }
    checkFileFormat( args.argInput, args.argInputFormat );

    // Special case of SAP ordering
    if ( args.argSapOrdering )
    {
        if ( args.argAlgorithm != "" && strcasecmp( args.argAlgorithm.c_str(), "ext" ) != 0 )
        {
            clog << "Warning: Forcing algorithm=ext for --sap-ordering" << endl;
            args.argAlgorithm = "ext";
        }
        if ( args.argIntermediateFormat != "" && strcasecmp( args.argIntermediateFormat.c_str(), "ascii" ) != 0 )
        {
            clog << "Warning: Forcing intermediate-format=ascii for --sap-ordering" << endl;
            args.argIntermediateFormat = "ascii";
        }
    }

    datasetMetadata.init( args.argInput, args.argInputFormat );
    hardwareConstraints.init( args.argHardwareConstraints );


    // Preparation of filters for resource estimation
    vector< pair< string, string > > filters;
    if ( !args.argInputFormat.empty() )
        filters.push_back( make_pair( "input format", args.argInputFormat ) );
    if ( !args.argOutputFormat.empty() )
        filters.push_back( make_pair( "output format", args.argOutputFormat ) );
    if ( !args.argAlgorithm.empty() )
        filters.push_back( make_pair( "algorithm", args.argAlgorithm ) );
    if ( !args.argIntermediateFormat.empty() )
        filters.push_back( make_pair( "intermediate format", args.argIntermediateFormat ) );
    if ( !args.argIntermediateMedium.empty() )
        filters.push_back( make_pair( "intermediate storage medium", args.argIntermediateMedium ) );
    filters.push_back( make_pair( "concatenate output", args.argConcatenateOutput?"on":"off" ) );
    filters.push_back( make_pair( "SAP ordering", args.argSapOrdering?"on":"off" ) );
    filters.push_back( make_pair( "generate endPosfile", args.argGenerateEndPosFile?"on":"off" ) );
    if ( args.argQualities == "ignore" )
        filters.push_back( make_pair( "permute qualities", "off" ) );
    else if ( args.argQualities == "permute" )
    {
        filters.push_back( make_pair( "permute qualities", "on" ) );
        if ( strcasecmp( args.argAlgorithm.c_str(), "ext" ) == 0 )
        {
            cerr << "Error: Quality permutation is not available with \"ext\" algorithm\n" << endl;
            printUsage();
            exit( 1 );
        }
    }
    else
    {
        cerr << "Error: Invalid --qualities=" << args.argQualities << " value\n" << endl;
        printUsage();
        exit( 1 );
    }

    // Resource estimation
    calculateResourceRequirements( filters );

    // Launch the fastest estimated configuration of Beetl
    if ( allEstimates.empty() )
    {
        cerr << "Error: Beetl doesn't seem to be able to run under these constraints. Try to relax some of them." << endl;
        exit( 2 );
    }
    launchBeetlBwt( allEstimates[0], args.argInput, args.argOutput );

    return 0;
}
