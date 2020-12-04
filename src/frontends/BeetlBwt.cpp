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

#include "BeetlBwt.hh"

#include "BCRext.hh"
#include "BCRexternalBWT.hh"
#include "Common.hh"
#include "DatasetMetadata.hh"
#include "parameters/BwtParameters.hh"
#include "config.h"
#include "libzoo/cli/Common.hh"
#include "libzoo/util/Logger.hh"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;
using namespace BeetlBwtParameters;


BwtParameters params;



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
    {
        unset();
    }

    bool operator==( const ResourceEstimate &rhs ) const
    {
        return ramRssMBytes == rhs.ramRssMBytes
               && ramVszMBytes == rhs.ramVszMBytes
               && timeSeconds == rhs.timeSeconds;
    }

    bool isSet() const
    {
        return ramRssMBytes != 0xFFFFFFFF
               || ramVszMBytes != 0xFFFFFFFF
               || timeSeconds != 0xFFFFFFFFFFFFFFFFull;
    }

    void unset()
    {
        ramRssMBytes = 0xFFFFFFFF;
        ramVszMBytes = 0xFFFFFFFF;
        timeSeconds = 0xFFFFFFFFFFFFFFFFull;
    }

    uint32_t ramRssMBytes;
    uint32_t ramVszMBytes;
    uint64_t timeSeconds;
};

typedef pair< BwtParameters, ResourceEstimate > BwtParamsAndEstimate;
vector< BwtParamsAndEstimate > allEstimates;

std::ostream &operator<<( std::ostream &os, const ResourceEstimate &obj )
{
    if ( !obj.isSet() )
    {
        os << "<undefined>";
        return os;
    }

    os << "RAM_RSS: " << obj.ramRssMBytes << "MB";
    os << ",\tRAM_VSZ: " << obj.ramVszMBytes << "MB";
    os << ",\ttime: " << obj.timeSeconds << "sec";
    if ( obj.timeSeconds >= 60 )
    {
        uint64_t time = obj.timeSeconds;
        // we don't show the seconds if we show the hours (and we round to the closest minutes)
        if ( time > 3600 )
        {
            time = ( ( time + 30 ) / 60 ) * 60;
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

std::ostream &operator<<( std::ostream &os, const BwtParamsAndEstimate &obj )
{
    os << obj.second << "\tfor ";
    obj.first.print( os, true, AUTOMATED );
    return os;
}


ResourceEstimate estimateResourcesFor( BwtParameters &paramValues );


/*
void recursiveFillOptionValuesAndRunEstimate( BwtParameters &paramValues )
{
        int depth = paramValues.size();
        if ( paramValues.getOptionName( depth ) != "" )
        {
            for ( int i = 0; paramValues.getOptionPossibleValue( depth, i ) != ""; ++i )
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
*/

void recursiveGenerateInputCombinationsAndRunEstimate( BwtParameters &intermediateParams, int unsigned depth )
{
    if ( depth < intermediateParams.entries_.size() )
    {
        if ( !( intermediateParams.entries_[depth].flags & AUTOMATED ) || intermediateParams.entries_[depth].isSet() )
        {
            recursiveGenerateInputCombinationsAndRunEstimate( intermediateParams, depth + 1 );
        }
        else
        {
            static const string switchPossibleValues[] = { "0", "1", "" };
            const string *possibleValues = ( intermediateParams.entries_[depth].flags & TYPE_SWITCH ) ? switchPossibleValues : intermediateParams.entries_[depth].possibleValues;
            for ( int i = 0; possibleValues[i] != ""; ++i )
            {
                BwtParameters newParams = intermediateParams;
                newParams.entries_[depth].silentSet( i );
                recursiveGenerateInputCombinationsAndRunEstimate( newParams, depth + 1 );
            }
        }
    }
    else
    {
        ResourceEstimate resourceEstimate = estimateResourcesFor( intermediateParams );
        if ( resourceEstimate.isSet() )
            allEstimates.push_back( make_pair( intermediateParams, resourceEstimate ) );
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

void filterEstimates( const vector< pair< string, string > > &filters, BwtParameters &paramValues )
{
    /*
        for ( vector< pair< string, string > >::const_iterator filter = filters.begin(); filter != filters.end(); ++filter )
        {
            Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << "Filtering with " << filter->first << " = " << filter->second << endl;
            vector< BwtParamsAndEstimate > filteredEstimates;

            // Identify resource num
            int resourceNum = 0;
            for ( resourceNum = 0; paramValues.getOptionName( resourceNum ) != ""; ++resourceNum )
            {
                if ( strcasecmp( filter->first.c_str(), paramValues.getOptionName( resourceNum ).c_str() ) == 0 )
                {
                    break;
                }
            }
            if ( paramValues.getOptionName( resourceNum ) == "" )
            {
                cerr << "Error: Impossible to identify resource " << filter->first << endl;
                exit( 1 );
            }

            // Identify resource value
            int resourceValue = 0;
            for ( resourceValue = 0; paramValues.getOptionPossibleValue( resourceNum, resourceValue ) != ""; ++resourceValue )
            {
                if ( strcasecmp( filter->second.c_str(), paramValues.getOptionPossibleValue( resourceNum, resourceValue ).c_str() ) == 0 )
                {
                    break;
                }
            }
            if ( paramValues.getOptionPossibleValue( resourceNum, resourceValue ) == "" )
            {
                cerr << "Error: Impossible to identify resource value " << filter->second << endl;
                exit( 1 );
            }

            // Do the filtering
            // cout << " identified filter: " << paramValues.getOptionName(resourceNum) << " = " << paramValues.getOptionPossibleValue(resourceNum,resourceValue) << endl;
            for ( unsigned int i = 0; i < allEstimates.size(); ++i )
            {
                bool keep = ( allEstimates[i].first[resourceNum] == resourceValue );
                if ( keep )
                {
                    filteredEstimates.push_back( allEstimates[i] );
                }
            }
            allEstimates.swap( filteredEstimates );
        }
    */
}

void calculateResourceRequirements()
{
    const unsigned int memoryLimitMB = params["memory limit MB"];

    recursiveGenerateInputCombinationsAndRunEstimate( params, 0 );

    // Sort by estimated processing time
    sort( allEstimates.begin(), allEstimates.end(), sortEstimatesByTime );

    /*
        // Merge identical estimates
        for ( unsigned int i = 1; i < allEstimates.size(); ++i )
        {
            if ( allEstimates[i - 1].second == allEstimates[i].second )
            {
                allEstimates[i - 1].first.mergeWith( allEstimates[i].first );
                allEstimates.erase( allEstimates.begin() + i );
                --i;
            }
        }
    */

    // Filter out estimates above the memory limit
    Logger::out() << "RAM constraint: " << memoryLimitMB << " MB" << endl;
    bool filteringOutTopPicks = true;
    for ( unsigned int i = 0; i < allEstimates.size(); ++i )
    {
        if ( allEstimates[i].second.ramRssMBytes > memoryLimitMB )
        {
            Logger_if( filteringOutTopPicks ? LOG_SHOW_IF_VERBOSE : LOG_FOR_DEBUGGING ) Logger::out() << "RAM too low for: " << allEstimates[i] << endl;
            allEstimates.erase( allEstimates.begin() + i );
            --i;
        }
        else
            filteringOutTopPicks = false;
    }

    // Print out /*top 10*/ estimates
    for ( unsigned int i = 0; /*i<10 &&*/ i < allEstimates.size(); ++i )
    {
        Logger_if( LOG_SHOW_IF_VERBOSE ) Logger::out() << allEstimates[i] << endl;
    }
}


ResourceEstimate estimateResourcesFor( BwtParameters &paramValues )
{
    struct ResourceEstimate result;
    result.ramRssMBytes = 0;
    result.ramVszMBytes = 0;
    result.timeSeconds = 0;

    //        cout << "Estimating resources for ";
    //        for (int i=0; i<paramValues.size(); ++i)
    //            cout << paramValues[i] << ", ";
    //        cout << endl;

    uint64_t nReads = datasetMetadata.nReads;
    if ( paramValues[ PARAMETER_ADD_REV_COMP ] == 1 )
        nReads *= 2;    // Reverse-complemented reads

    //    uint64_t dataSize = static_cast<uint64_t>( datasetMetadata.nBases * datasetMetadata.rleCompressibility / 8 );
    uint64_t dataRead = 0;
    uint64_t dataWritten = 0;
    for ( unsigned int i = 0; i < datasetMetadata.nCycles; ++i )
    {
        switch ( paramValues[PARAMETER_INTERMEDIATE_FORMAT] )
        {
            case INTERMEDIATE_FORMAT_ASCII:
            {
                dataRead += static_cast<uint64_t>( nReads * i );
                dataWritten += static_cast<uint64_t>( nReads * i );
                break;
            }
            case INTERMEDIATE_FORMAT_HUFFMAN:
                // unknown => special value so that it only gets selected if specified by the user
                result.timeSeconds = 10 * 60 * 60;
            case INTERMEDIATE_FORMAT_RLE:
            {
                float cycleRleCompressibility = 2 * datasetMetadata.rleCompressibility - datasetMetadata.rleCompressibility * i / datasetMetadata.nCycles;
                dataRead += static_cast<uint64_t>( nReads * i * cycleRleCompressibility / 8 );
                dataWritten += static_cast<uint64_t>( nReads * i * cycleRleCompressibility / 8 );
                break;
            }
            case INTERMEDIATE_FORMAT_MULTIRLE:
            {
                float cycleRleCompressibility = 2 * datasetMetadata.rleCompressibility - datasetMetadata.rleCompressibility * i / datasetMetadata.nCycles;
                dataRead += static_cast<uint64_t>( nReads * i * cycleRleCompressibility / 8 );
                dataWritten += static_cast<uint64_t>( 4 * nReads * cycleRleCompressibility / 8 );
                break;
            }
            default:
                assert( 0 && "should never reach here" );
        }
    }

    // Deactivate MultiRLE by estimating it slower (special value so that it only gets selected if specified by the user)
    if ( paramValues[PARAMETER_INTERMEDIATE_FORMAT] == INTERMEDIATE_FORMAT_MULTIRLE )
    {
        result.timeSeconds += 100 * 60 * 60;
    }

    switch ( paramValues[PARAMETER_ALGORITHM] )
    {
        case ALGORITHM_BCR:
            /*
                    // Beetl bcr -r
                PARAMETER_INPUT_FORMAT                     = not taken into consideration
                PARAMETER_ALGORITHM                        = bcr
                PARAMETER_INTERMEDIATE_FORMAT         = RLE
                PARAMETER_INTERMEDIATE_STORAGE_MEDIUM = disk
                PARAMETER_PARALLEL_PREFETCH                = OFF
                PARAMETER_PARALLEL_PROCESSING              = OFF
                PARAMETER_GENERATE_QUALITIES               = n/a
                PARAMETER_GENERATE_LCP                     = n/a
                PARAMETER_SINGLE_CYCLE                     = n/a

                    // Beetl bcr -r + parallel (branch027)
                PARAMETER_INPUT_FORMAT                     = not taken into consideration
                PARAMETER_ALGORITHM                        = bcr
                PARAMETER_INTERMEDIATE_FORMAT         = RLE
                PARAMETER_INTERMEDIATE_STORAGE_MEDIUM = disk
                PARAMETER_PARALLEL_PREFETCH                = ON
                PARAMETER_PARALLEL_PROCESSING              = OFF/2cores/4cores
                PARAMETER_GENERATE_QUALITIES               = n/a
                PARAMETER_GENERATE_LCP                     = n/a
                PARAMETER_SINGLE_CYCLE                     = n/a

                    // Beetl bcr -t (parallel)
                PARAMETER_INPUT_FORMAT                     = not taken into consideration
                PARAMETER_ALGORITHM                        = bcr
                PARAMETER_INTERMEDIATE_FORMAT         = MultiRLE
                PARAMETER_INTERMEDIATE_STORAGE_MEDIUM = RAM
                PARAMETER_PARALLEL_PREFETCH                = ON
                PARAMETER_PARALLEL_PROCESSING              = OFF/2cores/4cores
                PARAMETER_GENERATE_QUALITIES               = n/a
                PARAMETER_GENERATE_LCP                     = n/a
                PARAMETER_SINGLE_CYCLE                     = n/a
            */
        {
            float diskSpeed = min( hardwareConstraints.diskMaxSpeedPerProcess * ( 4/*1 << paramValues[PARAMETER_PARALLEL_PROCESSING]*/ )
                                   , hardwareConstraints.diskMaxTotalSpeed );
            float ramSpeed = min( hardwareConstraints.ramMaxSpeedPerProcess * ( 4/*1 << paramValues[PARAMETER_PARALLEL_PROCESSING]*/ )
                                  , hardwareConstraints.ramMaxTotalSpeed );

            // base RAM
            result.ramRssMBytes += 14 * nReads / 1024 / 1024; // main structures
            if ( true /*paramValues[PARAMETER_PARALLEL_PROCESSING]*/ )
            {
                result.ramRssMBytes += sizeof( struct sortElement ) * nReads / 1024 / 1024; // vectTriple duplicated for parallel action
                // buffer due to growing vectors: worst case is twice the data size of a well distributed data between the 16 buffers (multiplied by 4 as a precaution), max 1GB for each of the 16 buffers
                uint32_t bufferSize = min<uint32_t>( nReads / 16 * 4 * sizeof( struct sortElement ) / 1024 / 1024 * 2 / 16, 1024 );
                result.ramRssMBytes += 16 * bufferSize;
            }

            // TransposeFasta
            switch ( paramValues[PARAMETER_INPUT_FORMAT] )
            {
                case INPUT_FORMAT_FASTQ:
                    dataRead += static_cast<uint64_t>( 2 * nReads * datasetMetadata.nCycles );
                    dataWritten += static_cast<uint64_t>( nReads * datasetMetadata.nCycles );
                    break;
                case INPUT_FORMAT_FASTA:
                case INPUT_FORMAT_SEQ:
                    dataRead += static_cast<uint64_t>( nReads * datasetMetadata.nCycles );
                    dataWritten += static_cast<uint64_t>( nReads * datasetMetadata.nCycles );
                    break;
                case INPUT_FORMAT_CYC:
                case INPUT_FORMAT_BCL:
                    break;
            }

            uint64_t dataRw = dataRead + dataWritten;

            //#ifdef _OPENMP
            //always activated at the moment            if ( paramValues[PARAMETER_PARALLEL_PREFETCH] == true )
            {
                unsigned int prefetchRam = nReads / 1024 / 1024;
                result.ramRssMBytes += prefetchRam;
            }
            //else
            //#endif // ifdef _OPENMP
            {
                uint64_t prefetchTime = static_cast<uint64_t>( datasetMetadata.nBases /  hardwareConstraints.diskMaxSpeedPerProcess / 1024 / 1024 ) + 1;
                result.timeSeconds += prefetchTime;
            }

            if ( true /*paramValues[PARAMETER_INTERMEDIATE_STORAGE_MEDIUM] == INTERMEDIATE_STORAGE_MEDIUM_DISK*/ )
            {
                uint64_t dataRwTime = static_cast<uint64_t>( dataRw / diskSpeed / 1024 / 1024 ) + 16;
                result.timeSeconds += dataRwTime;
            }
            else // == INTERMEDIATE_STORAGE_MEDIUM_RAM
            {
                uint64_t dataRwTime = static_cast<uint64_t>( dataRw / ramSpeed / 1024 / 1024 ) + 16;
                result.timeSeconds += dataRwTime;

                result.ramRssMBytes += 14 * nReads / 1024 / 1024; //?
                result.ramRssMBytes += static_cast<uint64_t>( nReads * datasetMetadata.nCycles * datasetMetadata.rleCompressibility / 8 * 2 / 1024 / 1024 );
            }
            break;
        }

        case ALGORITHM_EXT:
        {
            // Beetl ext

            // Filter out unsupported options
            if ( paramValues[PARAMETER_INTERMEDIATE_FORMAT] == INTERMEDIATE_FORMAT_MULTIRLE )
            {
                result.unset();
                break;
            }

            // Fixed RAM consumption for all Beetl ext
            result.ramRssMBytes = 2;
            result.ramVszMBytes = 19;

            // TransposeFasta
            switch ( paramValues[PARAMETER_INPUT_FORMAT] )
            {
                case INPUT_FORMAT_CYC:
                case INPUT_FORMAT_BCL:
                    dataRead += static_cast<uint64_t>( nReads * datasetMetadata.nCycles );
                    dataWritten += static_cast<uint64_t>( nReads * datasetMetadata.nCycles );
                    break;
                case INPUT_FORMAT_FASTQ:
                case INPUT_FORMAT_FASTA:
                case INPUT_FORMAT_SEQ:
                    break;
            }

            // Main EXT formula
            uint64_t dataRw = dataRead + dataWritten;
            uint64_t dataSfiles_RW = static_cast<uint64_t>( nReads * datasetMetadata.nCycles );
            uint64_t dataPfiles_RW = static_cast<uint64_t>( nReads * 8 * datasetMetadata.nCycles );
            result.timeSeconds = static_cast<uint64_t>( ( dataRw + dataSfiles_RW + dataPfiles_RW ) / hardwareConstraints.diskMaxSpeedPerProcess / 1024 / 1024 ) + 19;
            break;
        }

        default:
            cerr << "Error: should never reach here" << endl;
            exit ( -1 );
    }

    Logger_if( LOG_FOR_DEBUGGING )
    {
        paramValues.print( Logger::out(), true, AUTOMATED );
        Logger::out() << " => " << result << endl;
    }
    return result;
}


void launchBeetlBwt( BwtParamsAndEstimate &config )
{
    const string &inputFilename = config.first["input filename"];
    const string &outputFilename = config.first["output filename"];

    Logger::out() << "\nLaunching the following configuration of Beetl-bwt:" << endl;
    config.first.print( Logger::out(), false );
    Logger::out() << endl;

    Logger::out() << "Estimated RAM consumption: " << config.second.ramRssMBytes << " MB" << endl;
    Logger::out() << endl;

    switch ( config.first[PARAMETER_ALGORITHM] )
    {
        case ALGORITHM_BCR:
        {
            int bcrMode = 0 ; // 0=build BWT
            CompressionFormatType outputCompression = compressionIncrementalRunLength; // todo
            BCRexternalBWT bwt( ( char * )inputFilename.c_str(), ( char * )outputFilename.c_str(), bcrMode, outputCompression, &config.first );
        }
        break;

        case ALGORITHM_EXT:
        {
            bool isAsciiOutput = ( config.first["output format"] == OUTPUT_FORMAT_ASCII );
            bool isHuffmanOutput = ( config.first["output format"] == OUTPUT_FORMAT_HUFFMAN );
            bool isRunlengthOutput = ( config.first["output format"] == OUTPUT_FORMAT_RLE );
            bool isImplicitSort = ( config.first["SAP ordering"] == 1 );
            bool isSeqInput = ( config.first["input format"] == INPUT_FORMAT_CYC );

            assert ( isRunlengthOutput || isAsciiOutput || isHuffmanOutput );

            Algorithm *pBCRext = new BCRext( isHuffmanOutput,
                                             isRunlengthOutput,
                                             isAsciiOutput,
                                             isImplicitSort,
                                             isSeqInput,
                                             inputFilename,
                                             outputFilename );

            // run main method
            pBCRext->run();

            // clean up
            delete pBCRext;
        }
        break;
    }
}


void printUsage()
{
    params.printUsage();

    cout << "Notes:" << endl;
    cout << "    BCR only : The following options force algorithm=bcr: --reverse, --pause-between-cycle, --qualities=permute, --add-rev-comp" << endl;
    cout << "    RLE      : run-length-encoded format" << endl;
//    cout << "    multiRLE : run-length-encoded using an incremental strategy with multiple files" << endl;
    cout << "    SAP      : implicit permutation to obtain more compressible BWT" << endl;
    cout << "    LCP      : length of Longest Common Prefix shared between a BWT letter and the next one. Stored using 4 bytes per BWT letter in files with -Lxx suffix." << endl;
    cout << "               (Note: forces algorithm=bcr, non-parallel and intermediate-format=ascii)" << endl;
    cout << "               (++++ Sorry, for computing the LCP array, you must set BUILD_LCP to 1 in src/shared/Tools.hh and compile again! ++++)" << endl;   
    cout << "    PBE      : prediction-based encoding" << endl;
#ifndef _OPENMP
    cout << endl;
    cout << "Warning:" << endl;
    cout << "    This version of BEETL hasn't been compiled with OpenMP parallelisation. Performance may be degraded." << endl;
#endif //ifndef _OPENMP
    cout << endl;
}


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
    for ( int i = 0; i < argc; ++i )
    {
        cout << " " << argv[i];
    }
    cout << "\n" << endl;


    if ( !params.parseArgv( argc, argv ) || params["help"] == 1 || !params.chechRequiredParameters() )
    {
        printUsage();
        exit( params["help"] == 0 );
    }

    // Auto-detection of missing arguments
    if ( !params["input format"].isSet() )
    {
        const string &filename = params["input filename"];
        string fileFormat = detectFileFormat( filename );
        if ( fileFormat.empty() )
        {
            cerr << "Error: file format not recognised for " << filename << endl;
            exit( -1 );
        }
        params["input format"] = fileFormat;
    }
    checkFileFormat( params["input filename"], params["input format"] );

    if ( params["memory limit MB"] <= 0 )
    {
        params["memory limit MB"] = detectMemoryLimitInMB();
    }

    // Special case of SAP ordering
    /*
        if ( params["SAP ordering"] == 1 )
        {
            if ( !params["algorithm"].isSet() || strcasecmp( params["algorithm"].userValue.c_str(), "ext" ) != 0 )
            {
                clog << "Warning: Forcing algorithm=ext for --sap-ordering" << endl;
                params["algorithm"] = "ext";
            }
            if ( !params["intermediate format"].isSet() || strcasecmp( params["intermediate format"].userValue.c_str(), "ascii" ) != 0 )
            {
                clog << "Warning: Forcing intermediate-format=ascii for --sap-ordering" << endl;
                params["intermediate format"] = "ascii";
            }
        }
    */

    // Special case of LCP generation
    if ( params["generate LCP"] == 1 )
    {
        if ( !params["algorithm"].isSet() || strcasecmp( params["algorithm"].userValue.c_str(), "bcr" ) != 0 )
        {
            clog << "Warning: Forcing algorithm=bcr for --generate-lcp" << endl;
            params["algorithm"] = "bcr";
        }
        if ( !params["intermediate format"].isSet() || strcasecmp( params["intermediate format"].userValue.c_str(), "ascii" ) != 0 )
        {
            clog << "Warning: Forcing intermediate-format=ascii for --generate-lcp" << endl;
            params["intermediate format"] = "ascii";
        }
    }

    // Switches only available with the BCR algorithm
    if ( params["reverse"] == 1
         || params["pause between cycles"] == 1
         || params["process qualities"] == "permute"
         || params["add reverse complement"] == 1
       )
    {
        if ( !params["algorithm"].isSet() || strcasecmp( params["algorithm"].userValue.c_str(), "bcr" ) != 0 )
        {
            clog << "Warning: Forcing algorithm=bcr for --reverse/--pause-between-cycle/--qualities=permute/--add-rev-comp" << endl;
            params["algorithm"] = "bcr";
        }
    }

    // Use default parameter values where needed
    params.commitDefaultValues();

    datasetMetadata.init( params["input filename"], params["input format"] );
    // TODO: hardwareConstraints.init( params["HardwareConstraints"] );

    // Resource estimation
    calculateResourceRequirements();

    // Launch the fastest estimated configuration of Beetl
    if ( allEstimates.empty() )
    {
        cerr << "Error: Beetl doesn't seem to be able to run under these constraints. Try to relax some of them." << endl;
        exit( 2 );
    }
    launchBeetlBwt( allEstimates[0] );

    return 0;
}
