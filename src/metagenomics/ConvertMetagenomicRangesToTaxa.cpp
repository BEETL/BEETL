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

#include "../shared/Alphabet.hh"
#include "../shared/Tools.hh"
#include "../shared/Types.hh"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

//#define USE_MMAP // de-activated by default and controlled from Makefile.am
#ifdef USE_MMAP
# include <sys/types.h>
# include <sys/stat.h>
# include <fcntl.h>
# include <sys/mman.h>
# include <unistd.h>
# include <string.h>
#endif

using namespace std;

#include "metaShared.hh"
#include "Krona.hh"
#include "OutputTsv.hh"

//#define DEBUG


TAXMAP taxInfo;
map<int, uint64_t> lcaCounts;
map< int, double > normalisationData;
map< int, double > normalisationDataTotalCount;
map< int, double > normalisationDataOccurrences;
map< int, map< int, uint64_t > > normalisationData2;


bool sortLcaCountByValue( const pair< int, uint64_t> &obj1, const pair< int, uint64_t> &obj2 )
{
    return obj1.second > obj2.second;
}

struct InputLine //metaBeetl.out.cycle23.subset2.5:BKPT+SGLT       5       T       22      0:0:1:0:0:0     0:0:1:0:0:0     36413   743151552       1       1
{
    string header;
    uint8_t pileNum;
    char pileChar;
    uint16_t cycle;
    string countA, countB;
    LetterNumber posA, posB;
    LetterNumber numA, numB;

    friend std::istream &operator>>( std::istream &is, InputLine &obj )
    {
        is
                >> obj.header
                >> obj.pileNum
                >> obj.pileChar
                >> obj.cycle
                >> obj.countA
                >> obj.countB
                >> obj.posA
                >> obj.posB
                >> obj.numA
                >> obj.numB;
        if ( obj.pileNum >= '0' && obj.pileNum <= '6' )
            obj.pileNum -= '0';
        return is;
    }

    friend std::ostream &operator<<( std::ostream &os, const InputLine &obj )
    {
        return os
               << "{ " << obj.header
               << ", " << ( obj.pileNum + '0' )
               << ", " << obj.pileChar
               << ", " << obj.cycle
               << ", " << obj.countA
               << ", " << obj.countB
               << ", " << obj.posA
               << ", " << obj.posB
               << ", " << obj.numA
               << ", " << obj.numB
               << " }";
    }
};

struct ShortInputLineWithId
{
    uint32_t id;
    //    uint8_t pileNum;
    LetterNumber posB;
    LetterNumber numB;
    LetterNumber numA;
    ShortInputLineWithId()
        : id( 0 )
        , posB( 0 )
        , numB( 0 )
        , numA( 0 )
    {}

    ShortInputLineWithId( const uint32_t id, const InputLine &inputLine, const int minKmerSize = 0, const int maxKmerSize = 0, const bool useWeights = true )
        : id( id )
        //        , pileNum( inputLine.pileNum )
        , posB( inputLine.posB )
        , numB( inputLine.numB )
        , numA( inputLine.numA )
    {
        if ( useWeights && maxKmerSize > 0 )
        {
            double weight = 1.0;
            if ( inputLine.header.find( "DIFF" ) == string::npos )
                weight = maxKmerSize - minKmerSize + 1;
            else
                weight = inputLine.cycle - minKmerSize + 1;

            weight = max( weight, 1.0 );
            numA *= weight;
        }
    }

    friend std::ostream &operator<<( std::ostream &os, const ShortInputLineWithId &obj )
    {
        return os
               << "{ " << obj.id
               //            << ", " << obj.pileNum
               << ", " << obj.posB
               << ", " << obj.numB
               << ", " << obj.numA
               << " }";
    }

    bool operator<( const ShortInputLineWithId &rhs ) const
    {
        //        return posB < rhs.posB;
        return posB < rhs.posB || ( posB == rhs.posB && numB > rhs.numB );
        //        return pileNum < rhs.pileNum || ( pileNum == rhs.pileNum && posB < rhs.posB );
    }
};

vector< vector< int> > fileNumToTaxIds;
void loadFileNumToTaxIds( const string &taxIdFilename )
{
    ifstream taxas( taxIdFilename );
    assert( taxas.good() );

    string line;
    //each line should contain the fileNumber followed by the taxIds split up with one space character
    while ( getline( taxas, line ) )
    {
        if ( line.empty() || line[0] == '#' )
            continue;

        istringstream iss( line );

        unsigned int fileNum;
        int taxId, lastTaxId = 0;
        vector<int> taxIds;

        iss >> fileNum;
        if ( fileNum != fileNumToTaxIds.size() )
            cout << "Wrong filenumber " << fileNum << " " << fileNumToTaxIds.size() << endl;

        while ( iss >> taxId )
        {
            if ( taxId )
                lastTaxId = taxId;
            else
                taxId = lastTaxId;

            taxIds.push_back( taxId );
        }
        //test if all TaxIds were found
        if ( taxIds.size() < taxLevelSize )
            cerr << "Tax Ids don't have enough taxonomic Information. Only " << taxIds.size() << " could be found (" << line << ")" << endl
                 << "Will add unknown taxa until size is right" << endl;
        else if ( taxIds.size() > taxLevelSize )
            cerr << "Tax Ids have too much taxonomic information: " << taxIds.size() << " instead of " << taxLevelSize << "(" << line << ")" << endl
                 << "Please note, that the taxonomic information about one file should be shown as: " << endl
                 << "FileNumber Superkingdom Phylum Order Family Genus Species Strain " << endl;
        taxIds.resize( taxLevelSize );

        fileNumToTaxIds.push_back( taxIds );
    }
}

void restrictCommonTaxaToIncludeFileNum( vector<int> &commonTaxa, const MetagFileNumRefType fileNum )
// 0=undefined, -1=not match
{
    vector< int > &taxIds = fileNumToTaxIds[fileNum];
    assert( taxIds.size() == taxLevelSize );

    if ( taxIds[0] == 0 ) return; // Special case when the genome taxonomy isn't known (usually due to an incorrect genbank identifier upstream)

    for ( unsigned int i = 0; i < taxLevelSize; ++i )
    {
        switch ( commonTaxa[i] )
        {
            case -1:
                break;
            case 0:
                commonTaxa[i] = taxIds[i];
                break;
            default:
                if ( commonTaxa[i] != taxIds[i] )
                    commonTaxa[i] = -1;
                break;
        }
    }
}


void populateTaxInfoCountsUsingLcaCounts()
{
    // Resetting total counts including children
    cerr << "Resetting total counts including children" << endl;
    for ( TAXMAP::iterator iter = taxInfo.begin() ; iter != taxInfo.end(); ++iter )
        for ( unsigned int s ( 0 ); s < /*wordMinSize.size()*/1; ++s )
        {
            iter->second.wordCountPerSizeOfChildren_[s] = 0;
            const auto &lcaCount = lcaCounts.find( iter->first );
            if ( lcaCount != lcaCounts.end() )
                iter->second.wordCountPerSize_[s] = lcaCounts[iter->first];
            else
                iter->second.wordCountPerSize_[s] = 0;
        }

    // Calculating total counts including children
    cerr << "Calculating total counts including children" << endl;
    for ( int level = taxLevelSize ; level >= 0 ; --level )
    {
        for ( TAXMAP::iterator iter = taxInfo.begin() ; iter != taxInfo.end(); ++iter )
        {
            if ( ( *iter ).second.taxLevel_ == level )
            {
                //                    cerr << iter->first << " " << iter->second.parentId_ << endl;
                TAXMAP::iterator parent = taxInfo.find( iter->second.parentId_ );
                assert( parent != taxInfo.end() );
                for ( unsigned int s ( 0 ); s < /*wordMinSize.size()*/1; ++s )
                    parent->second.wordCountPerSizeOfChildren_[s] += iter->second.wordCountPerSize_[s] + iter->second.wordCountPerSizeOfChildren_[s];
            }
        }
    }
}


void kronaOutput( const string &filename )
{
    clog << "Krona output" << endl;

    ofstream output( filename );
    printKronaHeader( output );
    //        printKronaDatasets( output, wordMinSize );

    TAXMAP::iterator topLevel = taxInfo.find( 1 ); // top level has taxonomy Id 1
    assert( topLevel != taxInfo.end() );
    printKronaChildren( topLevel, output, 0, taxInfo, /*wordMinSize.size()*/1 );

    printKronaFooter( output );
    output.close();
}

void tsvOutput( const string &filename )
{
    clog << "TSV output" << endl;

    ofstream output( filename );
    printTsvHeader( output );

    TAXMAP::iterator topLevel = taxInfo.find( 1 ); // top level has taxonomy Id 1
    assert( topLevel != taxInfo.end() );
    printTsvChildren( topLevel, output, 0, taxInfo, /*wordMinSize.size()*/1 );

    output.close();
}

void parseNormalisationFile( const string &filename )
{
    // File parsing
    ifstream ifs( filename );
    string line;

    while ( getline( ifs, line ) )
    {
        vector<string> lineSplit = splitString( line, " " );
        if ( lineSplit.size() < 6 || lineSplit[4].empty() ) continue;

        clog << lineSplit[2] << ", " << lineSplit[3] << ", " << lineSplit[4] << endl;
        int taxId = stoi( lineSplit[2] );
        double count = stod( lineSplit[3] );
        double totalCount = stod( lineSplit[4] );
        normalisationData[ taxId ] += totalCount ? count / totalCount : 0;
        normalisationDataTotalCount[ taxId ] += totalCount;
        normalisationDataOccurrences[ taxId ]++;

        assert( lineSplit[5] == "Counts:" );
        for ( uint i = 6; i < lineSplit.size(); ++i )
        {
            //clog << "Splitting " << lineSplit[i] << endl;
            vector<string> tokenSplit = splitString( lineSplit[i], ":" );
            if ( tokenSplit.size() == 2 )
            {
                int token_taxId = stoi( tokenSplit[0] );
                double token_count = stod( tokenSplit[1] );
                normalisationData2[ taxId ][ token_taxId ] += token_count;
            }
        }
    }
    clog << "Parsed " << normalisationData.size() << " lines" << endl;
}

void normaliseLcaCounts( const string &normalisationFilename )
{
    parseNormalisationFile( normalisationFilename );

    // Normalisation: transfer of counts from ancestor to main
    cout << "Counts: ";
    for ( auto & item : lcaCounts )
    {
        int taxId = item.first;
        clog << "Normalisation: taxId=" << taxId << endl;
        auto it = normalisationData2.find( taxId );
        if ( it != normalisationData2.end() )
        {
            uint64_t &countInData = item.second;
            uint64_t countInDb = it->second[ taxId ];
            if ( countInDb == 0 ) // TODO: use threshold instead of zero
            {
                clog << " Normalisation database entry for " << taxId << " has a value of 0, making it impossible to normalise this genome" << endl;
                continue;
            }
            double coef = ( double )countInData / ( double )countInDb;

            clog << "Normalisation: countInData=" << countInData << ", countInDb=" << countInDb << " => coef=" << coef << endl;
            if ( coef == 0 ) continue; // TODO: use threshold instead of zero

            for ( const auto & dbItem : it->second )
            {
                int ancestorTaxId = dbItem.first;
                if ( ancestorTaxId != taxId )
                {
                    auto ancestorIter = lcaCounts.find( ancestorTaxId );
                    if ( ancestorIter != lcaCounts.end() )
                    {
                        uint64_t &ancestorCountInData = ancestorIter->second;
                        uint64_t ancestorCountInDb = dbItem.second;

                        uint64_t countsToMoveFromAncestorToMain = ancestorCountInDb * coef;

                        clog << "Normalisation: taxId=" << taxId << " ancestorTaxId=" << ancestorTaxId << " ancestorCountInDb=" << ancestorCountInDb << " => countsToMoveFromAncestorToMain=" << countsToMoveFromAncestorToMain << endl;
                        clog << "   before: countInData=" << countInData << " ancestorCountInData=" << ancestorCountInData << endl;

                        countsToMoveFromAncestorToMain = min( countsToMoveFromAncestorToMain, ancestorCountInData );
                        ancestorCountInData -= countsToMoveFromAncestorToMain;
                        countInData += countsToMoveFromAncestorToMain;

                        clog << "   after:  countInData=" << countInData << " ancestorCountInData=" << ancestorCountInData << endl;
                    }
                }
            }
        }
        else
            cerr << "hmm " << item.first << endl;
    }
}


void normaliseLcaCounts_part2()
{
    // Normalisation relatively to genome size
    uint64_t countsBeforeScaling = 0;
    uint64_t countsAfterScaling = 0;

    // Normalisation of leaf entries + computation of average scaling factors for next step
    for ( auto & item : lcaCounts )
    {
        int taxId = item.first;
        clog << "Normalisation2: taxId=" << taxId << endl;
        auto it = normalisationData2.find( taxId );
        if ( it != normalisationData2.end() )
        {
            uint64_t &countInData = item.second;
            uint64_t countInDb = it->second[ taxId ];
            uint64_t totalCountInDb = normalisationDataTotalCount[ taxId ] / normalisationDataOccurrences[ taxId ];
            if ( countInDb == 0 && countInData != 0 ) // TODO: use threshold instead of zero, for countInDb
            {
                clog << " Normalisation database entry for " << taxId << " has a value of 0, making it impossible to normalise this genome" << endl;
                continue;
            }
            clog << "countInData=" << countInData << ", countInDb=" << countInDb << ", totalCountInDb=" << totalCountInDb;
            countsBeforeScaling += countInData;
            countInData *= 1000; // rescaling
            countInData /= totalCountInDb;
            countsAfterScaling += countInData;
            clog << " => " << countInData << endl;
        }
        else
            cerr << "hmm2 " << item.first << endl;
    }

    clog << "countsBeforeScaling=" << countsBeforeScaling << endl;
    clog << "countsAfterScaling=" << countsAfterScaling << endl;

    // Normalisation of non-leaf entries (i.e. ancestors) using average scaling factors
    for ( auto & item : lcaCounts )
    {
        int taxId = item.first;
        clog << "Normalisation2: taxId=" << taxId << endl;
        uint64_t &countInData = item.second;
        auto it = normalisationData2.find( taxId );
        if ( it != normalisationData2.end() )
        {
            uint64_t countInDb = it->second[ taxId ];
            if ( countInDb == 0 && countInData != 0 ) // TODO: use threshold instead of zero, for countInDb
            {
            }
            else
            {
                clog << " 2b: Not re-normalising database entry for " << taxId << ", as we already did it before" << endl;
                continue;
            }
        }
        clog << " 2b: Normalising database entry for " << taxId << " has a value of 0, making it impossible to normalise this genome" << endl;
        clog << "countInData=" << countInData;
        countInData *= countsAfterScaling;
        countInData /= countsBeforeScaling?:1;
        clog << " => " << countInData << endl;
    }
}


void pruneUnnormalisableCounts()
{
    // Remove counts for entries where the normalisation database has a zero count
    // This should only affect leaf entries (real genomes), as ancestor entries shouldn't appear at all as normalisationData entries
    for ( auto & item : lcaCounts )
    {
        int taxId = item.first;
        clog << "Normalisation pruning: taxId=" << taxId << endl;
        auto it = normalisationData2.find( taxId );
        if ( it != normalisationData2.end() )
        {
            uint64_t &countInData = item.second;
            uint64_t countInDb = it->second[ taxId ];
            if ( countInDb == 0 ) // TODO: use threshold instead of zero
            {
                clog << " Pruning entry " << taxId << ", due to the database containing a count value of 0, making it impossible to normalise this genome" << endl;
                countInData = 0;
                continue;
            }
        }
        else
            cerr << "hmm3 " << item.first << endl;
    }
}


int main( int argc, char **argv )
{
    bool useWeights = true;
    if ( argc == 9 && string(argv[8]) == "--without-weights" )
    {
        useWeights = false;
        --argc;
    }
    if ( argc != 8 )
    {
        cerr << "Usage: " << argv[0] << " fileNumToTaxTree databaseBwtPrefix names.dmp normalisationMetadata.txt minKmerSize maxKmerSize mainInputFilename/- <--without-weights>" << endl;
        exit( -1 );
    }

    int minKmerSize = atoi( argv[5] );
    int maxKmerSize = atoi( argv[6] );
    assert( minKmerSize > 0 );
    assert( maxKmerSize >= minKmerSize );

    loadFileNumToTaxIds( argv[1] );

    InputLine inputLine;
    //vector<InputLine> inputLines;
    vector<ShortInputLineWithId> shortInputLines[alphabetSize];
    uint32_t id = 0;

    shared_ptr<istream> mainInputStream = openInputFileOrDashAsCin( string( argv[7] ) );
    while ( *mainInputStream >> inputLine )
    {
        //cout << inputLine << endl;
        //inputLines.push_back( inputLine );
        shortInputLines[inputLine.pileNum].push_back( ShortInputLineWithId( id++, inputLine, minKmerSize, maxKmerSize, useWeights ) );
    }

    #pragma omp parallel for
    for ( int i = 0; i < alphabetSize; ++i )
    {
        ShortInputLineWithId lastItem;

        if ( !shortInputLines[i].empty() )
        {
            map<int, uint64_t> lcaCountsPerThread;
            sort( shortInputLines[i].begin(), shortInputLines[i].end() );

            ostringstream cFilename;
            cFilename << argv[2] << "-C0" << i;
#ifdef USE_MMAP
            assert( sizeof( off_t ) >= 8 && "64 bits architecture required to hold large files" );
            off_t mmappedCFileSize = 0;
            char *mmappedCFile = NULL;
            {
                int fd = open( cFilename.str().c_str(), O_RDONLY );
                #pragma omp critical (IO)
                if ( fd != -1 )
                    cout << "Memory-mapped file \"" << cFilename.str() << "\"" << endl;
                else
                    cout << "ERROR: Could not open File \"" << cFilename.str() << "\"" << endl;
                mmappedCFileSize = lseek( fd, 0, SEEK_END );
                mmappedCFile = ( char * )mmap( NULL, mmappedCFileSize, PROT_READ, MAP_SHARED, fd, 0 );
                if ( mmappedCFile == ( void * ) - 1 )
                {
                    perror ( "Error mmap " );
                    mmappedCFile = NULL;
                }
                close( fd );
            }
#else
            ifstream cFile( cFilename.str() );
#endif
            int lowestCommonAncestor = 0;

            for ( auto & item : shortInputLines[i] )
            {
                if ( lastItem < item )
                {
#ifdef DEBUG
                    ostringstream outStream;
                    outStream << "Pile " << i << ": " << item << " => ";
#endif

#ifdef USE_MMAP
                    char *mmappedAddr = mmappedCFile + item.posB * sizeof( MetagFileNumRefType );
#else
                    cFile.seekg( item.posB * sizeof( MetagFileNumRefType ) );
#endif
                    vector<int> commonTaxa( taxLevelSize, 0 );
                    for ( LetterNumber j = 0; j < item.numB; ++j )
                    {
                        MetagFileNumRefType fileNum;
#ifdef USE_MMAP
                        memcpy( reinterpret_cast<char *>( &fileNum ), mmappedAddr, sizeof( MetagFileNumRefType ) );
                        mmappedAddr += sizeof( MetagFileNumRefType );
#else
                        cFile.read( reinterpret_cast<char *>( &fileNum ), sizeof( MetagFileNumRefType ) );
#endif
#ifdef DEBUG
                        outStream << fileNum << " ";
#endif

                        restrictCommonTaxaToIncludeFileNum( commonTaxa, fileNum );
                    }

#ifdef DEBUG
                    cout << outStream.str() << " => ";
#endif
                    lowestCommonAncestor = 0;
                    for ( auto & item : commonTaxa )
                    {
#ifdef DEBUG
                        cout << item << " ";
#endif
                        if ( item != 0 && item != -1 )
                            lowestCommonAncestor = item;
                    }

                    lastItem = item;
                }
                else
                {
#ifdef DEBUG
                    cout << "same";
#endif
                }

#ifdef DEBUG
                cout << " => " << lowestCommonAncestor << endl;
#endif
                if ( lowestCommonAncestor == 0 )
                    lowestCommonAncestor = 1;

                lcaCountsPerThread[lowestCommonAncestor] += item.numA;
            }

#ifdef USE_MMAP
            munmap( mmappedCFile, mmappedCFileSize );
            mmappedCFile = NULL;
#endif

            // Gather all lcaCountPerThread together
            #pragma omp critical
            for ( auto & item : lcaCountsPerThread )
                lcaCounts[item.first] += item.second;
        }
    }

    string taxInfoFilename = argv[1];
    string ncbiNames = argv[3];
    string normalisationFilename = argv[4];
    taxInfo = loadTaxInformationForDatabase( taxInfoFilename, /*wordSize.size()*/1, ncbiNames );

    populateTaxInfoCountsUsingLcaCounts();

    cout << "Counts: ";
    for ( auto & item : lcaCounts )
        cout << item.first << ":" << item.second << " ";
    cout << endl;
    kronaOutput( "metaBeetl_krona.html" );
    tsvOutput( "metaBeetl.tsv" );

    clog << "Normalisation" << endl;
    normaliseLcaCounts( normalisationFilename );
    pruneUnnormalisableCounts();
    populateTaxInfoCountsUsingLcaCounts();
    /*
        cout << "Normalised counts: ";
        for ( auto & item : lcaCounts )
            cout << item.first << ":" << item.second << " ";
        cout << endl;
    */
    kronaOutput( "metaBeetl_krona_normalised.html" );
    tsvOutput( "metaBeetl_normalised.tsv" );

    normaliseLcaCounts_part2();
    populateTaxInfoCountsUsingLcaCounts();
    kronaOutput( "metaBeetl_krona_normalised2.html" );
    tsvOutput( "metaBeetl_normalised2.tsv" );

    // Final ranking
    vector< pair< int, uint64_t> > sortedLcaCounts;
    for ( auto & item : lcaCounts )
        sortedLcaCounts.push_back( item );
    sort( sortedLcaCounts.begin(), sortedLcaCounts.end(), sortLcaCountByValue );

    cout << "Ordered normalised counts: ";
    uint64_t sum = 0;
    for ( auto & item : sortedLcaCounts )
        if ( normalisationData.find( item.first ) != normalisationData.end() )
            sum += item.second;
    for ( auto & item : sortedLcaCounts )
        if ( normalisationData.find( item.first ) != normalisationData.end() )
            cout << item.first << ":" << item.second << " (" << ( 100 * item.second ) / sum << "%) ";
    cout << endl;

    return 0;
}
