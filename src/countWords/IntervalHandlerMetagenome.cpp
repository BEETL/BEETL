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

#include "IntervalHandlerMetagenome.hh"

#include "LetterCount.hh"
#include "libzoo/util/Logger.hh"

#include <fcntl.h>
#include <map>

using namespace std;


IntervalHandlerMetagenome::IntervalHandlerMetagenome( unsigned int minOcc,
        vector<string> &filenamesCSet,
        const vector<char *> &mmappedCFiles,
        vector< vector< int> > &fileNumToTaxIds,
        bool testDB,
        uint minWordLength,
        uint maxWordLength )
    : minOcc_( minOcc ),
      mmappedCFiles_( mmappedCFiles ),
      fileNumToTaxIds_( fileNumToTaxIds ), testDB_( testDB ),
      minWordLength_( minWordLength ),
      maxWordLength_( maxWordLength )
{
    cSetFileDescs_.resize( alphabetSize );
    posInFile_.resize( alphabetSize );

    for ( unsigned int i = 0; i < alphabetSize; i++ )
    {
        if ( mmappedCFiles_.size() <= i || mmappedCFiles_[i] == NULL )
        {
            int fd = open( filenamesCSet[i].c_str(), O_RDONLY );
            if ( fd == -1 )
            {
                #pragma omp critical (IO)
                {
                    cout << "ERROR: Could not open File \"" << filenamesCSet[i] << "\"" << endl;
                }
            }
            cSetFileDescs_[i] = fd;
            posInFile_[i] = 0;
        }
    }
}

IntervalHandlerMetagenome::~IntervalHandlerMetagenome()
{
    for ( int i = 0; i < alphabetSize; i++ )
    {
        if ( cSetFileDescs_[i] )
            close( cSetFileDescs_[i] );
    }
}

/*  @ int pileNum the number in which pile this suffix can be found. So it is known which pile can be used to look up the filenumbers
    @ ulong bwtPosition at this position the information in the file starts from which files the suffix came from
    @ unsigned int num Count of how many times the suffix is found in the database. This is the same number as fileNumber which have to read
    @ vector<MetagFileNumRefType> fileNumbers Vector to fill up with unique fileNumbers in which the suffix can be found.
    If the size of the vector is the same as the num the suffix occours only once in each file
    Returns nothing, but fills up the fileNumbers vector
    get the fileNumbers in a certain Range
*/
void IntervalHandlerMetagenome::getFileNumbersForRange( const unsigned int &pileNum, const LetterNumber &bwtPosition, const unsigned int &num, vector<MetagFileNumRefType> &fileNumbers )
{
    assert( fileNumbers.empty() );
    //fast fix because merged genomes should have no pile 6
    if ( pileNum == 6 )
        return;

    off_t wantedPos = ( bwtPosition & matchMask ) * sizeof( MetagFileNumRefType );
    MetagFileNumRefType *fileNum = NULL;

    const bool useMmap = ( mmappedCFiles_.size() > pileNum && mmappedCFiles_[pileNum] != NULL );
    if ( useMmap )
    {
        fileNum = ( MetagFileNumRefType * ) & ( mmappedCFiles_[pileNum][wantedPos] );
    }
    else
    {
        //go to the position in the File where the fileNumbers for each Bwt positions are indicated
        off_t currentPos = posInFile_[pileNum]; //lseek( cSetFileDescs_[pileNum], 0, SEEK_CUR );
        // We only use fseek if the new position is far from the current one. Otherwise we read data into a discarded buffer
        const off_t maxDistance = 32 * 1024; // 32KB
        if ( wantedPos >= currentPos && wantedPos < ( currentPos + maxDistance ) )
        {
            const off_t offset = wantedPos - currentPos;
            void *jumpOverData = malloc ( offset );
            assert( read( cSetFileDescs_[pileNum], jumpOverData, offset ) == offset );
            free( jumpOverData );
        }
        else
        {
            lseek( cSetFileDescs_[pileNum] , wantedPos, SEEK_SET );
        }
        //reserve the space for the fileNumbers. there are as many fileNumbers as the suffix count
        ssize_t fileNumbersBufSize = num * sizeof( MetagFileNumRefType );
        fileNum = ( MetagFileNumRefType * ) malloc ( fileNumbersBufSize );
        //read out the fileNumbers which are corresponding to this BWT-postion
        if ( read( cSetFileDescs_[pileNum], fileNum, fileNumbersBufSize ) != fileNumbersBufSize )
        {
            cerr << "Error reading " << fileNumbersBufSize << " bytes from cSetFileDescs pile " << pileNum << " at position " << wantedPos << endl;
            assert( false );
        }
        posInFile_[pileNum] = wantedPos + fileNumbersBufSize;
    }

    //to test if the suffix is a singleton in each of the different files the found fileNumbers
    //will be added to an vector, only if the fileNumber is not already there
    fileNumbers.push_back( fileNum[0] );
    for ( unsigned int i( 1 ); i < num; ++i )
    {
        if ( fileNumToTaxIds_[fileNum[i]][2] >= 0 )
        {
            bool alreadyAdded( false );
            for ( unsigned int j( 0 ); j < fileNumbers.size() ; ++j )
            {
                //if the fileNumber is already in the vector this means this suffix was found twice in one file and is no singleton
                //late on the size of the fileNumbers and the range number can be tested
                if ( fileNumbers[j] == fileNum[i] )
                {
                    alreadyAdded = true;
                    break;
                }
            }
            //add only fileNumbers which are not already in there
            if ( !alreadyAdded )
                fileNumbers.push_back( fileNum[i] );
        }
    }

    if ( !useMmap )
    {
        //free memory
        free( fileNum );
    }
    ///  Logger::out() << "got fileNumbers" <<endl;
}

/*
   @vector<unsigned int> sharedTaxIds vector of length taxLevelSize to fill up with the taxIds which are the same between file Numbers
   will contain one taxId for each level of the taxonomy, or nothing on the taxonomic levels if the files have no common taxa
   @vector<MetagFileNumRefType> fileNumbers file numbers where a suffix was found
   Returns vector<bool> intervalInSameTaxa vector of length taxLevelSize, contains true or false for each taxonomic Level,
   depending if the file numbers have a shared taxonomy

   TODO: could be made easier with just giving back the sharedTaxIds and a zero in the vector for each return value taxlevel where there is no shared taxa
*/
vector<bool> IntervalHandlerMetagenome::intervalInSameTaxa( vector<unsigned int> &sharedTaxIds, vector<MetagFileNumRefType> &fileNumbers )
{
    //first get the matching fileNumbers out of the file with the file numbers corresponding to bwt positions of the merging.
    vector<bool> taxSame;
    taxSame.resize( taxLevelSize );

    // if there was only one file it has to have a common taxa in all taxLevels
    if ( fileNumbers.size() == 1 )
    {
        for ( unsigned int i( taxLevelSize - 1 ); i > 0; i-- )
        {
            int taxIdForFile0 = fileNumToTaxIds_[fileNumbers[0]][i];
            if ( taxIdForFile0 != -1 )
            {
                taxSame[i] = true;
                //      Logger::out() << "FileNumbers size 1 " << fileNumbers[0] << ' ' << fileNumToTaxIds_[fileNumbers[0]][taxLevelSize-1] <<endl;
                sharedTaxIds[i] = taxIdForFile0;
            }
        }
        return taxSame;
    }

    //look if the files have at each point of the taxonomic tree different tax level or if they are the same
    //  Logger::out() << "intervalInSame " << fileNumbers.size() <<endl;
    for ( unsigned int i = 0 ;  i < taxLevelSize ; ++i )
    {
        // to remove outlier, first count all taxas
        map<int, int> taxaCount;
        //for each taxa for each level count how many were found
        for ( unsigned int j( 0 ); j < fileNumbers.size(); ++j )
        {
            MetagFileNumRefType &fileNumber = fileNumbers[j];
            int taxId = fileNumToTaxIds_[fileNumber][i];
            if ( taxId != -1 )
                taxaCount[taxId] += 1 ; // automatically creates element if it doesn't exist
        }
#ifdef MET_DEBUG
        #pragma omp critical (IO)
        {
            Logger::out() << i << endl;
            for ( map<int, int>::iterator it = taxaCount.begin(); it != taxaCount.end(); ++it )
                Logger::out() << ( *it ).first << ' ' << ( *it ).second << endl;
            Logger::out() << i << " taxa " << taxaCount.size() << endl;
        }
#endif
        //if there was only one possible taxa no outliner is possible
        if ( taxaCount.size() == 1 )
        {
            int taxId = taxaCount.begin()->first;
            assert( taxId >= 0 );
            taxSame[i] = true;
            sharedTaxIds[i] = taxId;
            continue;
        }
        else
            // Current algorithm: Don't select taxa at levels where multiple are present
            // Previously: Select any taxa that is present more than 50% of the time
            // Previously: first remove all taxa which only make up less than 10% of the possible taxas
            //   this could be made much better if it is taking in account of how many different taxa on one level are even possible
            //   this also heightens the amount of taxa found in extreme cases
        {
            /*
                        float countThreshold = 0.50 * fileNumbers.size();
                        vector<int> biggestTaxa;
                        for ( map<int, int>::iterator it = taxaCount.begin() ; it != taxaCount.end(); it++ )
                        {
                            if ( it->second > countThreshold && it->first != -1 )
                            {
                                biggestTaxa.push_back( ( *it ).first );
                            }
                        }
                        //    Logger::out() << i << " remaining taxa " << taxaCount.size() <<endl;
                        //if this leaves only one remaining taxa use that
                        if ( biggestTaxa.size() == 1 )
                        {
                            taxSame[i] = true;
                            sharedTaxIds[i] = biggestTaxa[0];
                        }
                        // else there is no dominant taxa in the lot
                        // so it is assumed that the files have no common taxonomy
                        else
            */
            {
                taxSame[i] = false;
                break;
            }
        }
    }
#ifdef MET_DEBUG
    for ( unsigned int i( 0 ); i < fileNumbers.size(); i++ )
        cerr << fileNumbers[i] << ' ' ;
    cerr << endl;
    for ( unsigned int i( 0 ); i < taxSame.size() ; i++ )
        cerr << sharedTaxIds[i] << ' ' ;
    cerr << endl;
    for ( unsigned int i( 0 ); i < taxSame.size(); i++ )
        cerr << taxSame[i] << ' ' ;
    cerr << endl;
#endif
    return taxSame;
}


void IntervalHandlerMetagenome::foundInBoth
( const int pileNum,
  const LetterCount &countsThisRangeA,
  const LetterCount &countsThisRangeB,
  const Range &thisRangeA,
  const Range &thisRangeB,
  AlphabetFlag &propagateIntervalA,
  AlphabetFlag &propagateIntervalB,
  bool &isBreakpointDetected,
  const int cycle
)
{
    //propagate A only if the Range also exists in B
    //print only if the word can't be propagated anymore,
    //or the start of the read is reached
    //only print if reads and referenceSequences can be only differently propagated
    //for example if the propagation A and C in the database is possible and in the reads the propagation A G is possible
    //print the taxonomic information because the different propagation between C and G could be a breakpoint between taxas
    bool differentProp( false );
    bool belowMinDepthBecauseOfEndOfReads( false );
    //    bool endOfRead( false );
    bool atLeastOneProp = false;
    bool singletonInRef = ( thisRangeB.num_ == 1 );

    for ( int l ( 1 ); l < alphabetSize; l++ )
    {
        if ( countsThisRangeB.count_[l] > 0 )
        {
            const bool propA = countsThisRangeA.count_[l] >= minOcc_;
            propagateIntervalA[l] = propA;
            propagateIntervalB[l] = propA || testDB_;
            atLeastOneProp |= propA;
            // If k-mer is not propagated because it slowly falls below threshold due to end-of-reads, we set this flag to output it
            if ( !propA && countsThisRangeA.count_[l] + countsThisRangeA.count_[0] >= minOcc_ )
                belowMinDepthBecauseOfEndOfReads = true;
        }
    }
    for ( int l ( 1 ); l < alphabetSize; l++ )
    {
        if ( countsThisRangeB.count_[l] == 0 )
        {
            //if the range in B can't be propagated anymore test if the propagation in A would have been possible
            //if yes print the taxonomic information but don't propagate anymore
            if ( countsThisRangeA.count_[l] >= minOcc_ && !atLeastOneProp )
                differentProp = true;
            propagateIntervalB[l] = false;
            propagateIntervalA[l] = false;
        }
    }
    // If some propagation is present, don't bother taking into account the potential that dollar signs cover a good call
    // Also, if some propagation is present, don't actually try to print anything
    if ( atLeastOneProp )
    {
        belowMinDepthBecauseOfEndOfReads = false;
        differentProp = false;
    }

    //if word is a singleton in the reference files don't propagate anymore
    if ( singletonInRef && cycle >= ( int )minWordLength_ )
    {
        for ( int l ( 1 ); l < alphabetSize; l++ )
        {
            propagateIntervalB[l] = false;
            propagateIntervalA[l] = false;
        }
    }

    //if at least one read reached the end print the taxonomic information (end of read detected by count['$']>0)
    //    if ( countsThisRangeA.count_[0] > 0 )
    //        endOfRead = true;
    //flag indicating if information about the database should be given
    if ( testDB_ )
    {
        //if the database should be tested propagate the intervals in B as long as they can be propagated
        //this should have been done anyway with the normal propagation, this is just to make sure
        for ( int l( 1 ); l < alphabetSize; l++ )
        {
            propagateIntervalB[l] = countsThisRangeB.count_[l] > 0;
        }
    }

    //don't bother about the Ns
    propagateIntervalB[whichPile[( int )dontKnowChar]] = false;
    propagateIntervalA[whichPile[( int )dontKnowChar]] = false;
    //testing if the word is in the reads and in certain genomes
    //only look the fileNumbers up if the word can't be propagated anymore.
    //so for each BWT Position only the longest word will be taken into account.
    //unique file numbers from which the suffixes came from
    vector<MetagFileNumRefType> fileNumbers;
    //shared taxa between the files, for each taxonomic level there should be one shared taxonomic id indicating if the
    vector<unsigned int> sharedTaxa;
    sharedTaxa.resize( taxLevelSize );
    //first look if there is interesting information about the database to be found
    if ( testDB_ )
    {
        //only print information about a certain singleton lengths
        //this is mostly to sorten the output of the database information as well as reduce running time
        //the numbers are completely randomly chosen
        if ( cycle == 35 || ( cycle % 50 ) == 0 )
        {
            // get the unique file numbers where this Range can be found
            getFileNumbersForRange( pileNum, thisRangeB.pos_, thisRangeB.num_, fileNumbers );
            // test if word is singleton in the files
            #pragma omp critical (IO)
            Logger::out() << "file numbers size " << fileNumbers.size() << " range " << thisRangeB.num_ << endl;
            if ( fileNumbers.size() <= thisRangeB.num_ )
            {
                //if the word is a singleton in the files, look if there is a shared taxonomic Id on the different taxonomic levels
                vector<bool> sameTaxa = intervalInSameTaxa( sharedTaxa, fileNumbers );
                //go through the taxonomic levels starting with the lowest (either strain or species)
                for ( unsigned int i( taxLevelSize - 1 ); i > 1; i-- )
                {
                    //look for each taxonomic level if there is a shared taxa and print the minimal information about this leven
                    //break if a shared information was found.
                    if ( sameTaxa[i] && sharedTaxa[i] != 0 )
                    {
                        //taxaCountPerWordLength[cycle][sharedTaxa[i]] += 1;
                        #pragma omp critical (IO)
                        Logger::out() << "BWord" << cycle << "|" << sharedTaxa[i] << "|";
                        break;
                    } //~same taxa
                }// ~ for loop for the tax levels
            }//~test of singletons
        }// ~if wordLengthtest
    }//~testDB


    //if the database should not be tested on its own the normal comparison between the two datasets takes place
    //the suffix of the range can only be once in each file if the amount of suffixes is smaller than all possible files. so start testing if it is smaller
    else if ( thisRangeB.num_ < fileNumToTaxIds_.size()
              //test file numbers and print information only if the suffix is longer than a certain minimal length chosen from the user
              //this reduces the amount of output and speeds up the count word algorithm because the filenumbers don't need to be looked at
              &&  cycle >= ( int )minWordLength_
              //test file numbers and print taxonomic information only if a possible tax break point or the end of the read is reached
            )
    {
        sharedTaxa.resize( taxLevelSize );

        const bool printBKPT = false;
        //#define DEBUG_DONT_CALCULATE_TAXA

        if ( printBKPT && ( differentProp || belowMinDepthBecauseOfEndOfReads || singletonInRef ) )
        {
            isBreakpointDetected = true;
            #pragma omp critical (IO)
            {
                Logger::out() << "BKPT ";
                if ( thisRangeB.word_.empty() )
                    Logger::out() << alphabet[pileNum] << string( cycle - 1, 'x' ); // No propagated sequence => Print what we know of the sequence
                else
                    Logger::out() << thisRangeB.word_;
                Logger::out()
                        << ' ' << countsThisRangeA.count_[0]
                        << ':' << countsThisRangeA.count_[1]
                        << ':' << countsThisRangeA.count_[2]
                        << ':' << countsThisRangeA.count_[3]
                        << ':' << countsThisRangeA.count_[4]
                        << ':' << countsThisRangeA.count_[5]
                        << ' ' << countsThisRangeB.count_[0]
                        << ':' << countsThisRangeB.count_[1]
                        << ':' << countsThisRangeB.count_[2]
                        << ':' << countsThisRangeB.count_[3]
                        << ':' << countsThisRangeB.count_[4]
                        << ':' << countsThisRangeB.count_[5]
                        << ' ' << ( thisRangeA.pos_ & matchMask )
                        << ' ' << ( thisRangeB.pos_ & matchMask )
                        << '\n';
            }
        }

#ifndef DEBUG_DONT_CALCULATE_TAXA
        bool maxLengthReached = ( cycle >= ( int )maxWordLength_ );
        if ( differentProp || belowMinDepthBecauseOfEndOfReads || maxLengthReached || singletonInRef )
        {
#define SIMPLE_OUTPUT_WITHOUT_FILE_NUMBERS
#ifdef SIMPLE_OUTPUT_WITHOUT_FILE_NUMBERS
            outFile_ << "BKPT";
            if ( differentProp )
                outFile_ << "+DIFF";
            if ( singletonInRef )
                outFile_ << "+SGLT"; // singleton in ref
            if ( thisRangeA.num_ == countsThisRangeA.count_[0] )
                outFile_ << "+EOR"; // end of read
            if ( maxLengthReached )
                outFile_ << "+MAX";

            outFile_
                    << ' ' << pileNum
                    << ' ' << alphabet[pileNum]
                    << ' ' << ( cycle - 1 )
                    << ' ' << countsThisRangeA.count_[0]
                    << ':' << countsThisRangeA.count_[1]
                    << ':' << countsThisRangeA.count_[2]
                    << ':' << countsThisRangeA.count_[3]
                    << ':' << countsThisRangeA.count_[4]
                    << ':' << countsThisRangeA.count_[5]
                    << ' ' << countsThisRangeB.count_[0]
                    << ':' << countsThisRangeB.count_[1]
                    << ':' << countsThisRangeB.count_[2]
                    << ':' << countsThisRangeB.count_[3]
                    << ':' << countsThisRangeB.count_[4]
                    << ':' << countsThisRangeB.count_[5]
                    << ' ' << ( thisRangeA.pos_ & matchMask )
                    << ' ' << ( thisRangeB.pos_ & matchMask )
                    << ' ' << thisRangeA.num_
                    << ' ' << thisRangeB.num_
                    << '\n';

#else // SIMPLE_OUTPUT_WITHOUT_FILE_NUMBERS

            //get the unique file number for the files where the suffix can be found in
            getFileNumbersForRange( pileNum, thisRangeB.pos_, thisRangeB.num_, fileNumbers );

            const bool printBKPTdetails = false;
            if ( printBKPTdetails )
                #pragma omp critical (IO)
            {
                Logger::out() << "BKPTfiles ";
                for ( unsigned int f( 0 ) ; f < fileNumbers.size(); f++ )
                    Logger::out() << fileNumbers[f] << ':';
                Logger::out() << '\n';
            }

            //only get taxonomic information if the word is a singleton in each file
            //unique file numbers must be the same as the numbers of suffixes, so that each file has exactly one suffix
            assert ( thisRangeB.num_ >= fileNumbers.size() );
            //#define IGNORE_KMERS_FROM_REPEATS
#ifdef IGNORE_KMERS_FROM_REPEATS
            if ( thisRangeB.num_ == fileNumbers.size() )
#endif // IGNORE_KMERS_FROM_REPEATS
            {
                vector<bool> sameTaxa = intervalInSameTaxa( sharedTaxa, fileNumbers );

                //print the information about a shared superkingdom only if a wordlength higher than 50 is reached
                //printing this information about shorter words inflates the output incredibly without giving much of useful information
                int smallestTaxLevel = ( cycle > 50 ) ? 0 : 1;

                //if there are shared taxonomic ids between the files where this suffix can be found
                //print the information about the deepest tax level
                for ( int i( taxLevelSize - 1 ) ; i >= smallestTaxLevel; --i )
                {
                    if ( sameTaxa[i] && sharedTaxa[i] != 0 )
                    {
                        //                        if ( !( differentProp || belowMinDepthBecauseOfEndOfReads ) )
                        //                            Logger::out() << "ending"; // == This MTAXA is due to maxLengthReached
                        /*
                                                if ( singletonInRef )
                                                    outFile_
                                                            << "MTAXA" << cycle
                                                            << ' ' << i
                                                            << ' ' << sharedTaxa[i]
                                                            << ' ' << thisRangeA.num_
                                                            << '\n';
                                                else
                        */
                        {
                            outFile_
                                    << "MTAXA " << i
                                    << ' ' << sharedTaxa[i]
                                    << ' ';
                            if ( thisRangeB.word_.empty() )
                                outFile_ << alphabet[pileNum] << string( cycle - 1, 'x' ); // No propagated sequence => Print what we know of the sequence, as downstream tool 'metabeetl-parseMetagenomeOutput' is just using first char and length of this string
                            else
                                outFile_ << thisRangeB.word_;
                            outFile_
                                    << ' ' << ( thisRangeB.pos_ & matchMask )
                                    << ' ' << countsThisRangeA.count_[0]
                                    << ':' << countsThisRangeA.count_[1]
                                    << ':' << countsThisRangeA.count_[2]
                                    << ':' << countsThisRangeA.count_[3]
                                    << ':' << countsThisRangeA.count_[4]
                                    << ':' << countsThisRangeA.count_[5]
                                    << ' ' << countsThisRangeB.count_[0]
                                    << ':' << countsThisRangeB.count_[1]
                                    << ':' << countsThisRangeB.count_[2]
                                    << ':' << countsThisRangeB.count_[3]
                                    << ':' << countsThisRangeB.count_[4]
                                    << ':' << countsThisRangeB.count_[5]
                                    << ' ';
                            for ( unsigned int f( 0 ) ; f < fileNumbers.size(); f++ )
                                outFile_ << fileNumbers[f] << ':';
                            outFile_ << '\n';
                        }

                        /*
                                                #pragma omp critical (IO)
                                                {
                                                    Logger::out() << oss.str();
                                                }
                        */
                        break;
                    }//~if shared taxa
                }//~for loop for the taxonomic levels
            }
#endif // SIMPLE_OUTPUT_WITHOUT_FILE_NUMBERS
        }
#endif //ifndef DEBUG_DONT_CALCULATE_TAXA
    }//~ if matches printing criteria
}

void IntervalHandlerMetagenome::foundInBOnly (
    const int pileNum,
    const LetterCount &countsSoFarB,
    const LetterCount &countsThisRangeB,
    const char *bwtSubstring,
    Range &thisRangeB,
    AlphabetFlag &propagateIntervalB,
    const int cycle
)
{
    //this hander is only interesting if you want to get information about the alignability of the database
    if ( testDB_ )
    {
        //propagate as long as possible
        for ( int l( 1 ); l < alphabetSize; l++ )
        {
            propagateIntervalB[l] = countsThisRangeB.count_[l] > 0;
            if ( pileNum == 0 )
                cerr << countsThisRangeB.count_[l] << ' ' ;
        }
        if ( pileNum == 0 )
            cerr << endl;
        //don't bother about the Ns they are creating to many false positives
        propagateIntervalB[whichPile[( int )dontKnowChar]] = false;
        //only print information about a certain word length to keep the output smaller
        if ( cycle == 35 || ( cycle % 50 ) == 0 )
        {
            //only try this if there even is a range like that
            //not sure why this happended a few times
            if ( thisRangeB.num_ > 0 )
            {
                vector<MetagFileNumRefType> fileNumbers;
                vector<unsigned int> sharedTaxa;
                sharedTaxa.resize( taxLevelSize );
                //get the file numbers in which the suffix can be found
                getFileNumbersForRange( pileNum, thisRangeB.pos_, thisRangeB.num_, fileNumbers );
                //test if the suffix is a singleton in the files
                if ( ( unsigned int )thisRangeB.num_ == fileNumbers.size() )
                {

                    vector<bool> sameTaxa;
                    //get the shared taxonomic information
                    sameTaxa = intervalInSameTaxa( sharedTaxa, fileNumbers );
                    //if there is some shared taxonomic information, print the information about the lowest taxonomic level
                    //in which taxonomic information can be found
                    for ( unsigned int i( taxLevelSize - 1 ); i > 1; i-- )
                    {
                        if ( sameTaxa[i] && sharedTaxa[i] != 0 )
                        {
                            // print the smallest possible information to keep the output small
                            #pragma omp critical (IO)
                            Logger::out() << "BWord" << cycle << "|" << sharedTaxa[i] << "|";
                            break;
                        }
                    }
#ifdef OLD
                    //printing all this created an output bigger than 500GB
                    for ( unsigned int i( 5 ) ; i != 0; i-- )
                    {
                        if ( sameTaxa[i] )
                            #pragma omp critical (IO)
                        {
                            Logger::out() << "BTAXA " << taxLevelNames[i] <<  ' ' << sharedTaxa[i] << ' ' << thisRangeB.word_ << ' ' ;
                            Logger::out() << ( thisRangeB.pos_ & matchMask ) << ' ' << thisRangeB.num_ << '\n';
                        }
                    }
                    if ( sameTaxa[taxLevelSize - 1] )
                        #pragma omp critical (IO)
                    {
                        speciesFound = true;

                        Logger::out()
                                << "BSPECIES"
                                << ' ' << thisRangeB.word_
                                << ' ' << ( thisRangeB.pos_ & matchMask )
                                << ' ' << sharedTaxa[taxLevelSize - 1]
                                << ' ' << countsThisRangeB.count_[0]
                                << ':' << countsThisRangeB.count_[1]
                                << ':' << countsThisRangeB.count_[2]
                                << ':' << countsThisRangeB.count_[3]
                                << ':' << countsThisRangeB.count_[4]
                                << ':' << countsThisRangeB.count_[5]
                                << ' ';
                        for ( unsigned int j( 0 ) ; j < fileNumbers.size() ; j++ )
                            Logger::out() << fileNumbers[j] << ':' ;
                        Logger::out() << '\n';
                    }
                    if ( fileNumbers.size() == 1 )
                        #pragma omp critical (IO)
                    {
                        Logger::out()
                        "BSINGLE"
                                << ' ' << thisRangeB.word_
                                << ' ' << ( thisRangeB.pos_ & matchMask )
                                << ' ' << countsThisRangeB.count_[0]
                                << ':' << countsThisRangeB.count_[1]
                                << ':' << countsThisRangeB.count_[2]
                                << ':' << countsThisRangeB.count_[3]
                                << ':' << countsThisRangeB.count_[4]
                                << ':' << countsThisRangeB.count_[5]
                                << ' ';
                        for ( unsigned int j( 0 ) ; j < fileNumbers.size() ; j++ )
                            Logger::out() << fileNumbers[j] << ':';
                        Logger::out() << '\n';
                    }
#endif
                }//~test on singletons
            }//~test if there is a number of ranges
        }//~if certain word length
    }//~testdb
    //if the database should not be tested don't propagate ranges only found there
    else
    {
        for ( unsigned int l ( 1 ); l < alphabetSize; l++ )
        {
            propagateIntervalB[l] = false;
        }
    }
}

void IntervalHandlerMetagenome::foundInAOnly
( const int pileNum,
  const LetterCount &countsSoFarA,
  const LetterCount &countsThisRangeA,
  const char *bwtSubstring,
  Range &thisRangeA,
  AlphabetFlag &propagateIntervalA,
  const int cycle
)
{
    //this should be interesting if one wants to find high accounted reads in A which have no singletons in the database
    // special interesting for the comparison between metatranscriptomes and metagenomes
    //maybe having a second minOcc, or using the minOcc only here and not in the handler FoundInBoth, because there one wants every occurence
    //  cerr <<"found in A only is that really needed?" <<endl;
    //Is that really interesting?
    //should probably always propagate A if higher than minOcc

    for ( int l( 0 ); l < alphabetSize; l++ )
        propagateIntervalA[l] = false;
}
