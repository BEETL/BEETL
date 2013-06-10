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

#include "IntervalHandlerMetagenome.hh"

#include "LetterCount.hh"

#include <map>

using namespace std;


struct TaxInformation
{
    unsigned short taxLevel_;
    vector<int> files_;
    vector < vector <unsigned long> > bwtPositions_;
    vector < vector <double> > wordCounts_;
    vector < vector <unsigned short> >wordLengths_;
    vector <vector <bool> > genomeSingletons_;
    unsigned int parentId_;
    unsigned long seqLengthSum_;
};

map< unsigned int, TaxInformation> TAXMAP;

/*  @ int pileNum the number in which pile this suffix can be found. So it is known which pile can be used to look up the filenumbers
    @ ulong bwtPosition at this position the information in the file starts from which files the suffix came from
    @ unsigned int num Count of how many times the suffix is found in the database. This is the same number as fileNumber which have to read
    @ vector<MetagFileNumRefType> fileNumbers Vector to fill up with unique fileNumbers in which the suffix can be found.
    If the size of the vector is the same as the num the suffix occours only once in each file
    Returns nothing, but fills up teh fileNumbers vector
    get the fileNumbers in a certain Range
*/
void IntervalHandlerMetagenome::getFileNumbersForRange( const int &pileNum, const ulong &bwtPosition, const unsigned int &num, vector<MetagFileNumRefType> &fileNumbers )
{
    //fast fix because merged genomes should have no pile 6
    if ( pileNum == 6 )
        return;
    //go to the position in the File where the fileNumbers for each Bwt positions are indicated
    fseek( mergeCSet_[pileNum] , ( ( bwtPosition & matchMask ) * sizeof( MetagFileNumRefType ) ), SEEK_SET );
    //reserve the space for the fileNumbers. there are as many fileNumbers as the suffix count
    MetagFileNumRefType *fileNum = ( MetagFileNumRefType * ) malloc ( num * ( sizeof( MetagFileNumRefType ) ) );
    //read out the fileNumbers which are corresponding to this BWT-postion
    fread( fileNum, sizeof( unsigned short ), num, mergeCSet_[pileNum] );

    //to test if the suffix is a singleton in each of the different files the found fileNumbers
    //will be added to an vector, only if the fileNumber is not already there
    fileNumbers.push_back( fileNum[0] );
    for ( unsigned int i( 1 ); i < num; ++i )
    {
        bool alreadyAdded( false );
        for ( unsigned int j( 0 ); j < fileNumbers.size() ; ++j )
        {
            //if the fileNumber is already in the vector this means this suffix was found twice in one file and is no singleton
            //late on the size of the fileNumbers and the range number can be tested
            if ( fileNumbers[j] == fileNum[i] || fileNumToTaxIds_[fileNum[i]][2] == -1 )
            {
                alreadyAdded = true;
                break;
            }
        }
        //add only fileNumbers which are not already in there
        if ( !alreadyAdded )
            fileNumbers.push_back( fileNum[i] );
    }
    //free memory
    delete fileNum;
    ///  cout << "got fileNumbers" <<endl;
}

/*
   @vector<unsigned int> sharedTaxIds vector of length taxLevelSize to fill up with the taxIds which are the same between file Numbers
   will contain one taxId for each level of the taxonomy, or nothing on the taxonomic levels if the files have no common taxa
   @vector<unsigned short >fileNumbers file numbers where a suffix was found
   Returns vector<bool> intervalInSameTaxa vector of length taxLevelSize, contains true or false for each taxonomic Level,
   depending if the file numbers have a shared taxonomy

   TODO: could be made easier with just giving back the sharedTaxIds and a zero in the vector for each return value taxlevel where there is no shared taxa
*/
vector<bool> IntervalHandlerMetagenome::intervalInSameTaxa( vector<unsigned int> &sharedTaxIds, vector<MetagFileNumRefType> &fileNumbers )
{
    //first get the matching fileNumbers out of the file with the file numbers corresponding to bwt positions of the merging.
    vector<bool> taxSame;
    taxSame.resize( taxLevelSize );
    //look if the files have at each point of the taxonomic tree different tax level or if they are the same
    bool sameTaxa( false );
    //  cout << "intervalInSame " << fileNumbers.size() <<endl;
    for ( int i = taxLevelSize - 1  ;  i >= 0 ; i-- )
    {
        // to remove outlier, first count all taxas
        map<int, int> taxaCount;
        //for each taxa for each level count how many where found
        for ( unsigned int j( 0 ); j < fileNumbers.size(); ++j )
        {
            if ( fileNumToTaxIds_[j][i] == -1 )
            {
                if ( taxaCount.find( fileNumToTaxIds_[j][i] ) == taxaCount.end() )
                {
                    taxaCount[fileNumToTaxIds_[fileNumbers[j]][i]] = 1;
                }
                else
                {
                    taxaCount[fileNumToTaxIds_[fileNumbers[j]][i]] += 1 ;
                }
            }
        }
#ifdef MET_DEBUG
        cout << i << endl;
        for ( map<int, int>::iterator it = taxaCount.begin(); it != taxaCount.end(); it++ )
            cout << ( *it ).first << " " << ( *it ).second << endl;
        cout << i << " taxa " << taxaCount.size() << endl;
#endif
        taxSame[i] = sameTaxa;
        //if there was only one possible taxa no outliner is possible
        bool testOutliner = true;
        if ( taxaCount.size() == 1 && fileNumToTaxIds_[fileNumbers[0]][i] != -1 )
        {
            taxSame[i] = true;
            sharedTaxIds[i] = fileNumToTaxIds_[fileNumbers[0]][i];
            break;
            //      testOutliner = false;
        }
        //    else if(taxaCount.size() < 3) {
        //first remove all taxa which only make up less than 10% of the possible taxas
        // this could be made much better if it is taking in account of how many different taxa on one level are even possible
        // this also heightens the amount of taxa found in extreem
        if ( testOutliner )
        {
            vector<int> biggestTaxa;
            for ( map<int, int>::iterator it = taxaCount.begin() ; it != taxaCount.end(); it++ )
            {
                //if there are two taxa, where one has 9 files and the other has 1. than the first one has more than 90% of the results, the second one less than 10%
                //if there are
                if ( ( ( float )( *it ).second / ( float ) fileNumbers.size() ) > 0.50 && ( ( *it ).first != -1 ) )
                {
                    biggestTaxa.push_back( ( *it ).first );
                }
                //}
            }
            //    cout << i << " remaining taxa " << taxaCount.size() <<endl;
            //if this leaves only one remaining taxa use that
            if ( biggestTaxa.size() == 1 && biggestTaxa[0] != -1 )
            {
                taxSame[i] = true;
                sharedTaxIds[i] = biggestTaxa[0];
            }
            // else there are more than one taxa getting more than 10% of the filenumbers
            // so it is assumed that the files have no common taxonomy
            else
            {
                sharedTaxIds[i] = 0;
                taxSame[i] = false;
            }
        }
    }
    // if there was only one file it has to have a common taxa in all taxLevels
    // only the lowest taxlevel needs to be set because only that will be printed
    if ( fileNumbers.size() == 1 )
    {
        for ( unsigned int i( taxLevelSize - 1 ); i > 0; i-- )
        {
            if ( fileNumToTaxIds_[fileNumbers[0]][i] != -1 )
            {
                taxSame[i] = true;
                //      cout << "FileNumbers size 1 " << fileNumbers[0] << " " << fileNumToTaxIds_[fileNumbers[0]][taxLevelSize-1] <<endl;
                sharedTaxIds[i] = fileNumToTaxIds_[fileNumbers[0]][i];
            }
        }
    }
#ifdef MET_DEBUG
    for ( unsigned int i( 0 ); i < fileNumbers.size(); i++ )
        cerr << fileNumbers[i] << " " ;
    cerr << endl;
    for ( unsigned int i( 0 ); i < taxSame.size() ; i++ )
        cerr << sharedTaxIds[i] << " " ;
    cerr << endl;
    for ( unsigned int i( 0 ); i < taxSame.size(); i++ )
        cerr << taxSame[i] << " " ;
    cerr << endl;
#endif
    return taxSame;
}


void IntervalHandlerMetagenome::foundInBoth
( const int pileNum,
  const LetterCount &countsThisRangeA, const LetterCount &countsThisRangeB,
  const Range &thisRangeA, const Range &thisRangeB,
  AlphabetFlag &propagateIntervalA, AlphabetFlag &propagateIntervalB )
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
    for ( int l ( 1 ); l < alphabetSize; l++ )
    {
        //if word is a singleton in the reference files don't propagate anymore
        if ( countsThisRangeB.count_[l] > 0 )
        {
            propagateIntervalB[l] = true;
            const bool propA = countsThisRangeA.count_[l] > minOcc_;
            propagateIntervalA[l] = propA;
            // If k-mer is not propagated because it slowly falls below threshold due to end-of-reads, we set this flag to output it
            if ( !propA && countsThisRangeA.count_[l] + countsThisRangeA.count_[0] >= minOcc_ )
                belowMinDepthBecauseOfEndOfReads = true;
        }
        else
        {
            //if the range in B can't be propagated anymore test if the propagation in A would have been possible
            //if yes print the taxonomic information but don't propagate anymore
            if ( countsThisRangeA.count_[l] >= minOcc_ )
                differentProp = true;
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
    propagateIntervalA[whichPile[( int ) dontKnowChar]] = false;
    //testing if the word is in the reads and in certain genomes
    //only look the fileNumbers up if the word can't be propagatet anymore.
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
        if ( thisRangeB.word_.length() == 35
             || thisRangeB.word_.length() == 50
             || thisRangeB.word_.length() == 100
             || thisRangeB.word_.length() == 150
             || thisRangeB.word_.length() == 200
             || thisRangeB.word_.length() == 250 )
        {
            // get the unique file numbers where this Range can be found
            getFileNumbersForRange( pileNum, thisRangeB.pos_, thisRangeB.num_, fileNumbers );
            // test if word is singleton in the files
            cout << "file numbers size " << fileNumbers.size() << " range " << thisRangeB.num_ << endl;
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
                        //taxaCountPerWordLength[thisRangeB.word_.length()][sharedTaxa[i]] += 1;
                        cout << "BWord" << thisRangeB.word_.length() << "|" << sharedTaxa[i] << "|";
                        break;
                    } //~same taxa
                }// ~ for loop for the tax levels
            }//~test of singletons
        }// ~if wordLengthtest
    }//~testDB


    //if the database should not be tested on its own the normal comparison between the two datasets takes place
    //the suffix of the range can only be once in each file if the amount of suffixes is smaller than all possible files. so start testing if it is smaller
    else if ( ( thisRangeB.num_ < fileNumToTaxIds_.size() )
              //test file numbers and print information only if the suffix is longer than a certain minimal length chosen from the user
              //this reduces the amount of output and speeds up the count word algorithm because the filenumbers don't need to be looked at
              &&  thisRangeB.word_.length() > minWordLength_
              //test file numbers and print taxonomic information only if a possible tax break point or the end of the read is reached
            )
    {

        sharedTaxa.resize( taxLevelSize );
        //put back breackpoint information, if
        const bool printBKPT = true;
        const bool printBKPTdetails = false;
        if ( printBKPT && ( differentProp || belowMinDepthBecauseOfEndOfReads ) )
        {
            printf(
                "BKPT %s %llu:%llu:%llu:%llu:%llu:%llu %llu:%llu:%llu:%llu:%llu:%llu %llu ",
                thisRangeB.word_.c_str(),
                countsThisRangeA.count_[0],
                countsThisRangeA.count_[1],
                countsThisRangeA.count_[2],
                countsThisRangeA.count_[3],
                countsThisRangeA.count_[4],
                countsThisRangeA.count_[5],
                countsThisRangeB.count_[0],
                countsThisRangeB.count_[1],
                countsThisRangeB.count_[2],
                countsThisRangeB.count_[3],
                countsThisRangeB.count_[4],
                countsThisRangeB.count_[5],
                ( thisRangeB.pos_ & matchMask ) );
            if ( printBKPTdetails )
                for ( unsigned int f( 0 ) ; f < fileNumbers.size(); f++ )
                    cout << fileNumbers[f] << ":";
            cout << endl;
        }

        assert( thisRangeA.word_.size() == thisRangeB.word_.size() );
        bool maxLengthReached = ( thisRangeA.word_.size() >= maxWordLength_ + 1 );
        if ( differentProp || belowMinDepthBecauseOfEndOfReads || maxLengthReached )
        {
            //get the unique file number for the files where the suffix can be found in
            getFileNumbersForRange( pileNum, thisRangeB.pos_, thisRangeB.num_, fileNumbers );

            //only get taxonomic information if the word is a singleton in each file
            //unique file numbers must be the same as the numbers of suffixes, so that each file has exactly one suffix
            if ( thisRangeB.num_ <= fileNumbers.size() )
            {
                vector<bool> sameTaxa = intervalInSameTaxa( sharedTaxa, fileNumbers );
                //print the information about a shared superkingdom only if a wordlength higher than 50 is reached
                //printing this information about shorter words inflates the output incredibly without giving much of useful information
                int smallestTaxLevel = ( thisRangeB.word_.length() > 50 ) ? 0 : 1;
                //if there are shared taxonomic ids between the files where this suffix can be found
                //print the information about the deepest tax level
                for ( int i( taxLevelSize - 1 ) ; i >= smallestTaxLevel; i-- )
                {
                    if ( sameTaxa[i] && sharedTaxa[i] != 0 )
                    {
                        //                        if ( !( differentProp || belowMinDepthBecauseOfEndOfReads ) )
                        //                            cout << "ending"; // == This MTAXA is due to maxLengthReached
                        cout << "MTAXA " << i <<  " " << sharedTaxa[i] << " " << thisRangeB.word_ ;
                        cout << " " << ( thisRangeB.pos_ & matchMask ) << " " ;
                        printf( "%llu:%llu:%llu:%llu:%llu:%llu %llu:%llu:%llu:%llu:%llu:%llu ",
                                countsThisRangeA.count_[0],
                                countsThisRangeA.count_[1],
                                countsThisRangeA.count_[2],
                                countsThisRangeA.count_[3],
                                countsThisRangeA.count_[4],
                                countsThisRangeA.count_[5],
                                countsThisRangeB.count_[0],
                                countsThisRangeB.count_[1],
                                countsThisRangeB.count_[2],
                                countsThisRangeB.count_[3],
                                countsThisRangeB.count_[4],
                                countsThisRangeB.count_[5]
                              );
                        for ( unsigned int f( 0 ) ; f < fileNumbers.size(); f++ )
                            cout << fileNumbers[f] << ":";
                        cout << endl;
                        break;
                    }//~if shared taxa
                }//~for loop for the taxonomic levels
            }
        }//~test on singletons in files
    }//~ if matches printing criteria
}

void IntervalHandlerMetagenome::foundInBOnly
( const int pileNum,
  const LetterCount &countsSoFarB,
  const LetterCount &countsThisRangeB,
  const Range &thisRangeB,
  AlphabetFlag &propagateIntervalB )
{
    //this hander is only interesting if you want to get information about the alignability of the database
    if ( testDB_ )
    {
        //propagate as long as possible
        for ( int l( 1 ); l < alphabetSize; l++ )
        {
            propagateIntervalB[l] = countsThisRangeB.count_[l] > 0;
            if ( pileNum == 0 )
                cerr << countsThisRangeB.count_[l] << " " ;
        }
        if ( pileNum == 0 )
            cerr << endl;
        //don't bother about the Ns they are creating to many false positives
        propagateIntervalB[whichPile[( int )dontKnowChar]] = false;
        //only print information about a certain word length to keep the output smaller
        if ( thisRangeB.word_.length() == 35
             || thisRangeB.word_.length() == 50
             || thisRangeB.word_.length() == 100
             || thisRangeB.word_.length() == 150
             || thisRangeB.word_.length() == 200
             || thisRangeB.word_.length() == 250 )
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
                            // print the smalles possible information to keep the output small
                            cout << "BWord" << thisRangeB.word_.length() << "|" << sharedTaxa[i] << "|";
                            break;
                        }
                    }
#ifdef OLD
                    //printing all this created an output bigger than 500GB
                    for ( unsigned int i( 5 ) ; i != 0; i-- )
                    {
                        if ( sameTaxa[i] )
                        {
                            cout << "BTAXA " << taxLevel[i] <<  " " << sharedTaxa[i] << " " << thisRangeB.word_ << " " ;
                            cout << ( thisRangeB.pos_ & matchMask ) << " " << thisRangeB.num_ << endl;
                        }
                    }
                    if ( sameTaxa[taxLevelSize - 1] )
                    {
                        speciesFound = true;

                        printf(
                            "BSPECIES %s %llu %d %llu:%llu:%llu:%llu:%llu:%llu ",
                            thisRangeB.word_.c_str(),
                            thisRangeB.pos_ & matchMask,
                            sharedTaxa[taxLevelSize - 1],
                            countsThisRangeB.count_[0],
                            countsThisRangeB.count_[1],
                            countsThisRangeB.count_[2],
                            countsThisRangeB.count_[3],
                            countsThisRangeB.count_[4],
                            countsThisRangeB.count_[5]
                        );
                        for ( unsigned int j( 0 ) ; j < fileNumbers.size() ; j++ )
                            cout << fileNumbers[j] << ":" ;
                        cout << endl;
                    }
                    if ( fileNumbers.size() == 1 )
                    {
                        printf(
                            "BSINGLE %s %llu %llu:%llu:%llu:%llu:%llu:%llu ",
                            thisRangeB.word_.c_str(),
                            thisRangeB.pos_ & matchMask,
                            countsThisRangeB.count_[0],
                            countsThisRangeB.count_[1],
                            countsThisRangeB.count_[2],
                            countsThisRangeB.count_[3],
                            countsThisRangeB.count_[4],
                            countsThisRangeB.count_[5]
                        );
                        for ( unsigned int j( 0 ) ; j < fileNumbers.size() ; j++ )
                            cout << fileNumbers[j] << ":" ;
                        cout << endl;
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
  const Range &thisRangeA,
  AlphabetFlag &propagateIntervalA )
{
    //this should be interesting if one wants to find high accounted reads in A which have no singletons in the database
    // special interesting for the comparison between metatranscriptomes and metagenomes
    //maybe having a second minOcc, or using the minOcc only here and not in the handler FoundInBoth, because there one wants every occurence
    //  cerr <<"found in A only is that really needed?" <<endl;
    //Is that really interesting?
    //should probably always propagate A if higher than minOcc

    //at the moment propagte es long as possible
    for ( int l( 0 ); l < alphabetSize; l++ )
        propagateIntervalA[l] = false;
}

