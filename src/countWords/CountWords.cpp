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

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "Types.hh"
#include "Config.hh"
#include "Timer.hh"

#include "LetterCount.hh"
#include "ReadBuffer.hh"
#include "BwtWriter.hh"
#include "BwtReader.hh"

#include "CountWords.hh"



#ifdef USE_POSIX_FILE_OPTIMIZATIONS
#define _XOPEN_SOURCE 600
#endif
#include <fcntl.h>

//#define DEBUG 1

#define COUNTWORDS_ID "@(#) $Id: CountWords.cpp,v 1.3 2011/11/17 15:52:37 tjakobi Exp $"

using namespace std;

countWords::countWords(bool bothSetsCompressed, bool inputACompressed,
        bool inputBCompressed, bool inputBIsReferenceGenome,
		       int paramN, int paramK, const vector<string>& setA, 
		       const vector<string>& setB) :
  setA_(setA), setB_(setB)
{

#ifdef XXX
    // get memory allocated
    setA_ = new char[setA.length()+1];
    setB_ = new char[setB.length()+1];

    // copt to char * in order to get valid c strings
    setA.copy(setA_, setA.length());
    setB.copy(setB_, setB.length());
    
    // append \0 to obtain a valid escaped c string
    setA_[setA.length()] = '\0';
    setB_[setB.length()] = '\0';
#endif

    // set tool flags
    bothSetsCompressed_ = bothSetsCompressed;
    inputACompressed_ = inputACompressed;
    inputBCompressed_ = inputBCompressed;
    inputBIsReferenceGenome_ = inputBIsReferenceGenome;
    paramN_ = paramN;
    paramK_ = paramK;
}

//
// IntervalHandlerSplice member function declarations
//

void IntervalHandlerSplice::foundInBoth
( const LetterCount& countsThisRangeA, const LetterCount& countsThisRangeB,
  const Range& thisRangeA, const Range& thisRangeB,
  AlphabetFlag& propagateIntervalA, AlphabetFlag& propagateIntervalB)
{
  bool sharedPath(false);
  LetterCountType maxSignalAOnly(0), maxSignalBOnly(0);

  for (int l(1);l<alphabetSize;l++) 			
  {
    if ((countsThisRangeB.count_[l]==0)
	&&(countsThisRangeA.count_[l]>maxSignalAOnly))
      maxSignalAOnly=countsThisRangeA.count_[l];

    if ((countsThisRangeA.count_[l]==0)
	&&(countsThisRangeB.count_[l]>maxSignalBOnly))
      maxSignalBOnly=countsThisRangeB.count_[l];
	
    sharedPath|=((countsThisRangeA.count_[l]>0)&&(countsThisRangeB.count_[l]>0));

    propagateIntervalA[l]=(countsThisRangeA.count_[l]>=minOcc_);
    propagateIntervalB[l]=(countsThisRangeB.count_[l]>0);
  } // ~for l
  if ((sharedPath==false)
      &&(maxSignalAOnly>=minOcc_)
      &&(maxSignalBOnly>=minOcc_))
  {
    cout << "BKPT " << thisRangeA.word_;
    for (int l(0);l<alphabetSize;l++) 
      cout << ((l==0)?" ":":") << countsThisRangeA.count_[l];
    for (int l(0);l<alphabetSize;l++) 
      cout << ((l==0)?" ":":") << countsThisRangeB.count_[l];	
    cout << endl;
  }

  // don't bother with Ns
  propagateIntervalA[whichPile[(int)dontKnowChar]]=false;
  propagateIntervalB[whichPile[(int)dontKnowChar]]=false;

} // ~foundInBoth



void IntervalHandlerSplice::foundInAOnly
( const LetterCount& countsThisRangeA,
  const Range& thisRangeA,
  AlphabetFlag& propagateIntervalA )
{
  if (countsThisRangeA.count_[0]>0)
  {
    cout << "READ " << thisRangeA.word_;
    cout << " " << thisRangeA.pos_;
    for (int l(0);l<alphabetSize;l++) 
      cout << ((l==0)?" ":":") << countsThisRangeA.count_[l];
    cout << endl;
  }
  // TBD print out IDs of discovered reads

  for (int l(1);l<alphabetSize;l++) 			
  {
    propagateIntervalA[l]=(countsThisRangeA.count_[l]>0);
  } // ~for l

  // don't bother with Ns
  propagateIntervalA[whichPile[(int)dontKnowChar]]=false;
} // ~foundInBoth

//
// IntervalHandlerReference member function declarations
//

void IntervalHandlerReference::foundInBoth
( const LetterCount& countsThisRangeA, const LetterCount& countsThisRangeB,
  const Range& thisRangeA, const Range& thisRangeB,
  AlphabetFlag& propagateIntervalA, AlphabetFlag& propagateIntervalB)
{
  bool significantNonRef(false);
  //  LetterCountType maxSignalAOnly(0), maxSignalBOnly(0);

  if (thisRangeB.num_>1)
  { // if k-mer is not unique in B (reference), 
    // propagate all B and all A that matches B  
    for (int l(1);l<alphabetSize;l++) 			
    {
      if (countsThisRangeB.count_[l]>0)
      {
	propagateIntervalB[l]=true;
	if (countsThisRangeA.count_[l]>0)
	  propagateIntervalA[l]=true;
	else
	  propagateIntervalA[l]=false;
      } // ~if
      else
      {
	  propagateIntervalA[l]=false;
	  propagateIntervalB[l]=false;
      } // ~else
    } // ~for
  } // ~if
  else
  {
    if (thisRangeB.num_!=1)
    {
        cerr << "thisRandB is no 1. Aborting." << endl;
    }
    for (int l(1);l<alphabetSize;l++) 			
    {
      if (countsThisRangeB.count_[l]>0)
      {
	propagateIntervalB[l]=true;
	if (countsThisRangeA.count_[l]>0)
	  propagateIntervalA[l]=true;
	else
	  propagateIntervalA[l]=false;
      } // ~if
      else
      {
	  propagateIntervalB[l]=false;
	  if (countsThisRangeA.count_[l]>minOcc_)
	  {
	    propagateIntervalA[l]=true;
	    significantNonRef=true;
	  }
	  else
	    propagateIntervalA[l]=false;
      } // ~else
    } // ~for
  

  } // ~else

  if (significantNonRef==true)
  {
    cout << "BKPT " << thisRangeB.word_;
    cout << " " << thisRangeB.pos_;
    for (int l(0);l<alphabetSize;l++) 
      cout << ((l==0)?" ":":") << countsThisRangeA.count_[l];
    for (int l(0);l<alphabetSize;l++) 
      cout << ((l==0)?" ":":") << countsThisRangeB.count_[l];	
    cout << endl;
  }

  // don't bother with Ns
  propagateIntervalA[whichPile[(int)dontKnowChar]]=false;
  propagateIntervalB[whichPile[(int)dontKnowChar]]=false;

} // ~foundInBoth

#ifdef OLD
{
  for (int l(1);l<alphabetSize;l++) 			
  {
    propagateIntervalA[l]=(countsThisRangeA.count_[l]>=minOcc_);
    propagateIntervalB[l]=(countsThisRangeB.count_[l]>0);
  } // ~for l
} // ~foundInBoth
#endif

void IntervalHandlerReference::foundInAOnly
( const LetterCount& countsThisRangeA,
  const Range& thisRangeA,
  AlphabetFlag& propagateIntervalA )
{

  bool significantPath(false);

  for (int l(1);l<alphabetSize;l++) 			
  {
    if (countsThisRangeA.count_[l]>=minOcc_)
    {
      significantPath=true;
      propagateIntervalA[l]=true;
    } // ~if
    else
    {
      propagateIntervalA[l]=false;
    } // ~else

  } // ~for l

  if (significantPath==false)
  {
    cout << "READ " << thisRangeA.word_;
    cout << " " << thisRangeA.pos_;
    for (int l(0);l<alphabetSize;l++) 
      cout << ((l==0)?" ":":") << countsThisRangeA.count_[l];
    cout << endl;
  }


#ifdef OLD
  // For now this is same as for Splice - continue until all reads found
  if (countsThisRangeA.count_[0]>0)
  {
    cout << "READ " << thisRangeA.word_;
    cout << " " << thisRangeA.pos_;
    for (int l(0);l<alphabetSize;l++) 
      cout << ((l==0)?" ":":") << countsThisRangeA.count_[l];
    cout << endl;
  }
  // TBD print out IDs of discovered reads

  for (int l(1);l<alphabetSize;l++) 			
  {
    propagateIntervalA[l]=(countsThisRangeA.count_[l]>0);
  } // ~for l
#endif

  // don't bother with Ns
  propagateIntervalA[whichPile[(int)dontKnowChar]]=false;
} // ~foundInBoth

//
// BackTracker member function declarations
//


void BackTracker::skipIfNecessary( const Range& thisRange, 
		      LetterCountType& currentPos,
		      BwtReaderBase& inBwt,
		      LetterCount& countsSoFar )
{
#ifdef DEBUG
  cout << "Want 2 skip " << thisRange.pos_-currentPos << ": " << currentPos << " to " << thisRange.pos_<< endl;
#endif
  if (thisRange.pos_>currentPos)
  {	
#ifdef DEBUG
    cout << "Skipping " << thisRange.pos_-currentPos << ": " << currentPos << " to " << thisRange.pos_<< endl;
#endif
    inBwt.readAndCount(countsSoFar,thisRange.pos_-currentPos);
    currentPos=thisRange.pos_;
  } // ~if	
  if (thisRange.pos_<currentPos){
      cerr << "thisRange is to low. Should be > " << currentPos << "." <<endl;
  }
} // ~BackTracker::skipIfNecessary 


template<class IntervalHandlerParam> 
void BackTracker::operator()
( int i, string& thisWord, IntervalHandlerParam& intervalHandler_  )
{
  LetterCount countsThisRangeA,	countsThisRangeB;	
  Range thisRangeA,thisRangeB;
  //  string thisWord;
  bool hasChild;
	bool notAtLastB(rB.getRange(thisRangeB)); 
	while (rA.getRange(thisRangeA))
	{
	  while ((notAtLastB)&&(thisRangeB.word_<thisRangeA.word_))
	    notAtLastB=rB.getRange(thisRangeB);
	  if ((notAtLastB)&&(thisRangeB.word_==thisRangeA.word_))
	  { // A-range has a counterpart B-range, propagate them both 
#ifdef DEBUG
	    cout << "RangeA: " << i << " " << j << " " 
	  	<< thisRangeA.word_ << " " << thisRangeA.pos_ << " " << thisRangeA.num_
		 << " -- " << currentPosA << " = ";
	    cout << "RangeB: " << i << " " << j << " " 
	  	<< thisRangeB.word_ << " " << thisRangeB.pos_ << " " << thisRangeB.num_
		 << " -- " << currentPosB << endl;
#endif


	    skipIfNecessary(thisRangeA,currentPosA,(*inBwtA),countsSoFarA);
	    // count children
	    countsThisRangeA.clear();
	    inBwtA->readAndCount(countsThisRangeA,thisRangeA.num_);
#ifdef DEBUG
	    countsThisRangeA.print();
#endif

	    skipIfNecessary(thisRangeB,currentPosB,(*inBwtB),countsSoFarB);
	    // count children
	    countsThisRangeB.clear();
	    inBwtB->readAndCount(countsThisRangeB,thisRangeB.num_);
#ifdef DEBUG
	    countsThisRangeB.print();
#endif

	    intervalHandler_.foundInBoth
	      ( countsThisRangeA, countsThisRangeB,
		thisRangeA, thisRangeB,
		propagateIntervalA_, propagateIntervalB_);


	    hasChild=false;
	    for (int l(1);l<alphabetSize;l++) 			
	    {
	      //	      assert ((propagateIntervalA_[l]==true)==(countsThisRangeA.count_[l]>=minOcc));
	      //  assert ((propagateIntervalB_[l]==true)==(countsThisRangeB.count_[l]>0));
	      //  if (countsThisRangeA.count_[l]>=minOcc)
	      if (propagateIntervalA_[l]==true)
	      {
		if (hasChild==false)
		{
		  // assert(thisWord.size()==thisRangeA.word_.size()+1);
		  thisWord.replace(1,thisRangeA.word_.size(),thisRangeA.word_);
		  hasChild=true;
		}
		thisWord[0]=alphabet[l];
		//		thisWord+=thisRangeA.word_;
		rA.addRange(l,i,thisWord,
			    countsSoFarA.count_[l],
			    countsThisRangeA.count_[l]);
		// if (countsThisRangeB.count_[l]>0)
		if (propagateIntervalB_[l]==true)
		  rB.addRange(l,i,thisWord,
			    countsSoFarB.count_[l],
			    countsThisRangeB.count_[l]);
	      } // ~if
	    } // ~for 

	    countsSoFarA+=countsThisRangeA;
	    currentPosA+=thisRangeA.num_;
	    countsSoFarB+=countsThisRangeB;
	    currentPosB+=thisRangeB.num_;

	    numRanges++;
	  } // ~if
	  else
	  { // A-range has no counterpart B-range, propagate while count threshold continues to be met
#ifdef DEBUG
	    cout << "RangeA: " << i << " " << j << " " 
	  	<< thisRangeA.word_ << " " << thisRangeA.pos_ << " " << thisRangeA.num_
		 << " -- " << currentPosA << " < " << thisRangeB.word_ << endl;
#endif

	    skipIfNecessary(thisRangeA,currentPosA,(*inBwtA),countsSoFarA);
	    // count children
	    countsThisRangeA.clear();
	    inBwtA->readAndCount(countsThisRangeA,thisRangeA.num_);
#ifdef DEBUG
	    countsThisRangeA.print();
#endif

	    intervalHandler_.foundInAOnly
	      ( countsThisRangeA,
		thisRangeA,
		propagateIntervalA_ );

	    
	    // add ranges for any children
 	    hasChild=false;
	    for (int l(1);l<alphabetSize;l++) 			
	    {
	      //	      if (countsThisRangeA.count_[l]>=minOcc)
	      if (propagateIntervalA_[l]==true)
	      {
		if (hasChild==false)
		{
		  // assert(thisWord.size()==thisRangeA.word_.size()+1);
		  thisWord.replace(1,thisRangeA.word_.size(),thisRangeA.word_);
		  hasChild=true;
		} // ~if
		thisWord[0]=alphabet[l];

		//		hasChild=true;
		//		thisWord=alphabet[l];
		//	thisWord+=thisRangeA.word_;
		rA.addRange(l,i,thisWord,
			    countsSoFarA.count_[l],
			    countsThisRangeA.count_[l]);
	      } // ~if
	    } // ~for l
	    if (hasChild==false)
	    { //  if no children, print word itself
	      cout << "GOLD " << thisRangeA.word_ << " " << thisRangeA.num_ << endl;
	      numSingletonRanges++;
	    } // ~if

	    countsSoFarA+=countsThisRangeA;
	    currentPosA+=thisRangeA.num_;

	    numRanges++;
	  } // ~else
	} // ~while
} // ~BackTracker::operator()

// by Tobias:
// replaces previous main method, params mapped by constructor to member vars

void countWords::run(void) {
  Timer  timer;


  // input streams for BWT from previous iter - 1 per pile
  // Used to compute the counts to do backward search
  vector <BwtReaderBase*> inBwtA(alphabetSize);
  vector <BwtReaderBase*> inBwtB(alphabetSize);



  
  //  RangeStoreRAM rA,rB;
  RangeStoreExternal rA("count_A_1","count_A_2"),rB("count_B_1","count_B_2");


  LetterCountEachPile countsPerPileA, countsCumulativeA;
  LetterCountEachPile countsPerPileB, countsCumulativeB;



  bool compressedInputA(false), compressedInputB(false);
  bool referenceMode(false); // if false, use splice junction mode

  // Tobias:
  // all params now directly set through constructor, no need to do CLI parsing
  // here anymore
  
//  while(1)
//  {
//    if (nextArg==numArgs) break;
//    else if (strcmp(args[nextArg],"-c")==0)
//    {
//      cerr << "Assuming BWT files are in compressed form" << endl;
//      compressedInputA=true;
//      compressedInputB=true;
//      nextArg++;
//    }
//    else if (strcmp(args[nextArg],"-cA")==0)
//    {
//      cerr << "Assuming BWT files for set A are in compressed form" << endl;
//      compressedInputA=true;
//      nextArg++;
//    }
//    else if (strcmp(args[nextArg],"-cB")==0)
//    {
//      cerr << "Assuming BWT files for set B are in compressed form" << endl;
//      compressedInputB=true;
//      nextArg++;
//    }
//    else if (strcmp(args[nextArg],"-ref")==0)
//    {
//      cerr << "Assuming BWT files for set B are from reference genome" << endl;
//      referenceMode=true;
//      nextArg++;
//    }
//    else break;
//  }

    referenceMode = inputBIsReferenceGenome_;
    compressedInputA = inputACompressed_;
    compressedInputB = inputBCompressed_;

    if (bothSetsCompressed_) {
        compressedInputA = true;
        compressedInputB = true;
    }

    int numCycles = paramK_;
    int minOcc = paramN_;

    cout << numCycles << " " << minOcc << " fred" << endl;
  
  // see above, needed for old own main method

//  if (nextArg!=numArgs)
//  {
//    numCycles=(atoi(args[nextArg++]));
//    if (nextArg!=numArgs)
//      minOcc=(atoi(args[nextArg++]));
//  } // ~if
//  
//  if ((numCycles==-1)||(minOcc==-1)||(numArgs-nextArg!=(2*alphabetSize)))
//  {
//    cerr << args[0] << ": find all words of length at least k that occur" 
//	 << endl
//	 << "at least n times in string set A and never in string set B"
//	 << endl << endl
//	 << "Usage: " << args[0] 
//	 << " [-c -cA -cB -ref] k n <BWT files for set A> <BWT files for set B>" << endl
//	 << "-cA:  assume BWT files for set A are in compressed format" << endl 
//	 << "-cB:  assume BWT files for set B are in compressed format" << endl 
//	 << "-c:   assume BWT files for sets A and B are in compressed format" << endl 
//	 << "-ref: assume set B is a reference genome" << endl 
//	 << endl;
//    assert(1==0);
//  }  

      //  int numCycles(atoi(args[nextArg++]));
      //  int minOcc(atoi(args[nextArg++]));

    cout << setA_.size() << " " << setB_.size() << endl;

  for (int i(0);i<alphabetSize;i++)
  {
    //#ifdef DEBUGX
    cout << i << " " << alphabet[i] << endl; 
    //#endif
    cout << setA_[i] << " " << setB_[i] << endl;
    if ((compressedInputA==true)&&(i!=0))
      inBwtA[i]= new BwtReaderRunLength(setA_[i]);
    else
      inBwtA[i]= new BwtReaderASCII(setA_[i]);
    inBwtA[i]->readAndCount(countsPerPileA[i]);

    if ((compressedInputB==true)&&(i!=0))
      inBwtB[i]= new BwtReaderRunLength(setB_[i]);
    else
      inBwtB[i]= new BwtReaderASCII(setB_[i]);
    //    inBwtB[i]= new BwtReaderASCII(args[i+alphabetSize+nextArg]);
    inBwtB[i]->readAndCount(countsPerPileB[i]);

#ifdef DEBUG
    countsPerPileA[i].print();
    countsPerPileB[i].print();
#endif
  }

  countsCumulativeA=countsPerPileA;
  countsCumulativeB=countsPerPileB;
  for (int i(1);i<alphabetSize;i++)
  {
    countsCumulativeA[i]+=countsCumulativeA[i-1];
    countsCumulativeB[i]+=countsCumulativeB[i-1];
  }

#ifdef DEBUG
  countsCumulativeA.print();
  countsCumulativeB.print();
#endif

  string thisWord;
  const int dontKnowIndex(whichPile[(int)dontKnowChar]);

  // sort out first iter
  for (int i(1);i<alphabetSize;++i)
  {
    for (int j(1);j<alphabetSize;++j)
    {  
#ifdef DEBUG
      cout << i << " " << j << endl;
#endif
      thisWord.clear();thisWord+=alphabet[j];thisWord+=alphabet[i];

      if ((i!=dontKnowIndex)&&(j!=dontKnowIndex))
      { // get rid of any ranges with N in them 
	if (countsPerPileA[i].count_[j]!=0)
	  rA.addRange(j,i,thisWord,
		      countsCumulativeA[i-1].count_[j],
		      countsPerPileA[i].count_[j]);
	if (countsPerPileB[i].count_[j]!=0)
	  rB.addRange(j,i,thisWord,
		      countsCumulativeB[i-1].count_[j],
		      countsPerPileB[i].count_[j]);
      } // ~if
    } // ~for j
  } // ~for i


  LetterCount countsSoFarA, countsSoFarB;	
  LetterCountType currentPosA, currentPosB;
  LetterCountType numRanges,numSingletonRanges;	
  //  Range thisRangeA,thisRangeB;
  rA.clear(); rB.clear(); 
  for (int c(0);c<numCycles;++c)
  { 
    cerr << "Starting iteration " << c << ", time now: " << timer.timeNow();
    cerr<< "Starting iteration " << c << ", usage: " << timer << endl;

    thisWord.resize(c+3);

    numRanges=0; numSingletonRanges=0;

    rA.swap(); rB.swap();

    


//    pTemp=pNext;pNext=pThis;pThis=pTemp;
    for (int i(1);i<alphabetSize;++i)
    {
      inBwtA[i]->rewindFile(); inBwtB[i]->rewindFile();
      currentPosA=0; currentPosB=0; 
      countsSoFarA.clear();countsSoFarA+=countsCumulativeA[i-1];
      countsSoFarB.clear();countsSoFarB+=countsCumulativeB[i-1];
      for (int j(1);j<alphabetSize;++j)
      {  
#ifdef DEBUG
	cout << "positions - A: " << currentPosA << " B: "
	     << currentPosB << endl;
#endif
	rA.setPortion(i,j);
	rB.setPortion(i,j);

	BackTracker backTracker( inBwtA[i],inBwtB[i],
				 currentPosA,currentPosB,
				 rA,rB,countsSoFarA,countsSoFarB,minOcc);
	if (referenceMode==true)
	{
	  IntervalHandlerReference intervalHandler(minOcc);
	  backTracker(i,thisWord,intervalHandler);
	}
	else
	{
	  IntervalHandlerSplice intervalHandler(minOcc);
	  backTracker(i,thisWord,intervalHandler);
	}
	numRanges+=backTracker.numRanges;
	numSingletonRanges+=backTracker.numSingletonRanges;


      } // ~for j
    }     // ~for i

    rA.clear(); rB.clear();


    cerr << "Finished cycle " << c << ": ranges=" << numRanges
	 << " singletons=" << numSingletonRanges 
	 << " usage: " << timer << endl;

    //    return 0; // %%%
    if (numRanges==0) break;
  } // ~for c

} // ~main


