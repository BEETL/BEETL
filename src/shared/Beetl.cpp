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

#include <cstdlib>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdio>

#include "Beetl.hh" // own declarations
#include "Algorithm.hh"   // framework class declaration

#include "BWTCollection.h" // interface to BCR
#include "CountWords.hh" // interface to countwords
#include "BCRext.hh"    // interface to BCRext


using namespace std;

#define BEETL_ID "$Id: Beetl.cpp,v 1.14 2011/11/28 16:38:32 acox Exp $"


int main(int numArgs, char** args) {

    if (numArgs < 2) {
        print_usage(args[0]);
    }

    if (strcmp(args[1],COMMAND_BCR) == 0) {

      CompressionFormatType bcrCompression(compressionASCII);
        // start at 2, 0 is the executable, 1 the command name parsed above
        for (int i = 2; i < numArgs; i++)

            if (args[i][0] == '-') { // only flags here "-X etc."

                switch (args[i][1]) {
                    case 'i':
                        // next param should be the filename, checking now
                        isArgumentOrExit(i + 1, numArgs);
                        fileIsReadableOrExit(args[i + 1]);
                        bcrFileIn = args[i + 1]; // this should be the name
                        cout << "-> input file is " << bcrFileIn << endl;
                        break;
                    case 'o':
                        isArgumentOrExit(i + 1, numArgs);
                        bcrFileOut = args[i + 1];
                        cout << "-> output file is " << bcrFileOut << endl;
                        break;
                    case 'm':
                        isArgumentOrExit(i + 1, numArgs);
                        bcrMode = atoi(args[i + 1]);
                        if (bcrMode > 2 || bcrMode < 0) {
                            cerr << bcrMode << " is no valid bcr mode " << endl;
                            exit(-1);
                        }
                        cout << "-> working mode set to \""
                                << bcrModes[bcrMode]
                                << "\""
                                << endl;
                        break;
                    case 'a':
		      bcrCompression=compressionASCII;;
                        cout << "-> writing ASCII encoded output"
                                << endl;
                        break;
                    case 'h':
		      bcrCompression=compressionHuffman;
                        cout << "Huffman encoding not yet supported, sorry."
                                << endl;
                        exit(-1);
                        //cout << "-> writing huffman encoded output"
                        //        << endl;                        
                        break;
                    case 'r':
		      bcrCompression=compressionRunLength;
                        cout << "-> writing runlength encoded output"
                                << endl;
                        break;
                    default:
                        cout << "!! unknown flag \""
                                << args[i][1] << "\"" << endl;
                        print_usage(args[0]);
                }
            }

        // check if all arguments are given
        if (bcrFileOut.length() > 0 && bcrFileIn.length() > 0 && bcrMode >= 0) {
            if (bcrMode == 0){
                isValidFastaFile(bcrFileIn.c_str());
            }

            // created new tool object
            Algorithm * pBCR 
	      = new BCR(bcrMode, bcrFileIn, bcrFileOut, bcrCompression);

            // run previous main method        
            pBCR->run();

            // clean up
            delete pBCR;

            // die
            exit(0);
        } else {
            // something wrong happened
            print_usage(args[0]);
        }

    } else if (strcmp(args[1], COMMAND_BCR_EXT) == 0) {


      // set defaults for BCRext mode
      bcrExtAsciiOutput=false; // use normal ASCII alphabet as output
      bcrExtHuffmanOutput=false; // use huffman encoding as compression
      bcrExtRunlengthOutput=true; // use RunLength encoding [default]
      bcrExtImplicitSort=false; // do implicit sort of input sequences


      for (int i = 2; i < numArgs; i++)
      {
	if (args[i][0] == '-') { // only flags here "-X etc."

	  std::string thisArg((string)args[i]);
	  if (thisArg=="-i")
	  {
	    // next param should be the filename, checking now
	    isArgumentOrExit(i + 1, numArgs);
	    fileIsReadableOrExit(args[i + 1]);
	    bcrExtFileIn = args[i + 1]; // this should be the name
	    cout << "-> input file is " << bcrExtFileIn << endl;
	  }
	  else if (thisArg=="-p")
	  {
	    isArgumentOrExit(i + 1, numArgs);
	    bcrFileOut = args[i + 1];
	    cout << "-> output prefix set to "
		 << bcrFileOut << endl;
	  }
	  else if (thisArg=="-a")
	  {
	    bcrExtAsciiOutput = true;
	    bcrExtRunlengthOutput = false;
	    bcrExtHuffmanOutput = false;
	    cout << "-> writing ASCII encoded output"
		 << endl;
	  }
	  else if (thisArg=="-h")
	  {
	    bcrExtAsciiOutput = false;
	    bcrExtRunlengthOutput = false;
	    bcrExtHuffmanOutput = true;
	    cout << "-> writing huffman encoded output"
		 << endl;                        
	  }
	  else if (thisArg=="-r")
	  {
	    bcrExtAsciiOutput = false;
	    bcrExtRunlengthOutput = true;
	    bcrExtHuffmanOutput = false;
	    cout << "-> writing runlength encoded output"
		 << endl;
	  }
	  else if (thisArg=="-sap")
	  {
	    bcrExtImplicitSort=true;
	    cout << "-> perform implicit sort of input sequences"
		 << endl;
	  }
	  else
	  {
	    cout << "!! unknown flag \""
		 << args[i][1] << "\"" << endl;
	    print_usage(args[0]);
	  }
	} // ~if begins with -
      } // ~for

      cout << bcrExtImplicitSort << bcrExtHuffmanOutput << bcrExtRunlengthOutput << endl;
      if ((bcrExtImplicitSort==true)&&
	  ((bcrExtHuffmanOutput==true)||(bcrExtRunlengthOutput==true)))
      {
	cout << "-> Note: -sap mode needs ASCII intermediate files," << endl
	     << "-> will revert to requested compression type for final output" << endl;

      }

        // check if all arguments are given
        if ((bcrExtRunlengthOutput || bcrExtAsciiOutput || bcrExtHuffmanOutput) // no huffman for now
             && isValidReadFile(bcrExtFileIn.c_str())) {
            if (bcrExtFileOutPrefix.length()==0){
                bcrExtFileOutPrefix=bcrExtFileOutPrefixDefault;
            }

            // created new tool object
            Algorithm * pBCRext = new BCRext(bcrExtHuffmanOutput,
                                             bcrExtRunlengthOutput,
                                             bcrExtAsciiOutput,
					     bcrExtImplicitSort,
                                             bcrExtFileIn,
                                             bcrExtFileOutPrefix);

            // run previous main method        
            pBCRext->run();

            // clean up
            delete pBCRext;

            // die
            exit(0);
        } else {
            // something wrong happened
            print_usage(args[0]);
        }

    } else if (strcmp(args[1], COMMAND_COUNTWORDS) == 0) {


      vector<string> filesA, filesB;

        for (int i = 2; i < numArgs; i++)

            if (args[i][0] == '-') { // only flags here "-X etc."
	      cout << "blah" << endl;
                switch (args[i][1]) {
                    case 'a':
		      while (args[++i][0]!='-')
		      {
			cout << args[i] << " fred " << endl;
                        fileIsReadableOrExit(args[i]);
			filesA.push_back(args[i]);
                        cout << "-> input file A is "
			     << filesA.back()
			     << endl;
		      }
		      i--;
		      break;
#ifdef OLD
                        // next param should be the filename, checking 
                        isArgumentOrExit(i + 1, numArgs);
                        fileIsReadableOrExit(args[i + 1]);
                        countWordsInputA = args[i + 1]; // should be the name
                        cout << "-> input file A is "
                                << countWordsInputA
                                << endl;
                        break;
#endif
                    case 'b':
		      while ((++i)!=numArgs)
		      {
                        fileIsReadableOrExit(args[i]);
			filesB.push_back(args[i]);
                        cout << "-> input file B is "
			     << filesB.back()
			     << endl;
		      }
		      //		      i--;
		      break;
#ifdef OLD
                        // next param should be the filename, checking now
                        isArgumentOrExit(i + 1, numArgs);
                        fileIsReadableOrExit(args[i + 1]);
                        countWordsInputB = args[i + 1]; // should be the name
                        cout << "-> input file B is "
                                << countWordsInputB
                                << endl;
                        break;
#endif
                    case 'k':
                        isArgumentOrExit(i + 1, numArgs);
                        minimalLengthK = atoi(args[i + 1]);
                        if (minimalLengthK < 0) {
                            cerr << "!! "
                                    << minimalLengthK
                                    << " is no valid length "
                                    << endl;
                            exit(-1);
                        }
                        cout << "-> minimal length k set to \""
                                << minimalLengthK
                                << "\""
                                << endl;
                        break;
                    case 'n':
                        isArgumentOrExit(i + 1, numArgs);
                        minimalOccurencesN = atoi(args[i + 1]);
                        if (minimalOccurencesN < 0) {
                            cerr << "!! "
                                    << minimalOccurencesN
                                    << " is no valid value "
                                    << endl;
                            exit(-1);
                        }
                        cout << "-> maximal occurences n set to \""
                                << minimalOccurencesN
                                << "\""
                                << endl;
                        break;
                    case 'r':
                        cout << "-> reference genome mode for set B " << endl;
                        ReferenceGenomeInputB = true;
                        break;
                    case 'A':
                        cout << "-> assuming set A is compressed" << endl;
                        compressedInputA = true;
                        break;
                    case 'B':
                        cout << "-> assuming set B is compressed" << endl;
                        compressedInputB = true;
                        break;
                    case 'C':
                        cout << "-> assuming set A&B are compressed" << endl;
                        compressedBoth = true;
                        break;
                    default:
                        cout << "!! unknown flag \"" << args[i][1]
                                << "\"" << endl;
                        print_usage(args[0]);
                }
            }

        // check for required arguments
        if ( (minimalLengthK>0) && (minimalOccurencesN>0) &&
	     (filesA.size()>0) && (filesA.size()==filesB.size())) {

            // create new tool object
            Algorithm * pcountWords = new countWords(compressedBoth, compressedInputA,
                    compressedInputB, ReferenceGenomeInputB, minimalOccurencesN,
                    minimalLengthK, filesA, filesB);
            
            // run the "main" method
            pcountWords->run();
            
            // clean up
            delete pcountWords;
            
            // closing time
            exit(0);
        } else {
            // oops
            print_usage(args[0]);
        }

    } else {
        cerr << "!! \"" << args[1] << "\" is no known command" << endl;
        print_usage(args[0]);
    }
    return 0;
}

void fileIsReadableOrExit(string filename) {

    FILE * pFile;
    // test file for read access
    pFile = fopen(filename.c_str(), "r");

    if (pFile != NULL) {
        fclose(pFile);
        return;
    } else {
        cerr << "!! \"" << filename << "\" is NOT readable!" << endl;
        exit(-1);
    }
}

void isArgumentOrExit(int num, int numArgs) {
    if ((num) > (numArgs - 1)) {
        cerr << "!! CLI parsing error. Wrong number of arguments?" << endl;
        exit(-1);
    }
}


void print_usage(char *args) {
    cerr << endl << "- This is the BEETL software library -" << endl
            << endl
      // Tony 13.6.12 - BEETL_ID is not informative now we have moved to git
      //            << "Framework version" << endl 
      //            << BEETL_ID << endl            
            << endl
            << "Included in this framework are the following algorithms" << endl
            << endl
            << endl
            << "-> BCRext - command \"" << COMMAND_BCR_EXT << "\"" << endl
            << "========================================================" << endl
            << "improved version of the original algorithm" << endl
            << "uses significantly less RAM (a.k.a. none) but depends heavily on I/O" << endl
            << endl
            << "Usage: " << args << " "
            << COMMAND_BCR_EXT <<" -i <read file> -p <output file prefix> [-r -a] [sap]" << endl
            // below: for huffman encoding uncomment when implemented
            //<< COMMAND_BCR_EXT <<" -i <read file> -p <output file prefix> [-h -r -a]" << endl
            << endl
            << "-i <file>:\tinput set of , 1 read per line, no fasta" << endl
            << "-p <string>:\toutput file names will start with \"prefix\"" << endl
            << "-a:\t\toutput ASCII encoded files" << endl
            << "-r:\t\toutput runlength encoded files [default]" << endl
      //            << "-h:\t\toutput Hufmann encoded files" << endl
            << "-sap:\t\tperform implicit permutation of collection to obtain more compressible BWT"  
            << endl
            << endl
            << "-> BCR - command \"" << COMMAND_BCR << "\"" << endl
            << "========================================================" << endl
            << "original algorithm to construct the BWT of a set of reads" << endl
            << "needs approximately 14GB of RAM for 1 billion reads" << endl
            << endl
            << "Usage: " << args << " " 
            << COMMAND_BCR <<" -i <fasta read file> -o <output file> -m <[0,1,2]>" << endl
            << endl
            << "-i <file>:\tinput set of reads" << endl
            << "-o <file>:\toutput file" << endl
            << "-m <n>:\t\tmode = 0 --> BCR " << endl
            << "\t\tmode = 1 --> unBCR " << endl
            << "\t\tmode = 2 --> Backward search + Locate SeqID " << endl
            << endl
            << endl
            << "-> countWords - command \"" << COMMAND_COUNTWORDS << "\"" << endl
            << "========================================================" << endl
            << "find all words of length at least k that occur" << endl
            << "at least n times in string set A and never in string set B" << endl
            << endl
            << "Usage: " << args << " " 
            << COMMAND_COUNTWORDS <<" [-A -B -C -r] -k <n> -n <n> -a <set A> -b <set B>" << endl
            << endl
            << "-A:\t\tassume BWT files for set A are in compressed format" << endl
            << "-B:\t\tassume BWT files for set B are in compressed format" << endl
            << "-C:\t\tassume BWT files for sets A and B are in compressed format" << endl
            << "-r:\t\tassume set B is a reference genome" << endl
            << "-k <n>:\t\tminimal length" << endl
            << "-n <n>:\t\tminimal occurences" << endl
            << "-a <file>:\tinput set A" << endl
            << "-b <file>:\tinput set B" << endl
            << endl
            << endl
            << "If you had fun using these algorithms you may cite:" << endl
            << "---------------------------------------------------" << endl          
            << "Markus J. Bauer, Anthony J. Cox and Giovanna Rosone" << endl
            << "Lightweight BWT Construction for Very Large String Collections. " << endl
	 << "Proceedings of CPM 2011, pp.219-231, doi: 10.1007/978-3-642-21458-5_20" << endl 
	 << "[Description of BWT construction algorithms in 'bcr' and 'ext' modes]" << endl << endl
	 << "Markus J. Bauer, Anthony J. Cox and Giovanna Rosone" << endl
	 << "Lightweight algorithms for constructing and inverting the BWT of string collections " 
	 << endl << "Theoretical Computer Science, doi: 10.1016/j.tcs.2012.02.002" << endl 
	 << "[As above plus description of BWT inversion in 'bcr' mode]" << endl << endl
	 << "Anthony J. Cox, Markus J. Bauer, Tobias Jakobi and Giovanna Rosone" << endl
	 << "Large-scale compression of genomic sequence databases with the Burrows-Wheeler transform"
	 << endl
	 << "Bioinformatics, doi: 10.1093/bioinformatics/bts173" << endl 
	 << "[Description of '-sap' compression boosting strategy in 'ext' mode]" << endl << endl
	 << "BEETL web page:" << endl
	 << "---------------" << endl
	 << "http://beetl.github.com/BEETL/" << endl << endl	 
	 << "BEETL user group:" << endl
	 << "-----------------" << endl
	 << "http://tech.groups.yahoo.com/group/BEETL/" << endl << endl;

    exit(0);
}
