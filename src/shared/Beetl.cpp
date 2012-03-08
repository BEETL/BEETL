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

        for (int i = 2; i < numArgs; i++)

            if (args[i][0] == '-') { // only flags here "-X etc."

                switch (args[i][1]) {
                    case 'i':
                        // next param should be the filename, checking now
                        isArgumentOrExit(i + 1, numArgs);
                        fileIsReadableOrExit(args[i + 1]);
                        bcrExtFileIn = args[i + 1]; // this should be the name
                        cout << "-> input file is " << bcrExtFileIn << endl;
                        break;
                    case 'p':
                        isArgumentOrExit(i + 1, numArgs);
                        bcrFileOut = args[i + 1];
                        cout << "-> output prefix set to "
                                << bcrFileOut << endl;
                        break;
                    case 'a':
                        if (!bcrExtRunlengthOutput && !bcrExtHuffmanOutput) {
                            bcrExtAsciiOutput = true;
                        }
                        cout << "-> writing ASCII encoded output"
                                << endl;
                        break;
                    case 'h':
                        if (!bcrExtRunlengthOutput && !bcrExtAsciiOutput) {
                            bcrExtHuffmanOutput = true;
                        }
                        cout << "-> writing huffman encoded output"
                                << endl;                        
                        break;
                    case 'r':
                        if (!bcrExtHuffmanOutput && !bcrExtAsciiOutput) {
                            bcrExtRunlengthOutput = true;
                        }
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
        if ((bcrExtRunlengthOutput || bcrExtAsciiOutput || bcrExtHuffmanOutput) // no huffman for now
             && isValidReadFile(bcrExtFileIn.c_str())) {
            if (bcrExtFileOutPrefix.length()==0){
                bcrExtFileOutPrefix=bcrExtFileOutPrefixDefault;
            }

            // created new tool object
            Algorithm * pBCRext = new BCRext(bcrExtHuffmanOutput,
                                             bcrExtRunlengthOutput,
                                             bcrExtAsciiOutput,
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

    } 
	else {
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
            << "Framework version" << endl 
            << BEETL_ID << endl            
            << endl
            << "Included in this framework are the following algorithms" << endl
            << endl
            << endl
            << "-> BCRext - command \"" << COMMAND_BCR_EXT << "\"" << endl
            << "========================================================" << endl
            << "improved version of the orginal algorithm" << endl
            << "uses significantly less RAM (a.k.a. none) but depends heavily on I/O" << endl
            << endl
            << "Usage: " << args << " "
            << COMMAND_BCR_EXT <<" -i <read file> -p <output file prefix> [-r -a]" << endl
            // below: for huffman encoding uncomment when implemented
            //<< COMMAND_BCR_EXT <<" -i <read file> -p <output file prefix> [-h -r -a]" << endl
            << endl
            << "-i <file>:\tinput set of , 1 read per line, no fasta" << endl
            << "-p <string>:\toutput file names will start with \"prefix\"" << endl
            << "-a:\t\toutput ASCII encoded files" << endl
            << "-r:\t\toutput runlength encoded files [default]" << endl
            << "-h:\t\toutput Hufmann encoded files" << endl
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
            << "If you had fun using these algorithms you may cite:" << endl
            << "---------------------------------------------------" << endl          
            << "Markus J. Bauer, Anthony J. Cox and Giovanna Rosone" << endl
            << "Lightweight BWT Construction for Very Large String Collections. " << endl
            << "Proceedings of CPM 2011, pp.219-231" << endl;
    exit(0);
}
