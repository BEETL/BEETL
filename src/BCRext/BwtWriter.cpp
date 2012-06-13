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

#include <cstring>
#include <cstdlib>

#include "LetterCount.hh"
#include "Tools.h"
#include "BwtWriter.hh"

//#define DEBUG 1
using namespace std;


//
// BwtWriterBase member function definitions
//


BwtWriterBase::BwtWriterBase( const std::string& fileName ) : pFile_(fopen(fileName.c_str(),"w"))
  {
#ifdef DEBUG
    cout << "BwtWriterBase opened file " << fileName << " " << pFile_ << endl;
#endif
   readWriteCheck(fileName.c_str(),1);    //    setvbuf( pFile_, NULL, _IOFBF, 262144);
  }
BwtWriterBase::~BwtWriterBase() 
{ 

fclose(pFile_); 
#ifdef DEBUG
 cout << "BwtWriterBase: closed file " <<pFile_ << endl;
#endif

}


//
// BwtWriterASCII member function definitions
//

BwtWriterASCII::BwtWriterASCII( const std::string& fileName ) : BwtWriterBase(fileName)
  {
#ifdef DEBUG
  cout << "BW ASCII ctor" << endl;
#endif
  }
BwtWriterASCII::~BwtWriterASCII() 
{
#ifdef DEBUG
  cout << "BW ASCII dtor" << endl;
#endif

}


void BwtWriterASCII::operator()( const char* p, int numChars )
{
#ifdef DEBUG
  cout << "BW ASCII () - " << *p << " " << numChars << endl;
#endif


  size_t bytesWritten = fwrite( p, sizeof(char), numChars, pFile_ );
  if (bytesWritten != (size_t)numChars) {
    cerr << "Unable to write "<< numChars
	 << " chars. Aborting." << endl;
    exit(-1);
  }
} // ~operator()

void BwtWriterASCII::sendRun( char c, int runLength )
{
  for (int i(0);i<runLength;i++) fprintf( pFile_, "%c", c );
}

//
// BwtWriterRunLength member function definitions
//

BwtWriterRunLength::~BwtWriterRunLength()
  { 
    if (runLength_!=0) encodeRun( lastChar_, runLength_); 


    if (pBuf_!=buf_)
    {

       size_t bytesWritten = fwrite( buf_, sizeof(char), (pBuf_-buf_), pFile_ );
       if (bytesWritten != (size_t)(pBuf_-buf_)) {
           cerr << "Unable to write "<< (pBuf_-buf_)
                   << " chars. Aborting." << endl;
                   exit(-1);
       }
      
#ifdef REPORT_COMPRESSION_RATIO 
      bytesWritten_+=(LetterCountType)(pBuf_-buf_);
#endif
    } 

#ifdef REPORT_COMPRESSION_RATIO 
    std::cout << "BwtWriterRunLength: received " 
	      << charsReceived_ << " chars, sent "
	      << bytesWritten_ << " bytes, compression " 
	      << ((double)8*bytesWritten_)/(charsReceived_) 
	      << " bits per char " << std::endl;
#endif


  }

void BwtWriterRunLength::operator()( const char* p, int numChars )
{
#ifdef DEBUG
  std::cout << "BW RL () - " << *p << " " << numChars << " state: " << lastChar_ << " " << runLength_ << endl;
#endif

    for (int i(0);i<numChars;i++)
    {
      if ((*p)==lastChar_) { runLength_++; }
      else
      {
	if (runLength_>0)
	{
	  encodeRun( lastChar_, runLength_ );
	} // ~if
	runLength_=1; lastChar_=*p;
      } // ~else 
      p++;
    } // ~for
    //    assert(fwrite( p, sizeof(char), numChars, pFile_ )==numChars);
  } // ~operator()

  void BwtWriterRunLength::sendChar( char c)
  {
      *pBuf_=c;
      if (++pBuf_==pBufMax_)
      {
#ifdef REPORT_COMPRESSION_RATIO 
	bytesWritten_+=ReadBufferSize;
#endif

       size_t bytesWritten = fwrite( buf_, sizeof(char), ReadBufferSize, pFile_ );
       if (bytesWritten != (size_t)ReadBufferSize) {
           cerr << "Unable to write "<< ReadBufferSize
                   << " chars. Aborting." << endl;
                   exit(-1);
       }
	pBuf_=buf_;
      }
  }

  void BwtWriterRunLength::encodeRun( char c, uint runLength )
  {
#ifdef DEBUG
    std::cout << "BW RL encodeRun - sending run " << c << " " << runLength << " " << pFile_ 
	      << std::endl;
#endif
#ifdef REPORT_COMPRESSION_RATIO 
    charsReceived_+=runLength;
#endif    
int charIndex(whichPile[(int)c]);

    if(charIndex==nv){
        cerr << "Char is not part of the alphabet. Aborting." << endl;
        exit(-1);
    }
    
    uchar outCode(0xF0|((uchar)charIndex));
    runLength--;
    const uint numMaxChars(runLength>>4);
    for (uint i(0);i<numMaxChars;i++)
    {
      sendChar(outCode);
    }
    runLength&=(uint)0xF;

    outCode=(((uchar)runLength)<<4);
    outCode|=charIndex;
    //    assert(((uint)outCode)<256);
    //    assert(fwrite( &outCode, sizeof(char), 1, pFile_ )==1);
    sendChar(outCode);
#ifdef DEBUG
    std::cout << "B sending " << (uint)outCode << " " << pFile_ << std::endl;
#endif
  } // ~encodeRun

  void BwtWriterRunLength::sendRun( char c, int runLength )
  {
#ifdef DEBUG
    std::cout << "BW RL sendRun - sending run " << c << " " << runLength << " " << endl; 
#endif
    if (runLength!=0)
    {
      if (c==lastChar_) { runLength_+=runLength; }
      else
      {
	if (runLength_!=0) encodeRun(lastChar_,runLength_);
	lastChar_=c;
	runLength_=runLength;
      }
    }
  } // ~sendRun



#ifdef ORIG
  void BwtWriterRunLength::sendRun( char c, int runLength )
  {
#ifdef DEBUG
    std::cout << "BW RL sendRun - sending run " << c << " " << runLength << " " << endl; 
#endif
    if (c==lastChar_) { runLength_+=runLength; }
    else
    {
      if (runLength_!=0) {encodeRun(lastChar_,runLength_);
      lastChar_=c;
      runLength_=runLength;}
    }
  } // ~sendRun
#endif


 // Huffman implementation
  
  BwtWriterHuffman::~BwtWriterHuffman() // destructor
   {
    sendRun(lastChar_,runLength_); //gets last normal chars from buffer
    sendRun(notInAlphabet,1);   // send termination char
    BwtWriterHuffman::emptyBuffer();
    while (bitsUsed_>0)
    {
      assert(fwrite(&soFar_.ui, sizeof(unsigned int), 1, pFile_)==1);
#ifdef REPORT_COMPRESSION_RATIO 
      bytesWritten_+=(LetterCountType) sizeof(unsigned int);
#endif
#ifdef DEBUG
      cout << endl;
      for (unsigned int i(1);i!=0;i<<=1)
	cout << ((soFar_.ui&i)?'1':'0');
      cout << endl;
#endif
      bitsUsed_-=32;
    }
#ifdef REPORT_COMPRESSION_RATIO 
    std::cout << "BwtWriterHuffman: received " 
	      << charsReceived_ << " chars, sent "
	      << bytesWritten_ << " bytes, compression " 
	      << ((double)8*bytesWritten_)/(charsReceived_) 
	      << " bits per char " << std::endl;
#endif
  } // ~BwtWriterHuffman() 

 void BwtWriterHuffman::operator()( const char* p, int numChars )
  { 
    for (int i(0);i<numChars;i++)
    {
      if ((*p)==lastChar_) { runLength_++; }
      else
      {
	if (runLength_>0)
	{
            sendRun( lastChar_, runLength_ );
	} // ~if

	runLength_=1; lastChar_=*p;
      } // ~else 
      p++;
    } // ~for
        sendRun(lastChar_, runLength_);
        runLength_ = 0;
  } // ~operator()

    void BwtWriterHuffman::sendToken(unsigned long long code, uint length)
    {
        toAdd_.ull = code;

        toAdd_.ull <<= bitsUsed_; // left shift to the next free position
        soFar_.ull |= toAdd_.ull; // update so far
        bitsUsed_ += length;

        if (bitsUsed_ > 32) // if we have more than 32bit / 4 byte
        {
            assert(fwrite(&soFar_.ui, sizeof (unsigned int), 1, pFile_) == 1);
#ifdef REPORT_COMPRESSION_RATIO 
      bytesWritten_+=(LetterCountType) sizeof(unsigned int);
#endif
#ifdef DEBUG
            for (unsigned int i(1); i != 0; i <<= 1)
                cout << ((soFar_.ui & i) ? '1' : '0');
            cout << endl;
#endif
            soFar_.ull >>= 32; // shift rest to the right border
            bitsUsed_ -= 32; // update bits used
        } // 

    } // ~sendToken()

    void BwtWriterHuffman::sendRun(char c, int runLength)
    {
#ifdef REPORT_COMPRESSION_RATIO 
    charsReceived_+=runLength;
#endif    
        if (runLength > 0) {

            for (int i(0); i < runLength; i++) {

                if (huffmanBufferPos == huffmanWriterBufferSize-1) {
                    symBuf[huffmanBufferPos] = c;                   
                    processBuffer(huffmanWriterBufferSize);
                } else {
                    symBuf[huffmanBufferPos] = c;
                    huffmanBufferPos++;
                } // ~else 

            } // ~for
        }
    } // ~sendRun

    void BwtWriterHuffman::emptyBuffer(void)
    {
        processBuffer(huffmanBufferPos);
    }
    
    void BwtWriterHuffman::processBuffer(int itemsToPrint)
    {
      //  cerr << pFile_ << "PROCESSING BUFFER" << endl;   

        if (itemsToPrint > 0){
        char localLastChar = 0;
        uint localRunLength=0;

        for (int i(0);i<itemsToPrint;i++)
            {
           //  cerr << pFile_ << " accessing char at " << i << endl;      
              if (symBuf[i]==localLastChar) { localRunLength++; }
              else
              {
                if (localRunLength>0)
            {
                    // get number of this char, 0-5
                    int charIndex(whichPile[(int) localLastChar]); 
                    assert(charIndex != nv); // crash if not from alphabet

                    if (localRunLength == 1) // single run only
                    {
                        sendToken(singleCharCode[charIndex],
                                  singleCharLength[charIndex]);
                    }// ~if
                            else 
                            {
                                sendToken(doubleCharCode[charIndex],
                                          doubleCharLength[charIndex]);
                                sendNum(localRunLength);
                            } //~else
                        } // ~if

                        localRunLength = 1;
                        localLastChar = symBuf[i];
                    } // ~else 
                } // ~for
        
       // process last entry of the buffer
         int charIndex(whichPile[(int) localLastChar]); // get number of this char, 0-5
         
         if ((int)localLastChar > 0 && whichPile[(int) localLastChar]< alphabetSize ){
         
         assert(charIndex != nv); // crash if not from alphabet
         
         if (localRunLength == 1) // single run only
        {
            sendToken(singleCharCode[charIndex], singleCharLength[charIndex]);
        }// ~if
                else {
                    sendToken(doubleCharCode[charIndex], doubleCharLength[charIndex]);
                    sendNum(localRunLength);
                } //~else       
            }
            huffmanBufferPos = 0; // reset counter
        }
    } // ~sendRun

    void BwtWriterHuffman::sendNum(uint runLength)
    {
        if (runLength < 17) // max 16 
        {
            runLength--;
            runLength--;
            numBuf_.ui = runLength; // set  new run length
            sendToken(numBuf_.ull, 4); // write one token, encoding for the runlength
        }// ~if
        else // larger than 16 -> 2 byte
        {
            runLength -= 17; // substract 16 + 1
            numBuf_.ui = 0xF; // set uint to 16 -> 1111
            sendToken(numBuf_.ull, 4); // send binary encoded 16 using 4 bits
            do {
                numBuf_.ui = runLength; // set uint to remaining runlength
                numBuf_.ui &= 0x7F; // AND with 111|1111
                // 1st: OR with 1000|0000 2nd multiply with 1 128 if RL > 127 or with 0
                numBuf_.ui |= 0x80 * (runLength > 0x7F);
                sendToken(numBuf_.ull, 8);
                runLength >>= 7;
            }// ~while
            while (runLength != 0);
        } // ~else      
    } // sendNum

//
// BwtWriterRunLength member function definitions
//

BwtWriterImplicit::~BwtWriterImplicit() 
{ 
  if (inSAP_==true)
  {
    flushSAP();
  }
  else if (lastChar_!=notInAlphabet)
  {
    pWriter_->sendRun(lastChar_, lastRun_);
  }
  delete pWriter_; 
}

void BwtWriterImplicit::flushSAP( void )
{
  assert(alphabet[firstSAP_]==lastChar_);
  if (countSAP_.count_[firstSAP_]>0) pWriter_->sendRun(alphabet[firstSAP_],countSAP_.count_[firstSAP_]);

  for (int i(0);i<alphabetSize;i++) 
  {
    if ((i!=firstSAP_)&&(countSAP_.count_[i]>0)) pWriter_->sendRun(alphabet[i],countSAP_.count_[i]);
  }
}

  
void BwtWriterImplicit::operator()( const char* p, int numChars )
{

  for (int i(0);i<numChars;i++,p++)
  {
    if (islower(*p))
    {
      if (inSAP_==false)
      {
	countSAP_.clear();
	assert (lastChar_!=notInAlphabet);
	firstSAP_=whichPile[(int)lastChar_];
	assert(firstSAP_!=nv);
	countSAP_.count_[firstSAP_]+=lastRun_;
	inSAP_=true;
      } // ~if
      countSAP_+=*p;
    } // ~if
    else
    {
      if (inSAP_==true)
      {
	flushSAP();
	inSAP_=false;
      } 
      else if (lastChar_!=notInAlphabet)
      {
	pWriter_->sendRun(lastChar_, lastRun_);
      }
      lastChar_=*p;
      lastRun_=1;
    }
  }
}

void BwtWriterImplicit::sendRun( char c, int runLength )
{
  if (islower(c))
  {
    if (inSAP_==false)
    {
      countSAP_.clear();
      assert (lastChar_!=notInAlphabet);
      firstSAP_=whichPile[(int)lastChar_];
      assert(firstSAP_!=nv);
      countSAP_.count_[firstSAP_]+=lastRun_;
      inSAP_=true;
    } // ~if
    countSAP_.count_[whichPile[(int)c]]+=runLength;
  }
  else
  {
    if (inSAP_==true)
      {
	flushSAP();
	inSAP_=false;
      } 
      else if (lastChar_!=notInAlphabet)
      {
	pWriter_->sendRun(lastChar_, lastRun_);
      }
      lastChar_=c;
      lastRun_=runLength;
  }

  //  (*pWriter_).sendRun(toupper(c), runLength);
}

#ifdef XXX  
void BwtWriterImplicit::operator()( const char* p, int numChars )
{
  // could be smarter about this
  char c;
  for (int i(0);i<numChars;i++,p++)
  {
    c=toupper(*p);
    (*pWriter_)(&c, 1);
  }
}

void BwtWriterImplicit::sendRun( char c, int runLength )
{
  (*pWriter_).sendRun(toupper(c), runLength);
}
#endif
