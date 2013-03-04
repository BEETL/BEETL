// tutorial about suffix arrays.
#include <iostream>
//#include <seqan/index.h>
#include <cassert>
#include <vector>
#include <queue>

#include "Alphabet.hh"
#include "Timer.hh"
#include <fstream>

//using namespace seqan;
using namespace std;
//#define DEBUG;
 	
void getFileName( const string& stem, const char code, const int pile,
		  string& fileName )
{
  fileName=stem;
  fileName+='-';
  fileName+=code;
  fileName+='0';
  assert(pile<=9);
  fileName+=(char)(48+pile);
  // cerr << "Made file name " << fileName << endl;
}


long getFileSize(FILE *file)
{
  long lCurPos, lEndPos;
  lCurPos = ftell(file);
  fseek(file, 0, 2);
  lEndPos = ftell(file);
  fseek(file, lCurPos, 0);
  return lEndPos;
}

struct SuffixInfo
{
  // FILE* pFile_;
  unsigned short fileNum_;
  unsigned char bwtSymbol_;
  const char* pSeq_;
  unsigned loc_;
};
 	
bool operator<(const SuffixInfo& lhs, const SuffixInfo& rhs)
{
  // assert(1==0);
  int compare(strcmp(lhs.pSeq_+lhs.loc_,rhs.pSeq_+rhs.loc_));

#ifdef DEBUG

    cout << "compare l " << lhs.fileNum_ <<" " << lhs.loc_ <<  "  " << strlen(lhs.pSeq_)<< endl;
    cout << "compare r " << rhs.fileNum_ <<" " << rhs.loc_ << "  " << strlen(rhs.pSeq_) << endl;
    cout << compare <<endl; 
    int i(0);
    while(lhs.pSeq_[lhs.loc_+i] == rhs.pSeq_[rhs.loc_+i]) i++;
    cout << i << " " << lhs.pSeq_[lhs.loc_+i] << " " << rhs.pSeq_[rhs.loc_+i]<<endl;
#endif DEBUG

  // if (compare==0) cout << "Hello" << endl;
  //       cout << lhs.fileNum_ << " " << lhs.loc_ << " " << rhs.fileNum_ << " " << rhs.loc_ << " " << lhs.pSeq_+lhs.loc_ << " " << (long)rhs.pSeq_ << endl;
  return ((compare==0)?(lhs.fileNum_>rhs.fileNum_):(compare>0));
  // return ((compare==0)?(lhs.fileNum_>rhs.fileNum_):(compare>0));
  
  // return (strcmp(lhs.pSeq_+lhs.loc_,rhs.pSeq_+rhs.loc_)>0);
} 	
 	
 	
void readFastaFile( const char* name, vector<char>& data ){
  ifstream fastaFile(name,ios::in);
  string line;
  while(fastaFile.good()){
    getline(fastaFile,line);
    if(line[0] != '>'){
      for(int i(0) ; i<line.length() ; i++){
	if( strchr("ACGNT",line[i])==NULL){
	  cout << "found wrong char " << line[i] << " in " << name << endl;
	  line[i] = 'N';
	}
	data.push_back(line[i]);
      }
    }
  }
  fastaFile.close();
  data.push_back('\0');
} // ~readFastaFile


int main ( int numArgs, const char* args[] )
{
  if (numArgs<4)
    {
      cerr << "Usage: " << args[0] << " pileNum outputPrefix files" << endl;
      exit(EXIT_FAILURE);
    }
  
  int pileNum(atoi(args[1]));
  cerr << "Will merge BWTs from pile " << pileNum << endl;  
  //assert((pileNum>=0)&&(pileNum<alphabetSize));
  
  cerr << "Output file will be named " << args[2] << endl;
  
  vector<vector<char> > chromSeqs;

  vector<queue<unsigned> > suffixArrays;
  //  vector<vector<unsigned>::iterator> suffixArraysEnd;
  vector<queue<char> > bwts;
  //  vector<string::iterator> bwtsEnd;
  
  //  vector<vector<unsigned> > suffixstore;
  //vector<string> bwtstore;
  
  FILE* arrayInputStream;
  FILE* bwtInputStream;
  FILE* pFileChrom;
  
  FILE* pFileArray;
  FILE* pFileBWT;
 
  string fileNameStem(args[2]), fileName;
 	
  getFileName(fileNameStem, 'A', pileNum, fileName);
  cerr << "Will send suffix array entries to " << fileName << endl;
  pFileArray=fopen(fileName.c_str(),"w");
  assert(pFileArray!=NULL);
  
  getFileName(fileNameStem, 'B', pileNum, fileName);
  cerr << "Will send BWT entries to " << fileName << endl;
  pFileBWT=fopen(fileName.c_str(),"w");
  assert(pFileBWT!=NULL);
  
  getFileName(fileNameStem, 'C', pileNum, fileName);
  cerr << "Will send genome information to " << fileName << endl;
  pFileChrom=fopen(fileName.c_str(),"w");
  assert(pFileChrom!=NULL);
 	

  SuffixInfo thisSuffix;
  priority_queue<SuffixInfo,vector<SuffixInfo> > pq;
  ofstream fileCounter("filecounter.csv", ios::out);
  int fileN(0);
  for (int i(3); i<numArgs;i++)
    {
      //      cout << i <<endl;
      //read the fasta file in
      chromSeqs.push_back(vector<char>());
      readFastaFile(args[i],chromSeqs.back());
     
      //cerr << "Sequence is of length " << chromSeqs.back().size() << endl;
      if(i %20 ==0) {
	cerr << "Sequence " << i <<  " has length " << chromSeqs.back().size() << endl;
      }
      fileNameStem="bwt_"+(string)args[i];
      getFileName(fileNameStem, 'A', pileNum, fileName);
      //      cerr << "Will open suffix array file " << fileName << endl;
      //read in the whole suffix array
      queue<unsigned> thisSuffArr;
      queue<char> bwtQue;
      suffixArrays.push_back(thisSuffArr);
      bwts.push_back(bwtQue);
      arrayInputStream = fopen(fileName.c_str(),"r");      
      if(arrayInputStream != NULL){
	long fileSize = getFileSize(arrayInputStream);
 	//allocate the space
	unsigned* arrayBuf = new unsigned[fileSize/4];
	assert(fread(arrayBuf,4,fileSize/4,arrayInputStream) == (fileSize/4));
	fileCounter << fileN << "," << args[i] <<endl;
	fileN++;
	//      cout << "read suffix\n";
	for(int k(0) ; k< (fileSize / 4); k++){
	  thisSuffArr.push(arrayBuf[k]);
	} 
	suffixArrays.back() = thisSuffArr;
	fclose(arrayInputStream);
      }
      if(!suffixArrays.back().empty()){
	//read the bwt
	getFileName(fileNameStem, 'B', pileNum, fileName);
	if(i%20 ==0){cerr << "Will open BWT file " << fileName << endl;}
	ifstream bwtIn(fileName.c_str(),ios::in);
       if(bwtIn.good()){
	  string line;
	  getline(bwtIn, line);
	  for(int s(0); s< line.length(); s ++)
	    bwtQue.push(line[s]);
	  bwts.back() = bwtQue;
	}
      }
      //      cerr <"got bwt" <<endl;
      else
	{ cerr << "Warning: could not open " << fileName << endl; }
    }
  //chromSeqs = whole sequence out of fasta file
  // for each bacterium initialise one suffix
  cerr << "Got " << chromSeqs.size() << " sequences "<< endl;
  cerr << "Got " << suffixArrays.size() << " suffixe" << endl;
  fileCounter.close();
  for (unsigned int i(0);i<chromSeqs.size();i++)
    { 
      if(!suffixArrays[i].empty()){
	//pointer to first char in vector<char> of the genome sequence
	thisSuffix.pSeq_=&(chromSeqs[i][0]);
	thisSuffix.fileNum_= (unsigned short) i;
	thisSuffix.loc_ = suffixArrays[i].front();
	suffixArrays[i].pop();
	thisSuffix.bwtSymbol_ = bwts[i].front();
	bwts[i].pop();
	//      bwtIndices[i] += 1;
	pq.push(thisSuffix);
      }
    }
 // ~for
  
  int j(0);
  
#ifdef DEBUG
  const char* prevSuff(NULL);
  const char* thisSuff(NULL);
  unsigned char prevNum;
#endif
  
  while(!pq.empty())
    {
      thisSuffix=pq.top();
       
#ifdef DEBUG
      for(int jj(0); jj<50 ;jj++)
	cout <<*(thisSuffix.pSeq_+thisSuffix.loc_+jj);
      cout << endl;

      thisSuff=thisSuffix.pSeq_+thisSuffix.loc_;
      if (prevSuff!=NULL)
 	{
	  assert((strcmp(prevSuff,thisSuff)<1)||
		 ((strcmp(prevSuff, thisSuff)==0)&&(prevNum<thisSuffix.fileNum_)));
 	}
      prevSuff=thisSuff; prevNum=thisSuffix.fileNum_;
#endif

      pq.pop();
      assert(fwrite(&thisSuffix.loc_,sizeof(unsigned),1,pFileArray)==1);
      assert(fwrite(&thisSuffix.bwtSymbol_,sizeof(unsigned char),1,pFileBWT)==1);
      assert(fwrite(&thisSuffix.fileNum_,sizeof(unsigned short),1,pFileChrom)==1);
      
      if ((j%1000000)==0) cout << "." << endl;
      j++;
#ifdef DEBUG
      //   for (int j(0);j<50;j++)
      //    	cout << *(thisSuffix.pSeq_+thisSuffix.loc_+j);
      // cout << endl;
#endif
     //read 4 bytes
      
      if(!suffixArrays[thisSuffix.fileNum_].empty()){
	thisSuffix.loc_= suffixArrays[thisSuffix.fileNum_].front();
	suffixArrays[thisSuffix.fileNum_].pop();
	//	cout << "suffix loc assigned " << suffixArrays[thisSuffix.fileNum_][suffixIndices[thisSuffix.fileNum_]];
	//	suffixIndices[thisSuffix.fileNum_] +=1;
	thisSuffix.bwtSymbol_ = bwts[thisSuffix.fileNum_].front(); 
	bwts[thisSuffix.fileNum_].pop();
	//	cout << "bwt assigned " << bwts[thisSuffix.fileNum_][bwtIndices[thisSuffix.fileNum_]];
	//	bwtIndices[thisSuffix.fileNum_] +=1;
       //open bwt

      //	  bwtInputStream = fopen(filesBWT[thisSuffix.fileNum_].c_str(),"r");
	  //set right position
      //  fseek(bwtInputStream, bwtStreamPositions[thisSuffix.fileNum_], SEEK_SET);
	  //read one byet
      //  assert(fread(&thisSuffix.bwtSymbol_,sizeof(unsigned char),1,
      //	       bwtInputStream)==1);
	  //remember new position
      //  bwtStreamPositions[thisSuffix.fileNum_] += sizeof(unsigned char);
	  //close bwt
      //	  fclose(bwtInputStream);
	pq.push(thisSuffix);
      }
    } // ~while
  
  
  /*  for (unsigned int i(0);i<files.size();i++)
    {
      if (files[i]!=NULL)
 	{
	  fclose(files[i]);
 	}
      
      if (filesBWT[i]!=NULL)
 	{
	  fclose(filesBWT[i]);
 	}
    }
  */
  fclose(pFileArray);
  fclose(pFileBWT);
  fclose(pFileChrom);
 	
} // ~main
