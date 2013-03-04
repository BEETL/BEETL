#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<cassert>
#include<sstream>
#include<Alphabet.hh>

using namespace std;


struct HitInformation
{
  unsigned bwtPosition_;
  unsigned short charCount_[6];
  string suffix_;
};

struct SequenceInformation
{
  unsigned short fileNum_;
  string tag_;
  string fileName_;
  unsigned int sequenceLength_;  
  vector<HitInformation> hits_;
  vector<unsigned> suffixStarts_;
  vector<unsigned> suffixEnds_;
  vector<unsigned int> suffixCoverage_;
};



vector<SequenceInformation> loadSequenceInfoToTest(vector<string> sequenceFiles, vector<string> mergingOut);

vector<SequenceInformation> loadSequenceInfoToParse(string sigmaSuffixOut, string headerFile, string fileCounter, string mergCZeroOut); 

unsigned short findFileNum(long BWTPosition, FILE* mergeCPileOut);

unsigned int findSuffixPosition(long BWTPosition, string mergAPileOut);

void testSigma(vector<SequenceInformation> &info, FILE* mergZeroA, FILE* mergSeroC);

void printUsage(string message);

vector<string> split(string s, string token);

void parseWordCountInformation(string wordCountOut, vector<SequenceInformation> &info, 
vector< vector <unsigned short> > fileCount , 
int minKMer); 

void getCoverageInformation(vector<SequenceInformation> &info, vector<string> mergingSuffixOut); 

void testSuffixe(vector<SequenceInformation> &info, vector<FILE*> suffixPositions, vector<FILE*> mergCOutput);

void testCountWordOutput(vector<SequenceInformation> &seqs, vector<FILE*> mergeCFiles, vector<string> mergeSuffixFiles);

vector < vector<unsigned short> > loadMergCInformation(vector<string> mergeCFiles);

//int whichPile(char c);

long getFileSize(FILE *file)
{
  long lCurPos, lEndPos;
  lCurPos = ftell(file);
  fseek(file, 0, SEEK_END);
  lEndPos = ftell(file);
  fseek(file, lCurPos, SEEK_SET);
  return lEndPos;
}

int main(int argc, char **argv){
  if(argc < 4){
    printUsage("Not enoug arguments given");
    return 1;
  }
  //     cout << "Please give a fileName of the countWords and what you want to do with it \n";
  string countWordOut;
  vector<string> sequenceFileNames;
  vector<string> mergingCOutput;
  vector<string> suffixFiles;
  string sigmaSuffixOut;
  string headerFile;
  string fileCounter;
  string suffixZeroFile;
  int minKMer;
  
  bool countFound(0), mCFound(0), mIFound(0), mAFound(0), testMerging(0), parseCount(0), parseCountForCoverage(0), mA0Found(0), fHFound(0); 
  
  //  cerr<< "Usage: met -t -mC <Files> -mI <Files> -mC <Files> -mA <Files> "<<endl
  //  << "or: met -c -cO File -mC <Files> -mF File, -mA0 File, -fH"<<endl
  //  << "or: met -c -t -cO File -mC <Files> -mF <File> -mA <Files>" <<endl<< endl;
  
  for(int i(1); i <argc;i++){
    if(argv[i][1] == 't'){
      testMerging=true;
    }
    else if(argv[i][1] == 'p'){
      parseCount =true;
    }
    else if(argv[i][1] == 'c'){
      parseCount =true;
      parseCountForCoverage =true;
    }
    if(strcmp(argv[i],"-o") == 0){
      countWordOut = argv[i+1];
      cout <<endl<<"count " << argv[i+1] <<endl;
      countFound =true;
    }
    else if(strcmp(argv[i], "-mF") ==0){
      fileCounter = argv[i+1];
    }
    else if(strcmp(argv[i], "-mA0")==0){
      mA0Found = true;
      suffixZeroFile = argv[i+1];
    }
    else if(strcmp(argv[i] , "-fH") == 0){
      fHFound =true;
      headerFile = argv[i+1];
    }
    else if(strcmp(argv[i],"-mI") == 0){
      mIFound =true;
      for(int j(i+1); j<argc ; j++){
	if(argv[j][0] == '-')
	  break;
	//	cout << "seq " << argv[j] << endl;
	sequenceFileNames.push_back(argv[j]);
      }
      cout <<"Got " << sequenceFileNames.size() << " sequenceFiles." <<endl;
    }
    else if(strcmp(argv[i],"-mC") == 0){
      mCFound =true;
      for(int j(i+1); j<argc; j++){
	if(argv[j][0]== '-')
	  break;
	mergingCOutput.push_back(argv[j]);
	cout << "mergC " << argv[j] <<endl;
      }
    }
    else if(strcmp(argv[i], "-mA") == 0){
      for(int j(i+1); j<argc; j++){
        if(argv[j][0] =='-')
          break;
        suffixFiles.push_back(argv[j]);
	mAFound =true;
	mA0Found =true;
      }
    }
    else if(strcmp(argv[i], "-k") ==0)
      minKMer = atoi(argv[i+1]);
  }
  
  if(testMerging){ 
    if(!mIFound){
      printUsage("no Sequence files given.");
      return 1;
    }
    if(!mCFound){
      printUsage("no merging C files are given");
      return 1;
    }
    if(!mAFound){
      printUsage("no Suffix files given" );
      return 1;
    }
  }
  if(parseCount){
    if(!fHFound ){
      printUsage("no HeaderFile is given");
      return 1;
    }
    if(!mA0Found){
      printUsage("no A0 file found");
      return 1;
    }
    if(!countFound) {
      printUsage("No countFile found");
      return 1;
    }
    if( !mCFound){
      printUsage("No mc found");
      return 1;
    }
  }
  if(parseCountForCoverage){
    if(!mAFound){
      printUsage("No Suffix files are given");
      return 1;
    }
    suffixZeroFile = suffixFiles[0];
  }

  

  
  vector<SequenceInformation> seq;
  if(testMerging) 
    seq= loadSequenceInfoToTest(sequenceFileNames,mergingCOutput);
  else
    seq = loadSequenceInfoToParse(suffixZeroFile,headerFile, fileCounter, mergingCOutput[0] );
  
  

  vector<FILE*> mergs;
  mergs.resize(mergingCOutput.size());
  for(unsigned int i(0); i< mergingCOutput.size(); i++){
    stringstream s;
    s<<"C0"<<i;
    if (mergingCOutput[i].find(s.str())!= string::npos){
      mergs[i] = fopen(mergingCOutput[i].c_str(),"r");
      //  cout <<i <<" " << s.str() << " " << mergingCOutput[i] <<endl;
    }
  }
  /*
  //test postions
  //A 522428 -> G_1441_rev "works"
  unsigned short fileNum = findFileNum(522428,mergs[1]);
  cout << fileNum << "  " << seq[fileNum].fileName_ <<endl;

  //T 522428 -> G_2051 "works" 
  fileNum = findFileNum(522428,mergs[5]);
  cout << "T " << fileNum << "  " << seq[fileNum].fileName_ << "\t should be G_2051" << endl;
  //G 522428   fileNum =G_495_rev "works" 
  fileNum = findFileNum(522428,mergs[3]);
  cout << "G " << fileNum << "  " << seq[fileNum].fileName_ << "\t should be G_495_rev" << endl;
  //C 401233 -> G_1166 "works" 
  fileNum = findFileNum(401233,mergs[2]);
  cout << "C " << fileNum << " " << seq[fileNum].fileName_ << "\t should be G_1166"<<endl;
  //C  623101 -> G_494 does NOT work
  fileNum = findFileNum(623101,mergs[2]);
  cout << "C " << fileNum << " " << seq[fileNum].fileName_ << "\t should be G_494"<<endl;
  //T 1069383 -> G_871 "does NOT  work"
  fileNum = findFileNum(1069383,mergs[5]);
  cout << "T " << fileNum << " " << seq[fileNum].fileName_ << "\t should be G_871 " <<endl;
  
  //to test one of error T 1069384 
  fileNum = findFileNum(1069384,mergs[5]);
  cout << "T " << fileNum << " " << seq[fileNum].fileName_ << "\t should be G_871 " <<endl;
  */

  //to test one of error T 1069382
  fileNum = findFileNum(1069382,mergs[5]);
 cout << "T " << fileNum << " " << seq[fileNum].fileName_ << "\t should be G_871 " <<endl;

  //A 1285939 -> G_1522_rev "works" 
  fileNum = findFileNum(1285939,mergs[1]);
  cout << "A " << fileNum << " " << seq[fileNum].fileName_ << "\t shoudl be G_1522_rev " << endl;

  //cout << "SuffixFiles " <<suffixFiles.size() << <<endl;
  if(testMerging){
    vector<FILE*> suffPositions;
    for(unsigned int i(0) ; i< suffixFiles.size() ; i++)
      suffPositions.push_back(fopen(suffixFiles[i].c_str(),"r"));
    if(suffPositions[0] == NULL)
      printUsage("Output A00 of merging is not readable");
    testSigma(seq, suffPositions[0], mergs[0]);
    
     testSuffixe(seq, suffPositions, mergs);
     cout << " start parsing " << endl;
     vector < vector< unsigned short > > fileNumbers = loadMergCInformation(mergingCOutput);

     parseWordCountInformation(countWordOut,seq,fileNumbers,minKMer);
     testCountWordOutput(seq, mergs, suffixFiles);
  }
  cout << "start parsing " <<endl;
  if(parseCount){
    //getCoverageinfo
    vector < vector< unsigned short > > fileNumbers = loadMergCInformation(mergingCOutput);
    parseWordCountInformation(countWordOut, seq, fileNumbers, minKMer);
    for(unsigned int i(0) ; i <seq.size(); i++){
      if(seq[i].hits_.size() > 0)
	cout << seq[i].fileName_ << " " << seq[i].tag_ << seq[i].hits_.size() <<endl;
    }
  }
  if(parseCountForCoverage){
    getCoverageInformation(seq, suffixFiles);
  }
}

vector <vector<unsigned short> > loadMergCInformation(vector<string> mergeCFiles) {

  vector< vector<unsigned short> > fileCountInfo;
  fileCountInfo.resize(alphabetSize);
  for(unsigned int i (0) ; i<mergeCFiles.size(); i++){
    cout << i << " open File " << mergeCFiles[i] << endl;
    FILE* merC = fopen(mergeCFiles[i].c_str(), "r" );
    unsigned short count;
    while(fread(&count, sizeof(unsigned short ), 1,merC) == 1)
      fileCountInfo[i].push_back(count);
    fclose(merC);
    cout <<i <<" has " << fileCountInfo[i].size() <<" bwtPositions " <<endl;
  }
  return fileCountInfo;
}



void getCoverageInformation(vector<SequenceInformation> &seqInfo, vector<string> mergAOutput){
  cout << "getCoverageInformation" << endl;
  //go through all the sequenceinformation
  for(unsigned int i(0) ; i<seqInfo.size(); i++){
    //    if(seqInfo[i].hits_.size() != 0 ){
      //      vector<unsigned> suffixStarts;
      // vector<unsigned> suffixEnds;
      //vector<unsigned int> suffixCoverage;
      //cout << "Seq got hits "<< seqInfo[i].fileName_ << "\t" << seqInfo[i].hits_.size() <<endl;
      //cout << "Rev Seq got hits "<< seqRev.fileName_ <<"\t" <<seqRev.hits_.size() <<endl;
      for(unsigned int k(0) ; k<seqInfo[i].hits_.size() ; k++){
	unsigned pos;
	int pile = whichPile[(int)seqInfo[i].hits_[k].suffix_[0]];	  
	long bwtPosition = seqInfo[i].hits_[k].bwtPosition_;
	//	cout << "bwtPos " << bwtPosition <<endl;
	pos = findSuffixPosition(bwtPosition, mergAOutput[pile]);
	//	cout << seqInfo[i].fileName_ <<endl;
	//	cout <<seqInfo[i].hits_[k].suffix_;
	//	cout << "S " << pos <<endl;
	//cout << "L " << seqInfo[i].sequenceLength_<<endl;
      }
      /*if(seqInfo[i].fileName_.find("rev") == string::npos){
	seqInfo[i].suffixStarts_.push_back(pos);
	  seqInfo[i].suffixEnds_.push_back(pos +seqInfo[i].hits_[k].suffix_.length());
	  }
	  else{
	  seqInfo[i].suffixStarts_.push_back(seqInfo[i].sequenceLength_ - (pos + seqInfo[i].hits_[k].suffix_.length()));
	  seqInfo[i].suffixEnds_.push_back(seqInfo[i].sequenceLength_ - pos);
	}
	int coverage(0);
	for(unsigned int c(0) ; c< 6 ; c++){
	coverage = coverage + seqInfo[i].hits_[k].charCount_[c];
	}
	seqInfo[i].suffixCoverage_.push_back(coverage);
	}*/
	//      seqInfo[i].suffixStarts_ = suffixStarts;
      //seqInfo[i].suffixEnds_ = suffixEnds;
      //seqInfo[i].suffixCoverage_ = suffixCoverage;
      // }
  }
}


void parseWordCountInformation(string countWordOut, vector<SequenceInformation> &seqInfo, 
			       vector< vector<unsigned short> > filePositions,
			       int minK)
{
  ifstream countWordsIn(countWordOut.c_str(), ios::in);
  string line;
  int i(0);
  //  HitInformation hit;
  while(countWordsIn.good()){
    i++;
    vector<string> lineVector;
    getline(countWordsIn,line);
    lineVector = split(line," " );
    if(lineVector[0].compare("BKPT") == 0){
      if(lineVector[1].length() > (unsigned) minK){
	//hit.suffix_ = lineVector[1];
	//hit.bwtPosition_ = atoi(lineVector[2].c_str());	
	//    cout << lineVector[1].lengthendl;                                                                                                       
	vector<string> positionsRead = split(lineVector[3],":");
	int readCount(0);
	for(unsigned int i (0) ; i< positionsRead.size() ; i++){
	  //hit.charCount_[i] = atoi(positionsRead[i].c_str());
	  readCount += atoi(positionsRead[i].c_str());
	}
	long bwtPosition = atol(lineVector[2].c_str());
	unsigned short fileNum =  filePositions[whichPile[(int)lineVector[1][0]]][bwtPosition];

	  // findFileNum(bwtPosition, mergCOutput[whichPile[(int)lineVector[1][0]]]);
	//	seqInfo[fileNum].hits_.push_back(hit);
	//	unsigned pos;
	cout << "------------------------------------------" <<endl;
      	//cout << "atoi " << atoi(lineVector[2].c_str()) <<endl;
	//	cout << "p " << whichPile[(int)lineVector[1][0]] <<endl;
	//	unsigned  pos = findSuffixPosition(hit.bwtPosition_, mergeAOutput[whichPile[(int)lineVector[1][0]]]);
	cout << line << endl;
	cout << "fielnum " << fileNum << endl;
	cout <<  " " << seqInfo[fileNum].fileName_ <<  "  " << seqInfo[fileNum].tag_ << endl;
	/*if(seqInfo[fileNum].fileName_.find("rev") == string::npos){
	  
	//	  seqInfo[fileNum].suffixStarts_.push_back(pos);
	  //seqInfo[fileNum].suffixEnds_.push_back(pos + hit.suffix_.length());
        }
        else{
	  //  seqInfo[fileNum].suffixStarts_.push_back(seqInfo[fileNum].sequenceLength_ - (pos + hit.suffix_.length()));
	  // seqInfo[fileNum].suffixEnds_.push_back(seqInfo[fileNum].sequenceLength_ - pos);
	  }*/
	//seqInfo[fileNum].suffixCoverage_.push_back(readCount);
	//	cout << line << endl;
	//	cout << hit.bwtPosition_ << " " << hit.suffix << endl;
      }
      
      // break;                                                                                                                                       
    }
    if(i%1000000==0) cerr << "." <<endl;                                                                                                      
  }
}


void testSigma(vector<SequenceInformation> &seqInfo, FILE* zeroMergA, FILE* zeroMergC){
  //  cout << "A "<< getFileSize(zeroMergA) <<endl;
  //cout << "C " << getFileSize(zeroMergC) << endl;
  //  cout << "testSigma" << zeroMergA.size()<< endl;
  unsigned short fNum;
  unsigned suffixPos;
  cout << "seqInf " << seqInfo.size() << endl;
  int count(0);
  
  cout << "unsigned int " << sizeof(unsigned int)<<endl;
  cout << "unsigned " << sizeof(unsigned) << endl;
  vector<unsigned> suffixPosit;
  vector<unsigned short> fileNums;

  while(fread(&fNum, sizeof(unsigned short), 1 , zeroMergC) == 1){
    fileNums.push_back(fNum);
  }
  while(fread(&suffixPos, sizeof(unsigned), 1, zeroMergA) == 1){
    suffixPosit.push_back(suffixPos);
  }
  // cerr << "there is something wrong with the suffixe at position " << fNum << endl;
  //break;
  
  cout << "sufixpos " << suffixPosit.size()<<endl;
  cout << "fileNums " << fileNums.size()<<endl;  
  count++;
  bool suffixOK(true);
  for(unsigned int i(0); i< fileNums.size() ; i++){
    //    cout <<endl << "------" <<endl;
    if(seqInfo[fileNums[i]].sequenceLength_ != suffixPosit[i]){
      cout <<"------------------------------------------------" <<endl;
      cout << "testSigma " << seqInfo[fileNums[i]].fileName_ << "\t\t" << seqInfo[fileNums[i]].sequenceLength_ <<endl;
      cout << "testSigma " << i << "\t\t" << suffixPosit[i]<< endl;
      suffixOK = false;
    }
  }
  if(!suffixOK) 
    cerr<< "There is something wrong with either the Input Files or the Sigma Suffix Positions" <<endl;
  else
    cerr << "Sigma suffixe seem to be ok" <<endl;
}


unsigned short findFileNum(long BWTPosition, FILE* mergeCPile){
  long byteOff = (BWTPosition) * sizeof(unsigned short);
  //  cout <<"offset " <<  byteOff <<endl;
  //  cout <<"fileSize " << getFileSize(mergeCPile)<<endl;
  fseek(mergeCPile, byteOff, SEEK_SET);
  unsigned short fNum;
  fread(&fNum, sizeof(unsigned short),1, mergeCPile);
  fseek(mergeCPile, 0, SEEK_SET);
  //  cout << fNum << " file " << seqInfo[fNum].fileName_ <<endl; 
  return fNum;
}

unsigned findSuffixPosition(long BWTPosition, string mergeAPile){
  //  cout << "find suffix " << BWTPosition <<endl;
  
  FILE* merg  = fopen(mergeAPile.c_str() , "r");
  fseek(merg,(BWTPosition * sizeof(unsigned)), SEEK_SET);
  unsigned sufPos;
  fread(&sufPos, sizeof(unsigned),1,merg);
  fclose(merg);
  //  cout << "return " <<sufPos <<endl;
  return sufPos;
}


vector<SequenceInformation> loadSequenceInfoToParse(string sigmaSuffixOut, string headerFile, string fileCounter, string mergCZeroOut){
  //first generate the SequenceInformation
  vector<SequenceInformation> seqsInfo;
  SequenceInformation oneSeq;
  FILE* cZeroOut= fopen(mergCZeroOut.c_str(), "r");
  FILE* aZeroOut= fopen(sigmaSuffixOut.c_str(),"r");
  ifstream fileCount(fileCounter.c_str(),ios::in);
  assert(cZeroOut != NULL);
  assert(aZeroOut != NULL);
  string line;
  string fName;
  unsigned short fileNum;
  unsigned suffixPosition; //suffix sigma position should be the same as the length of the sequence
  int fCount(0);
  while(fread(&fileNum, sizeof(unsigned short), 1, cZeroOut) == 1){
    if(fread(&suffixPosition,sizeof(unsigned),1,aZeroOut) == 1){
      if(fileCount.good()){
	getline(fileCount, line);
	unsigned short fileCNum = (unsigned short) strtoul(line.substr(0,line.find(",")).c_str(),NULL,0);
	if(fileCNum != fileNum|| fileCNum != fCount){
	  cerr<< "something wrong with the fileNumber" <<endl;
	  break;
	}
	fName = line.substr(line.find(",")+1 );
	oneSeq.fileNum_ =fileCNum;
	oneSeq.sequenceLength_ = suffixPosition;
	oneSeq.fileName_ = fName;
	seqsInfo.push_back(oneSeq);
	fCount++;
      }
      else{
	cerr <<" file count done before End of suffix pos "<<endl;
      }
    }else{
      cerr <<"suffixes done before c count " <<endl;
    }
  }
  fileCount.close();
  ifstream head(headerFile.c_str() ,ios::in);
  int tagCount(0);

  while(head.good()){
    getline(head,line);
    string fileCount = line.substr(0,line.find(","));
    string tag = line.substr(line.find(",") + 1);
    string f = "G_" + fileCount;
    string fr = "G_"+ fileCount + "_rev";
    int revCount(0);
    for(unsigned int i (0); i<seqsInfo.size(); i++) {      
      if(seqsInfo[i].fileName_.compare(f) ==0 || seqsInfo[i].fileName_.compare(fr) == 0){
        revCount++;
	//	cout << seqsInfo[i].fileName_ << " f " << f << " fr "  <<fr <<endl; 
	seqsInfo[i].tag_ = tag;
	tagCount++;
	if(revCount == 2)
	  break;
      }
    }
  }
  //  for(unsigned int i (0) ; i< seqsInfo.size(); i++)
  //cout << i << " s " << seqsInfo[i].fileName_ << " " << seqsInfo[i].tag_ <<endl;

  cout << "seqInf " << seqsInfo.size() << " tagCount " << tagCount <<endl;
  head.close();
  return seqsInfo;
}


vector<SequenceInformation> loadSequenceInfoToTest(vector<string> sequenceFiles, vector<string> mergingOutput){
  FILE* mergeZero = NULL;
  vector<SequenceInformation> seqInfo;
  seqInfo.resize(sequenceFiles.size());
  SequenceInformation oneSeq;
  // seqInfo.resize(sequenceFiles.size());
  cout << "sequence File count "<< sequenceFiles.size() << endl;
  vector<unsigned short> posVector;
  //just to test the amout of input files
  
  for(unsigned int i(0); i< mergingOutput.size(); i++) {
    FILE* merg = fopen(mergingOutput[i].c_str(), "r");
    //    cout <<"open " << mergingOutput[i]<<endl;
    assert(merg!=NULL);
    if(getFileSize(merg) % sizeof(unsigned short) != 0){
      cerr << "MergfileC " << mergingOutput[i] << " has the wrong size" << endl;
    }
    if(mergingOutput[i].find("-C00")!=string::npos){
      cout << "found it " << mergingOutput[i] << endl;
      mergeZero = fopen(mergingOutput[i].c_str(),"r");
      break;
    }
  } 
  assert(mergeZero!=NULL);
  unsigned short fNum ;
  while(fread(&fNum, sizeof(unsigned short),1,mergeZero) == 1){
    posVector.push_back(fNum);
  }
  cout << "positionSize " << posVector.size() << " seq length " << sequenceFiles. size() << endl; 

  assert(posVector.size() == sequenceFiles.size());

  fclose(mergeZero);
  for(unsigned int i(0); i< sequenceFiles.size(); i++){
    unsigned int sequencelength = 0; 
   ifstream fasFile(sequenceFiles[i].c_str(), ios::in);
    string line;
    
    while(fasFile.good()){
      getline(fasFile,line);
      if(line[0]=='>')
	oneSeq.tag_ = line.substr(1,line.length()-1);
      else
	sequencelength += line.length();
    }
    oneSeq.fileNum_ = (unsigned short) i;
    int found;
    found = sequenceFiles[i].find_last_of("/\\");
    
    oneSeq.fileName_ = sequenceFiles[i].substr(found+1);
    oneSeq.sequenceLength_ = sequencelength;
    
    seqInfo[(unsigned short) i] =  oneSeq;
  }
  cout << "SeqIno " << seqInfo.size()<< endl;
  return seqInfo;
}

void testSuffixe(vector<SequenceInformation> &seqs, vector<FILE*> suffixFiles, vector<FILE*> mergeCFiles){

  for(unsigned int i(0) ; i<mergeCFiles.size() ;i++){
    vector<unsigned> suffixPositions;
    vector<unsigned short> fileCounts;
    unsigned short fileCount;
    unsigned sufPos;
    while(fread(&fileCount,sizeof(unsigned short),1, mergeCFiles[i]) == 1){
      fileCounts.push_back(fileCount); 
      if(fread(&sufPos,sizeof(unsigned), 1, suffixFiles[i]) == 1){
	suffixPositions.push_back(sufPos);
      }
      else{
	cerr << "There are not enough Suffixpositions " << endl;
      }
    }
    cout << "Test suffixe " << suffixPositions.size() << endl;
    cout << "Test fileCount " << fileCounts.size() << endl;
    bool suffixOK(true);
    for(unsigned int j(0); j< suffixPositions.size(); j++){
      if(suffixPositions[j] > seqs[fileCounts[j]].sequenceLength_){
	cerr << "Position "<<suffixPositions[j] << " is higher than " << seqs[fileCounts[j]].sequenceLength_ << " in file " << seqs[fileCounts[j]].fileName_ << endl;
	suffixOK = false;
      }
    }
    if(!suffixOK) 
      cerr << "There is something wrong with the suffixe in " << i << endl;
    else
      cerr << "Everything ok with suffixe in "<< i <<endl;
    fseek(suffixFiles[i],0,SEEK_SET);
    fseek(mergeCFiles[i],0,SEEK_SET);
  }
  //second test this time test the fseek reading of single suffixe
}

void testCountWordOutput(vector<SequenceInformation> &seqs, vector<FILE*> mergeCFiles, vector<string> mergeSuffixFiles){
  vector<unsigned> allSuffixPositions[mergeCFiles.size()]; 
  vector<unsigned short> allFilePositions[mergeCFiles.size()] ;
  unsigned short fileCount;
  unsigned sufPos;
  for(unsigned int i(0) ; i<mergeCFiles.size() ;i++){
    cout << "testcoundword " << i <<endl;
    
    fseek(mergeCFiles[i], 0, SEEK_SET);
    FILE* suff = fopen(mergeSuffixFiles[i].c_str(),"r");
    fseek(suff, 0, SEEK_SET);
    cout << "size c " << getFileSize(mergeCFiles[i]) <<endl;
    cout << "size s " << getFileSize(suff) <<endl; 
    vector<unsigned> suffixPositions;
    vector<unsigned short> fileCounts;
    while(fread(&fileCount,sizeof(unsigned short),1 ,mergeCFiles[i]) == 1){
      fileCounts.push_back(fileCount);
      if(fread(&sufPos,sizeof(unsigned), 1, suff) == 1){
	suffixPositions.push_back(sufPos);
      }
      else{
	cerr << "There are not enough Suffixpositions " << endl;
      }
    }
    cout << "SuffixPositions " << i << " " << suffixPositions.size() << endl;
    allSuffixPositions[i] = suffixPositions;
    allFilePositions[i]  = fileCounts;
    fseek(mergeCFiles[i], 0, SEEK_SET);
    fseek(suff, 0, SEEK_SET);
  }
  bool allOk(true);
  cout << "got file and suffix information " << endl;
  for(unsigned int j(0); j< seqs.size(); j++){
    cout << seqs[j].fileName_ <<endl;
    for(unsigned int i(0) ; i<seqs[j].hits_.size(); i++){
      
      long bwtPos = seqs[j].hits_[i].bwtPosition_;

      int pile = whichPile[(int)seqs[j].hits_[i].suffix_[0]];

      unsigned int sufixPosFseek = findSuffixPosition(bwtPos, mergeSuffixFiles[pile]);
      cout << "S suffixStart " << sufixPosFseek << endl;
      if(sufixPosFseek != allSuffixPositions[pile][bwtPos]){
	cerr << "positions gained with fseek not the same as vector " << sufixPosFseek << " " << allSuffixPositions[pile][bwtPos] << endl;
	allOk = false;
      }
      if(bwtPos > allSuffixPositions[pile].size()){
	cerr << "BWTPos " << bwtPos << " to Big for " << pile <<" " << allSuffixPositions[pile].size() <<endl; 
	allOk = false;
	break;
      }
      unsigned sufixPos = allSuffixPositions[pile][bwtPos];
      if(sufixPos > seqs[j].sequenceLength_){
	cerr << "something wrong with position " << sufixPos << " in sequence " << seqs[j].fileName_;
	allOk = false;
      }
    }
  }
  if(allOk)
    cout << "Everything seems to be in order " << endl;
}
  

void printUsage(string message){
  cerr<< "Usage of metagenomics classifier" << endl
      << message<< endl;
  cerr<< "Usage: met -t -mC <Files> -mI <Files> -mA <Files> "<<endl
      << "or: met -p -o File -mC <Files> -mF File -mA0 File -fH File"<<endl 
      << "or: met -p -t -o File -mC <Files> -mF <File> -mC <Files> -mA <Files>" <<endl
      << "or: met -c -o File -mC <Files> -mA <Files> -mF File -fH File." << endl <<endl;

  
  cerr  << "-t Test if the merging of the files went correct and all sequences given in mI are in the right order."<<endl
	<< "-p Classify the countWord output corresponding to the reference sequence" <<endl
	<< "-c Parses the count Word output for the coverage of the Genome/Gene in the Database" <<endl
	<< " If both -t and -c flags are given, the algorithm will test the merging output and abbort when there is something wrong with it." << endl
	<< endl
	<< "-o outPut of countWords. Must be given with the -c option"<< endl
	<< "-mC ouPut C of the merging" << endl
	<< "-mI Input sequences of merging, must be in same order as merging. Must be given with the -t option." <<endl
	<< "-mF FileCounter Output of the merging " <<endl
	<< "-mA output A of the merging" <<endl 
	<< "-mA0 Output A00 of the merging files" <<endl
	<< "-fH Header output of renaming of the sequence files" <<endl  ; 
  
}


vector<string>  split (string s, string token){
  vector<string> vs;
  while(s.find(token) != string::npos){
    vs.push_back(s.substr(0, s.find(token)));
    s = s.substr(s.find(token)+1);
  }
  vs.push_back(s);
  return vs;
}
