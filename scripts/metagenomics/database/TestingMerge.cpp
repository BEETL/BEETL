#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<cassert>
#include<sstream>

using namespace std;

unsigned short findFileNum(long BWTPosition, FILE* mergeCPileOut);

unsigned int findSuffixPosition(long BWTPosition, FILE* mergAPileOut);



struct SequenceInformation
{
  unsigned short fileNum_;
  string sequence_;
  string fileName_;
};

vector<SequenceInformation> loadSequenceInformationWithTest( FILE* mergC0File, vector<string> sequenceFiles);

void testMerging(vector<SequenceInformation>& seqInfo, FILE* mergCFile, FILE* mergAFile);


int main(int argc, char** argv){
  vector<string> sequenceFiles;
  vector<FILE*> mergCFiles;
  vector<FILE*> mergAFiles;
  for(unsigned int i(1) ; i<argc;i++){
    if(strcmp(argv[i], "-sI") == 0){
      for(unsigned int j(i+1); j< argc; j++){
	if(argv[j][0] == '-')
	  break;
	sequenceFiles.push_back(argv[j]);
      }
    }
    else if(strcmp(argv[i], "-mC")==0){
      for(unsigned int j(i+1); j <argc;j++){
	if(argv[j][0] =='-')
	  break;
	FILE* mergC = fopen(argv[j], "r");
	if(mergC != NULL){
	  cout << "mergC " << argv[j] << endl;
	  mergCFiles.push_back(mergC);
	}
      }
    }
    else if(strcmp(argv[i], "-mA") == 0){
      for(unsigned int j(i+1); j<argc;j++){
	if(argv[j][0] == '-')
	  break;
	FILE* mergA =fopen(argv[j],"r");
	if(mergA != NULL){
	  mergAFiles.push_back(mergA);
	  cout <<"mergA "<< argv[j] << endl;
	}
      }
    }
  }
  vector <SequenceInformation> seqs = loadSequenceInformationWithTest(mergCFiles[0], sequenceFiles);  
  cout << "got seq info " << endl;
  for(unsigned int i (0); i < mergAFiles.size(); i++){
    cout << " Test Files " <<i <<endl;
    testMerging(seqs, mergCFiles[i], mergAFiles[i]);
  }
}

long getFileSize(FILE *file)
{
  long lCurPos, lEndPos;
  lCurPos = ftell(file);
  fseek(file, 0, SEEK_END);
  lEndPos = ftell(file);
  fseek(file, lCurPos, SEEK_SET);
  return lEndPos;
}

void testMerging(vector<SequenceInformation> &seqInfo, FILE* mergCFile, FILE* mergAFile){
  long fileSize = getFileSize(mergCFile);
  unsigned short fileNumOld = 30000;
  unsigned long suffixPosOld =-1;
  for(unsigned long i(0); i < fileSize /sizeof(unsigned short); i++){
    unsigned short fileNum = findFileNum(i,mergCFile);
    //    cout <<"f " <<  fileNum <<endl;
    unsigned long suffPos = findSuffixPosition(i, mergAFile);
    //cout << suffPos << " seq length " << seqInfo[fileNum].sequence_.length() << endl;
    if(fileNumOld != 30000){
      string current = seqInfo[fileNum].sequence_.substr(suffPos,seqInfo[fileNum].sequence_.length()-suffPos);
      string old = seqInfo[fileNumOld].sequence_.substr(suffixPosOld, seqInfo[fileNumOld].sequence_.length()- suffixPosOld); 
      //cout << "current " << current<< endl;
      //cout << "old     " <<old <<endl;
      if(current.compare(old) == 0 && fileNum < fileNumOld){
	cout << "File number bigger file " << fileNum << " sp " << suffPos << " bigger than " << fileNumOld << " sp " <<suffixPosOld << endl;
      }
      // AA.compare(AT) = -1
      else if(current.compare(old) < 0){ 
	cout << endl << "o " << old.substr(0,100) <<endl;
	cout << "n " << current.substr(0,100) <<endl;
	cout << "file " << fileNum << " sp " << suffPos << " bigger than " << fileNumOld << " sp " << suffixPosOld << endl;
      }
    }
    fileNumOld = fileNum;
    suffixPosOld = suffPos;
  }
}

vector<SequenceInformation> loadSequenceInformationWithTest(FILE* mergCFile, vector<string> sequenceFiles){
  FILE* mergeZero = NULL;
  vector<SequenceInformation> seqInfo;
  seqInfo.resize(sequenceFiles.size());
  SequenceInformation oneSeq;
  // seqInfo.resize(sequenceFiles.size());                                                                                                              
  cout << "sequence File count "<< sequenceFiles.size() << endl;
  vector<unsigned short> posVector;

  unsigned short fNum ;
  while(fread(&fNum, sizeof(unsigned short),1,mergCFile) == 1){
    posVector.push_back(fNum);
  }
  fseek(mergCFile, 0, SEEK_SET);
  cout <<"pos " << posVector.size() <<" " << sequenceFiles.size() <<endl;
  //  assert(posVector.size() == sequenceFiles.size());
  
  cout <<sequenceFiles.size()<<endl;

  for(unsigned int i(0); i< sequenceFiles.size(); i++){
    if(posVector[i] != (unsigned short) i)
      cerr << "There is something wrong with the files "<< sequenceFiles[i] << " sigma " <<posVector[i-1]  <<endl;
    
    string sequence = "";
    ifstream fasFile(sequenceFiles[i].c_str(), ios::in);
    string line;
    while(fasFile.good()){
      getline(fasFile,line);
      if(line[0]!='>'){
	sequence += line;
      }
    }
    oneSeq.fileNum_ = (unsigned short) i;
    int found;
    found = sequenceFiles[i].find_last_of("/\\");
    oneSeq.fileName_ = sequenceFiles[i].substr(found+1);
    oneSeq.sequence_ = sequence;
    //    cout << i << endl;
    seqInfo[i] = oneSeq;
    //    cout<< seqInfo[i].sequence_.length() <<endl;
  }
  cout << "return seq " <<endl;
  return seqInfo;
}

unsigned short findFileNum(long BWTPosition, FILE* mergeCPile){
  long byteOff = (BWTPosition) * sizeof(unsigned short);
  fseek(mergeCPile, byteOff, SEEK_SET);
  unsigned short fNum;
  fread(&fNum, sizeof(unsigned short),1, mergeCPile);
  fseek(mergeCPile, 0, SEEK_SET);
  //  cout << fNum << " file " << seqInfo[fNum].fileName_ <<endl;
  return fNum;
}

unsigned findSuffixPosition(long BWTPosition, FILE* merg){
  //  cout << "find suffix " << BWTPosition <<endl; 
  fseek(merg,(BWTPosition * sizeof(unsigned)), SEEK_SET);
  unsigned sufPos;
  fread(&sufPos, sizeof(unsigned),1,merg);                                                                       
  return sufPos;
}
