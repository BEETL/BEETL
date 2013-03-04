#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<cassert>
#include<sstream>

#include <string.h>
#include <cstdlib>

using namespace std;

int printUsage(string message);
 
struct Overlap
{
  unsigned int start_;
  unsigned int end_;
  unsigned int suffixCount_;
  unsigned int readCount_;

};

struct BWTInformation
{
  unsigned short pileNum_;
  double readCount_;
  //  unsigned double readAverage_;
  unsigned short wordLength_;
  unsigned int taxId_;
  unsigned short taxLevel_;
  vector< unsigned int > charBCount_;
  vector< unsigned int > charACount_;
  vector<unsigned short> fileNumbers_;
  vector<int> readIds;
};

struct TaxInformation
{
  unsigned short taxLevel_;
  vector<int> files_;
  double *wordCountPerSize_;
  string name_;
  vector <unsigned long>  bwtPositions_;
  vector <double>  wordCounts_;
  vector <unsigned short> wordLengths_;
  vector <unsigned short> pileNumbers_;
  unsigned int parentId_;
  double normalisedCount_;
  unsigned long seqLengthSum_;
};

struct FileInformation
{
  vector<unsigned> suffixPos_;
  vector< int > suffixCounts_;
  vector<unsigned short > suffixLengths_;
  vector<unsigned long > bwtPositions_;
  int sequenceLength_;
  vector<char> suffixChar_;
};

int countForTaxLevel[7];


typedef map<int,TaxInformation> TAXMAP;
typedef map<int,FileInformation> FILEMAP;

typedef map <unsigned long, BWTInformation > BWTMAP;

typedef map <unsigned int, string> READMAP;
long getFileSize(FILE *file)
{
  long lCurPos, lEndPos;
  lCurPos = ftell(file);
  fseek(file, 0, 2);
  lEndPos = ftell(file);
  fseek(file, lCurPos, 0);
  return lEndPos;
}

int whichPile(char c);
int getTaxonomicLevel(string s);
double readCount(0.0);

string levelNames[] = {"superkingdom","phylum" ,"class", "order", "family", "genus", "species", "strain"};

vector<string>  split (string s, string token){
  vector<string> vs;
  while(s.find(token) != string::npos){
    vs.push_back(s.substr(0, s.find(token)));
    s = s.substr(s.find(token)+(token.length()));
  }
  vs.push_back(s);
  return vs;
}

unsigned int alphabetSize(7);
unsigned int taxLevelSize(8);

TAXMAP loadTaxInformationForDatabase(string taxFile, int cycleSize, string ncbiNames);

READMAP loadReadInformation(string seqFile);

void getRefInfoOfCount(TAXMAP &taxInfo, string parseWordOut, vector<int> wordSize);
//void getRefInfoOfBinarys(TAXMAP &taxInfo, vector<string> wordCountOut, int cycleSize);
void printTaxInformation(TAXMAP &taxInfo, FILEMAP &fileInfo);

void parseForCoverage(string countWordOutput, vector<FILE*> mergeAOutput, FILEMAP &fileInfo);

vector<BWTMAP> getBWTInformation(string parseWordOut, int minWordLength, vector<int> wordSizes, string mergedZeroFile, READMAP readMap);

void getTaxCountThroughBWTInfo(vector<BWTMAP> &bwtInfo, TAXMAP &taxInfo, vector<int> cycleSize);

void countWordTest(TAXMAP &taxInfo, vector<int> taxIDsInTest, vector<double> testExpections);

void printTaxTree(TAXMAP &taxInfo, FILEMAP &fileInfo, vector<int> wordSize);

void loadFileNumToTaxIds(string taxIdNames);

vector<BWTMAP> getBWTInformationThroughOnceParsed(string countWordOutput, vector<int> WordSizes);

vector< vector< int> > fileNumToTaxIds;

vector<bool> intervalInSameTaxa(vector<unsigned int>& sharedTaxIds, vector<unsigned short>& fileNumbers);

void getSecondaryInformation(string parsedWordCountOutput);


int main(int argc, char **argv){
  if(argc < 2)
    return printUsage("not enough arguments");
  bool referenceInfo(false);
  bool fastAnalysis(false);
  bool testTaxCount(false);
  bool secondParsed(false);
  bool countTaxLevel(false);
  bool getCoverageInformation(false);

  string taxInfo;
  string parseCount;
  vector<FILE*> mergeCOutput;
  vector<FILE*> mergeAOutput;
  string countWordOutput;
  string ncbiTaxonomyNames;
  string outputFile;
  string mergedZeroFile;
  vector<int> taxIds;
  vector<double> expections;
  string ncbiNames;
  vector<int> wordSize;
  string sequenceFile;
  for(unsigned int i(1); i<argc; i++) {
    if(strcmp(argv[i], "-f") ==0)
      fastAnalysis = true;
    else if(strcmp(argv[i], "-c")==0)
      getCoverageInformation = true;
    else if(strcmp(argv[i], "-p")==0)
      secondParsed = true;
    if(strcmp(argv[i], "-t")==0)
      taxInfo = argv[i+1];
    if(strcmp(argv[i],"-b") ==0)
      countWordOutput = argv[i+1];
    if(strcmp(argv[i], "-ids") ==0) {
      for(unsigned int j(i+1); j< argc; j++){
	if(argv[j][0]== '-')
	  break;
	taxIds.push_back(atoi(argv[j]));
      }
    }
    if(strcmp(argv[i],"-exp" ) ==0){
      for(unsigned int j(i+1); j< argc; j++){
	if(argv[j][0] == '-')
	  break;
	expections.push_back(atof(argv[j]));
      }
    }
    if(strcmp(argv[i], "-w" ) == 0) {
      for(unsigned int j(i+1); j <argc ; j++){
	if(argv[j][0] =='-')
	  break;
	wordSize.push_back(atoi(argv[j]));
      }
    }
    if(strcmp(argv[i],"-m" ) == 0){
      for(unsigned int j(i+1); j < argc; j++){
	if(argv[j][0] == '-')
	  break;
	FILE* mergeAOut = fopen(argv[j],"r" );
	if(mergeAOut == NULL){
	  cerr << "Can not read " << argv[j] << endl;
	  return 1;
	}
	mergeAOutput.push_back(mergeAOut);
      }
    }
    
    if (strcmp(argv[i], "-n") == 0)
      ncbiNames = argv[i+1];
    if(strcmp(argv[i], "-z" ) == 0)
      mergedZeroFile = argv[i+1];
    if(strcmp(argv[i] ,"-s") == 0)
      sequenceFile = argv[i+1];
  }
  if(referenceInfo){
    cerr << "Parsing the reference information" <<endl;
    loadFileNumToTaxIds(taxInfo);
    map<int, TaxInformation> taxInformation = loadTaxInformationForDatabase(taxInfo, wordSize.size(), ncbiNames);
    // getSequenceLengths(taxInformation, mergAOutput[0]);
    
    cerr << "taxinfo "  <<taxInformation.size() <<endl;
    //    getRefInfoOfBinarys(taxInformation, countWordOutput);
    getRefInfoOfCount(taxInformation,countWordOutput, wordSize);
  }
  else if(getCoverageInformation){
    cerr << "getcoverageInformation " << endl;
    FILEMAP fileInfo;
    parseForCoverage(countWordOutput, mergeAOutput, fileInfo);
  }
  else if(countTaxLevel) {
    getSecondaryInformation(countWordOutput);
  }

  else if(fastAnalysis){
    TAXMAP taxInformation = loadTaxInformationForDatabase(taxInfo, wordSize.size(), ncbiNames);
    // for (TAXMAP::iterator it = taxInformation.begin(); it != taxInformation.end(); it++){
      //      if( (*it).second.taxLevel_ == 5)
      //cout << (*it).first << "\t" <<(*it).second.files_.size() << endl;
      //}
    int minWordLength = 300;
    for(int i (0); i <wordSize.size(); i++ ){
      cerr <<"Searching for length " << wordSize[i] << endl;
     minWordLength = (minWordLength > wordSize[i] ) ? wordSize[i] : minWordLength;
    }
    READMAP readMap =  loadReadInformation(sequenceFile);
    vector<BWTMAP> bwtInfo = getBWTInformation(countWordOutput, minWordLength, wordSize, mergedZeroFile, readMap);
    cerr <<"start getting taxa"<<endl;
    getTaxCountThroughBWTInfo(bwtInfo, taxInformation, wordSize);
    

    cerr << "Simple Analysis done " <<endl;
   
    cerr << "Got all file Infos " <<endl;
    FILEMAP fileInfo;
    printTaxTree(taxInformation, fileInfo, wordSize);
    cerr << "printed taxTree" << endl;
  }
  else if(secondParsed){
    cerr <<"Expecting already parsed BWT output" <<endl;
    TAXMAP taxInformation = loadTaxInformationForDatabase(taxInfo, wordSize.size(), ncbiNames);
   
    vector<BWTMAP> bwtInfo = getBWTInformationThroughOnceParsed(countWordOutput,  wordSize);
    cerr <<"got bwtinfo " << bwtInfo.size() << endl;
    getTaxCountThroughBWTInfo(bwtInfo, taxInformation, wordSize);
    FILEMAP fileInfo;
    printTaxTree(taxInformation, fileInfo, wordSize);
  }
  else if(testTaxCount){
    cout << "Test Tax Count" <<endl;    
    if(taxIds.size() != expections.size()) {
      cerr << "TaxIds and counts have different size Abbort " <<endl;
      return -1;
    }
    TAXMAP taxInformation = loadTaxInformationForDatabase(taxInfo, wordSize.size(), ncbiNames);
    cerr <<"TaxSize " << taxInformation.size() <<endl;
    int minWordLength =300;
    for(int i (0); i <wordSize.size(); i++ )
      minWordLength = (minWordLength > wordSize[i] ) ? wordSize[i] : minWordLength;
    READMAP readMap = loadReadInformation(sequenceFile); 
    vector<BWTMAP> bwtInfo = getBWTInformation(countWordOutput, minWordLength, wordSize, mergedZeroFile, readMap);
    
    getTaxCountThroughBWTInfo(bwtInfo, taxInformation, wordSize);
    cerr << "Got taxa Count " <<endl;
    countWordTest(taxInformation,taxIds,expections);
  }
}

void printTaxTree(TAXMAP &taxInfo, FILEMAP &fileInfo, vector<int> wordMinSize){
  cerr << "print taxa Tree" <<endl;
  cerr << "file info " << fileInfo.size() <<endl;
  for( int s (0); s < wordMinSize.size(); s++){
    cout  <<wordMinSize[s]<<endl;
    stringstream ss;
    ss<< wordMinSize[s];
    ofstream output(ss.str().c_str(), ios::out);
    for(int level =0 ; level <= taxLevelSize ; level++){
      cerr << "l" <<level << endl;
      for(TAXMAP::iterator iter = taxInfo.begin() ; iter != taxInfo.end(); iter++){
	
	if( (*iter).second.taxLevel_ == level){ 
	  output << "TAXA\t" << level << "\t" << (*iter).first  << "\t" << (*iter).second.wordCountPerSize_[s] << "\t" << (*iter).second.name_ << "\t" << (*iter).second.parentId_<<endl ;
	  
	  /*	  if(level == (taxLevelSize -1) && (*iter).second.wordCountPerSize_[s] >30)
	    {
	    vector<int> files = (*iter).second.files_;
	    for(unsigned int j(0);j < files.size(); j++){
	    output << "FILE " << (*iter).second.files_[j] <<endl;		
	    FileInformation file = fileInfo[files[j]];
	    vector<unsigned int> suffixPos = file.suffixPos_;
	    for(unsigned int k(0); k < suffixPos.size(); k++){
	    if(file.suffixLength_[k] >= wordMinSize[s]){
	    output << "SUFF " << suffixPos[k]<< " ";
	    if(suffixPos[k] > file.sequenceLength_)
		      cerr << "suffixPosition is wrong "<< suffixPos[k] << endl;
		      for(unsigned int c(0) ; c <file.charCounts_[k].size(); c++)
		      output << file.charCounts_[k][c] << " ";
		      output << endl;
		      }
		      }
		      }
		      }*/
	}
      }
    }
    output.close();
  }
}


/*
  look in the output of the parsing (for the different word lengths) 
   how much could be classified for each taxonomic level
   this can help to get a statistical analysis on how much could be found
*/
void getSecondaryInformation(string countWordOutput){
  ifstream parsedIn(countWordOutput.c_str(), ios::in);
  vector<unsigned long> taxLevelCount;
  for(unsigned int level (0); level <taxLevelSize; level++)
    taxLevelCount.push_back(0);
  string line;
  while(parsedIn.good()){
    //TAXA 1 976 5118577658
    getline(parsedIn,line);
    cout <<line<<endl;
    if(line.length() >8){
      vector<string> splitLine = split(line, " ");
      taxLevelCount[atoi(splitLine[1].c_str())] += atol(splitLine[3].c_str());
    }
  }
  
  for(unsigned int lev(0); lev <taxLevelSize; lev++)
    cout << levelNames[lev] <<  " " << taxLevelCount[lev] <<endl;
}

void getRefInfoOfCount(TAXMAP &taxInfo, string parseWordOut, vector<int> wordSizes){
  ifstream wordIn(parseWordOut.c_str(), ios::in);
  string line;
  int lineCount(0);
  int wordSize(0);
  int wordIndex(0);
  while(wordIn.good()){
    lineCount++;
    if(lineCount %1000000 ==0)
      cerr<<" . " <<  taxInfo[2].wordCountPerSize_[35];
    getline(wordIn,line);
    //    cout << line << endl;
    int taxId;
    vector<string> countBHit = split(line," " );
    //database information found
    if(line[0] == 'B'){
      int taxLevel = atoi(countBHit[1].c_str());
      if(wordSize != countBHit[4].length()){
	wordSize = countBHit[4].length();
	for(unsigned int i(0); i  < wordSizes.size();i++)
	  if(wordSize = (wordSizes[i])){
	    wordIndex = i;
	    break;
	  }
      }
      vector<unsigned short > fileNumbers ;
      vector<string> fileNumberString = split(countBHit[5], ":");
      //      cout << countBHit[4] <<endl;
      for(unsigned int i(0); i < fileNumberString.size(); i++){
	fileNumbers.push_back(atoi(fileNumberString[i].c_str()));
	//cout << fileNumbers[i]<<endl;
      }
      vector<unsigned int > sharedTaxIds;
      sharedTaxIds.resize(taxLevelSize);     
      
      vector<bool> intervalSameTaxa = intervalInSameTaxa(sharedTaxIds,fileNumbers);
      taxInfo[sharedTaxIds[taxLevel]].wordCountPerSize_[wordSize]++;
      //      cerr <<"level " << taxLevel <<endl;
      for(int i(taxLevelSize-1) ; i >=0 ; i--){
	
	//	cerr << i <<" "<< sharedTaxIds[i] <<endl;
	//	if(intervalInSameTaxa[i])
	  
	if(sharedTaxIds[i] != 0) 
	  break;
      }
      //      else
      //cerr <<"there is something wrong" <<endl;

    } 
  }
  //add the counts for the single taxLevel up
  double classWordSum = 0;
  map<int,int> taxIdToLength;
  for(unsigned int s(0); s<wordSizes.size() ; s++){
    for(int level (taxLevelSize -1) ; level >=-1; level--){
      for(TAXMAP::iterator it = taxInfo.begin(); it != taxInfo.end(); it++){
        //add the count of the children from the taxonomic node to the count of the taxonomic node
        double childrenCount = 0;
        for(TAXMAP::iterator ch = taxInfo.begin(); ch != taxInfo.end(); ch ++){
          if( (*ch).second.parentId_ == (*it).first)
            childrenCount += (*ch).second.wordCountPerSize_[s];
	  
	  (*it).second.wordCountPerSize_[s] += childrenCount;
	}
      }
    }
  }
  // print information
  for(int s(0); s<wordSizes.size() ; s++){
    stringstream ss;
    ss<< wordSizes[s] << "refInfo";
    ofstream output(ss.str().c_str(), ios::out);
    for(int level (taxLevelSize -1) ; level >=-1; level--){
      for(TAXMAP::iterator it = taxInfo.begin(); it != taxInfo.end(); it++){
	if((*it).second.taxLevel_ == level){
	  output << "REF " <<level << " " <<  (*it).first << " " << (*it).second.wordCountPerSize_[s] << " " << (*it).second.seqLengthSum_ << endl;
	}
      }
    }    
    output.close();
  }
}

vector<BWTMAP> getBWTInformationThroughOnceParsed(string parsedBWTOutput, vector<int> wordSizes){
  ifstream largestBWTOutput(parsedBWTOutput.c_str(),ios::in);
  vector<BWTMAP> bwtInfo;
  bwtInfo.resize(wordSizes.size());
  string line;
  int lineCount(0);
  while(largestBWTOutput.good()){
    if(lineCount%100000 == 0)
      cerr <<lineCount << " ";
    lineCount ++;
    getline(largestBWTOutput, line);
    if(line.length() > 10){
      vector<string> splitLine = split(line, " ");
      //      146065189 287 58 1 4 0:0:0:2:2:0: 0:0:0:4:0:0: 2544:2208:2974:4163:
      // bwt = 146065189, taxId = 287, wordLengt = 58, pileNum = 1, 4 = readCount, 0:0:0:2:2:0: = countA (reads), 0:0:0:4:0:0: = countB (reference), 2544:2208:2974:4163: = fileNumbers
      unsigned long bwtPosition = atol(splitLine[0].c_str());
      unsigned int wordLength = atoi(splitLine[2].c_str());
      bool firstBWT = true;
      BWTInformation bwt;
      //save all BWT Positions which are
      for(unsigned int s(0); s < wordSizes.size(); s ++){
	//take the smalles possible bwt information for a bwt
	if(wordLength >= wordSizes[s]
	   && bwtInfo[s].find(bwtPosition) == bwtInfo[s].end() ){
	  //set all the BWT information only once for each newly saved BWT position
	  if(firstBWT){
	    unsigned int taxId = atoi(splitLine[1].c_str());
	    unsigned short pileNum = atoi(splitLine[3].c_str());
	    vector<string> countAstring = split(splitLine[5], ":");
	    vector<unsigned int> countA;
	    int readCount = 0;
	    //      cout <<line <<endl;
	    for(unsigned int i(0); i < countAstring.size(); i++){
	      countA.push_back(atoi(countAstring[i].c_str()));
	      readCount += countA[i];
	    }
	    //cout << readCount <<endl;
	    vector<string> countsBString = split(splitLine[5], ":");
	    vector<unsigned int> countB;
	    for(unsigned int i(0); i< countsBString.size(); i++)
	      countB.push_back(atoi(countsBString[i].c_str()));
	    vector<string> fileNumbersString = split(splitLine[6], ":");
	    vector<unsigned short > fileNumbers;
	    for(unsigned int i(0);i<fileNumbersString.size(); i++)
	      fileNumbers.push_back(atoi(fileNumbersString[i].c_str()));
	    
	    bwt.pileNum_ = pileNum;
	    bwt.readCount_ = readCount;
	    bwt.wordLength_ = wordLength;
	    bwt.taxId_ = taxId;
	    bwt.charACount_ = countA;
	    bwt.charBCount_ = countB;
	    bwt.fileNumbers_ = fileNumbers;
	    firstBWT = false;
	  }
	  bwtInfo[s][bwtPosition] = bwt;
	}
      }
    }
  }
  cerr << "got all bwtInfo " << bwtInfo.size()<<endl;
  return bwtInfo;
}

READMAP loadReadInformation(string seqFile){
  ifstream seqs(seqFile.c_str(),ios::in);
  string line;
  READMAP readMap;
  int seqCount(0);
  while(seqs.good()){
    getline(seqs, line);
    readMap[seqCount] = line;
    seqCount++;
  }
  return readMap;
}

vector<BWTMAP> getBWTInformation(string countWordOutput, int minWordLength, vector<int> wordSizes, string mergedZeroFile, READMAP readMap){
  ifstream wordCount( countWordOutput.c_str(),ios::in);
  string line;
  vector<BWTMAP> bwtInfo;
  bwtInfo.resize(wordSizes.size());
  int lineCount(0);
  bool firstHit =false;
  cerr << "Min wordlength " << minWordLength<<endl;
  int wordSizeReached = wordSizes[0];
  int indexWordSize(0);
  FILE* mergeZero = fopen(mergedZeroFile.c_str(),"r");
  unsigned short fileNum(0);
  unsigned genomeLength;
  vector<unsigned> genomeLengths;
  if(mergeZero ==NULL)
    cerr << "File A-00 is not readable. The reads won't be normalized by the genome Length!" <<endl;
  
  if(mergeZero != NULL){
    while(fread(&genomeLength, sizeof(unsigned), 1, mergeZero) == 1){
      genomeLengths.push_back(genomeLength);
      fileNum++;
    }
  }
  else{
    for(unsigned int i (0); i < 100000; i++)
      genomeLengths.push_back(1);
  }

  while(wordCount.good()){
    lineCount++;
    if((lineCount%10000000) ==0){
      cerr <<lineCount  <<  " bwt  "<< wordSizeReached << " " << bwtInfo[indexWordSize].size() << endl;
    }
    getline(wordCount,line);
    if(line.substr(0,5).compare("MTAXA") ==0){
      vector<string> splitLine = split(line," ");
      if(firstHit || splitLine[3].length() > minWordLength){
	if(!firstHit)
	  cerr <<line << endl;
	firstHit = true;
	//MTAXA 6 562 AAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCC 1455 0:2:0:0:0:0 0:0:1:0:0:0 2247:
	// MTAXA 6 37734 ACTAGGGGTCCA 982822322 0:3:1:1:0:2 0:2:0:0:0:0  299:301:  
      	int readCount = 0;
	unsigned long bwtPosition = atol(splitLine[4].c_str());
	unsigned short taxLevel = (unsigned short) atoi(splitLine[1].c_str());
	int wordLength = splitLine[3].length();
	bool firstBWT = true;
	BWTInformation bwt;
	//save all BWT Positions which are 
	for(unsigned int s(0); s < wordSizes.size(); s ++){
	  //take the smallest possible BWT positions for each interested word length. 
	  //this also means the  highest possible count for these lengths
	  if(wordLength >= wordSizes[s] 
	     && bwtInfo[s].find(bwtPosition) == bwtInfo[s].end() ){
	    wordSizeReached = wordSizes[s];
	    indexWordSize = s;
	    //compute the BWT information only once for each line which fits.
	    if(firstBWT){
	      int taxId = atoi(splitLine[2].c_str());
	      vector<string> readCounts = split(splitLine[5],":");
	      vector<unsigned int> charACount;
	      for(unsigned int i(0); i < readCounts.size(); i++)
		{
		  readCount += atoi(readCounts[i].c_str());
		  charACount.push_back(atoi(readCounts[i].c_str()));
		}
	      vector<string> fileCounts = split(splitLine[6], ":");
	      vector<unsigned int> charBCount;
	      for(unsigned int i(0); i < fileCounts.size(); i++)
		charBCount.push_back(atoi(fileCounts[i].c_str()));
	      vector<string> fileNumbers = split (splitLine[7], ":");
	      vector<unsigned short> fileNums;
	      unsigned long genomeLengthsSum(0) ;
	      for(unsigned int i(0); i<fileNumbers.size() -1; i++){
		unsigned short fileNum = (unsigned short) atoi(fileNumbers[i].c_str());
		fileNums.push_back(fileNum);
		genomeLengthsSum += genomeLengths[fileNum];
	      }
	      double averageLengths = (double) genomeLengthsSum / (double) fileNums.size();
	      int pileNum = whichPile(splitLine[3][0]);
	      for(READMAP::iterator it = readMap.begin(); it != readMap.end(); it++){
		//find the sequence for the suffix
		if( (*it).second.find(splitLine[3]) != string::npos)
		  bwt.readIds.push_back((*it).first);
	      }
	      bwt.pileNum_ = pileNum;	      
	      bwt.readCount_ = (double) readCount/ averageLengths;
	      bwt.wordLength_ = splitLine[3].length();
	      bwt.taxId_ = taxId;
	      bwt.charACount_ = charACount;
	      bwt.charBCount_ = charBCount;
	      bwt.fileNumbers_ = fileNums;
	      bwt.taxLevel_ = taxLevel;
	      firstBWT =false;
	    }
	    bwtInfo[s][bwtPosition] = bwt;
	  }
	}
      
      }
    }
  }
  cerr << "got all bwtInformation " <<bwtInfo.size() << endl;
  stringstream ss; 
  ss << countWordOutput << "_largestBWT" ;
  ofstream bwtOut(ss.str().c_str(), ios::out);
  cerr << "print BWT information in " << countWordOutput<< "_largestBWT" <<endl;
  
  for(unsigned int b(0) ; b < bwtInfo.size();b++){
    cerr << b << " " << bwtInfo[b].size() <<endl;
    for(BWTMAP::iterator it = bwtInfo[b].begin(); it!= bwtInfo[b].end(); it ++){
      BWTInformation bwt = (*it).second;
      bwtOut << (*it).first << " " << bwt.taxId_ << " " << bwt.wordLength_ 
	   << " " << bwt.pileNum_ << " " << bwt.readCount_ << " " ;
      for(unsigned int i(0) ; i< bwt.charACount_.size(); i++)
	bwtOut << bwt.charACount_[i] << ":" ;
      bwtOut << " ";
      for(unsigned int i(0); i< bwt.charBCount_.size(); i++)
	bwtOut << bwt.charBCount_[i] << ":" ;
      bwtOut << " ";
      for(unsigned int i(0); i < bwt.fileNumbers_.size(); i++)
	bwtOut << bwt.fileNumbers_[i] << ":";
      bwtOut<<endl;
    }
  }
  
  bwtOut.close();
  return bwtInfo;
}


/*get the wordCounts for each wordSizes in the taxInformation
 */
void getTaxCountThroughBWTInfo(vector<BWTMAP> &bwtInfo, TAXMAP &taxInfo, vector<int> wordSize){
  vector<double> classWordCount;
  int words(0);
  for(int s (0) ; s < wordSize.size(); s++){
    cerr << "getting tax count for " << wordSize[s] << " " << bwtInfo[s].size() <<endl;
    int minWord = wordSize[s];
    double wordCLevel(0);
    for(BWTMAP::iterator bwtIt = bwtInfo[s].begin(); bwtIt != bwtInfo[s].end(); bwtIt++){
      for(TAXMAP::iterator taxIt = taxInfo.begin(); taxIt != taxInfo.end(); taxIt++){
	if( (*bwtIt).second.taxId_ == (*taxIt).first) {
	  (*taxIt).second.wordCountPerSize_[s] += (*bwtIt).second.readCount_;
	  wordCLevel += (*bwtIt).second.readCount_;
	  break;
	}
      }
    }
    cerr << "overall count for " << minWord << " " << wordCLevel <<endl;
  }
  //add the counts for the single taxLevel up
    /*  double classWordSum = 0;
	for(int s(0); s<wordSize.size() ; s++){
	for(int level (taxLevelSize -1) ; level >=-1; level--){
	for(TAXMAP::iterator it = taxInfo.begin(); it != taxInfo.end(); it++){
	//add the count of the children from the taxonomic node to the count of the taxonomic node
  int childrenCount = 0;
  for(TAXMAP::iterator ch = taxInfo.begin(); ch != taxInfo.end(); ch ++){
  if((*ch).first != 0){
  if( (*ch).second.parentId_ == (*it).first)
  childrenCount += (*ch).second.wordCountPerSize_[s];
  }
  }
  (*it).second.wordCountPerSize_[s] += childrenCount;
  }
  }
  }*/
      
  /*
    map<int,double> taxIdToCount;
  double testSuffixCount(0);
  for(int level (taxLevelSize-1) ; level >= 0; level--){
  double levelSum(0);
  for(TAXMAP::iterator it = taxInfo.begin(); it != taxInfo.end();it++){
  double suffixSum(0);
  int taxId = (*it).first;
  if((*it).second.taxLevel_ == level ){
  //get the tax count
        for(unsigned int bwtPos(0); bwtPos < (*it).second.wordCounts_.size(); bwtPos++){
          if((*it).second.wordLengths_[bwtPos] >= wordSize[0]){
	    suffixSum += (*it).second.wordCounts_[bwtPos];
	    testSuffixCount += (*it).second.wordCounts_[bwtPos];
	  }
        }
        //      cerr <<suffixSum <<endl;
	//get the children counts (which should already be processed
        for(TAXMAP::iterator its =taxInfo.begin(); its != taxInfo.end(); its ++){
          if((*its).second.parentId_ != 0 && (*its).second.parentId_ == taxId){
            suffixSum += taxIdToCount[(*its).first];
            //      cerr << "got parent " <<taxId << " : " << (*its).first<<endl;
          }
        }
        taxIdToCount[taxId] = suffixSum;
	//	if(taxId == 37734 || taxId == 1350)
	if(level ==0)
	  cerr << taxId << "\t" << suffixSum <<endl;
        double countDividedByReadNumber =suffixSum / classWordSum;
	levelSum += countDividedByReadNumber;
	taxInfo[taxId].normalisedCount_ = countDividedByReadNumber;
	cout << "TAX\t" <<taxId << "\t" << suffixSum <<endl;
      }
    }    
    
    }
  
    cerr << testSuffixCount <<" " <<endl;*/
}

void countWordTest(TAXMAP &taxInfo, vector<int> taxIDsInTest, vector<double> testExpections){
  map<int, double> expectedCountsPerId;
  for(unsigned int  i (0); i <taxIDsInTest.size(); i++)
    expectedCountsPerId[taxIDsInTest[i]] = testExpections[i];
  
  //get all possible taxIds and their expected Counts
  for(int level(taxLevelSize-1); level >=0; level--){
    for(TAXMAP::iterator it = taxInfo.begin(); it != taxInfo.end(); it++){
      if((*it).second.taxLevel_ == level){
	double expCount = (expectedCountsPerId.find((*it).first) != expectedCountsPerId.end()) 
	  ? expectedCountsPerId[(*it).first]
	  : 0 ;
	bool inTest (false);
	for(TAXMAP::iterator child = taxInfo.begin(); child != taxInfo.end(); child++){
	  if(((*child).second.parentId_ == (*it).first) 
	     && (expectedCountsPerId.find((*child).first) != expectedCountsPerId.end())){
	    expCount += expectedCountsPerId[(*child).first];
	    inTest =true;
	  }
	}
	if(inTest){
	  expectedCountsPerId[(*it).first] = expCount;
	}
      }
    }
  }
  cout << "expectedCounts " << expectedCountsPerId.size() <<endl;
  map<int,double>::iterator its;
  //  for(its = expectedCountsPerId.begin(); its != expectedCountsPerId.end(); its++)
  //  cout << (*its).first <<"\t" << (*its).second <<endl;
  
  for(int level (taxLevelSize-1); level >=0; level--){
    double wronglyClassified(0);
    double rightlyClassified(0);
    double allClassified(0);
    double overestimated(0);
    double underestimated(0);
    double plainWrong(0);
    for(TAXMAP::iterator it = taxInfo.begin(); it != taxInfo.end(); it++){
      if((*it).second.taxLevel_ == level){
	
	//if the id was in the testDataset
	if(expectedCountsPerId.find((*it).first) != expectedCountsPerId.end()){
	  //if more was found than expected, take the overestimation as wrong, the rest as right
	  if(expectedCountsPerId[(*it).first] <= (*it).second.normalisedCount_) {
	    rightlyClassified +=  expectedCountsPerId[(*it).first];
	    wronglyClassified += (*it).second.normalisedCount_ - expectedCountsPerId[(*it).first];
	    overestimated += (*it).second.normalisedCount_ - expectedCountsPerId[(*it).first];
	    allClassified += (*it).second.normalisedCount_;
	    //    cout << (*it).first << " exp smaller " << rightlyClassified << " " << wronglyClassified << endl;
	  }
	  //if there was less found than expected, take all that was found as correct and the rest of what was not found as wrong
	  else {
	    rightlyClassified += (*it).second.normalisedCount_;
	    wronglyClassified += expectedCountsPerId[(*it).first] - (*it).second.normalisedCount_;
	    underestimated += expectedCountsPerId[(*it).first] - (*it).second.normalisedCount_;
	    allClassified += (*it).second.normalisedCount_;
	    //cout << (*it).first << " exp smaller " << rightlyClassified << " " << wronglyClassified << endl;
	  }
	}
	//take all which was not expected as wrong
	else{
	  plainWrong += (*it).second.normalisedCount_;
	  wronglyClassified += (*it).second.normalisedCount_;
	  allClassified += (*it).second.normalisedCount_;
	  //	  cout << (*it).first << " exp not there " << rightlyClassified << " " << wronglyClassified << endl;
	}
      }
    }
    
    cout << "Level " << level << endl << "wrong " << wronglyClassified<< endl << "right " << rightlyClassified<<endl << "all   " << allClassified<<endl;
    cout << "over " << overestimated << endl << "under " << underestimated << endl << " wrong " << plainWrong <<endl;
  }  
}
   
int whichPile(char c){
  switch(c) 
    {
    case 'A':
      return 1;
      break;
    case 'C':
      return 2;
      break;
    case 'G':
      return 3;
      break;
    case 'N':
      return 4;
      break;
    case 'T':
      return 5;
      break;
    default:
      return 6 ;
    }
}

map<int, TaxInformation> loadTaxInformationForDatabase(string taxInfo, int cycleSize, string ncbiNamesDMP){
  ifstream dbTax(taxInfo.c_str(),ios::in);
  map<int,TaxInformation> taxInformation;
  string line;
  cerr << "loading NCBI taxonomy " <<endl;
  ifstream ncbiNames(ncbiNamesDMP.c_str(), ios::in);
  map<int,string> idToScientificName;
  while(ncbiNames.good()){
    getline(ncbiNames,line);
    vector<string> lineVector = split(line , "\t|\t");
    string name = lineVector[1];
    //    cout << "ncbiName >"<< name <<"<" <<endl;
    int id = atoi(lineVector[0].c_str());
    if(line.find("scientific") != -1)
      idToScientificName[id] = name;
  }
  
  cerr <<"Got NCBI Taxonomy " << idToScientificName.size() << endl;

  while(dbTax.good()){
    getline(dbTax,line);
    vector<string> splitTax = split(line, " ");
    int fileNum = atoi(splitTax[0].c_str());
    for(unsigned int i(1) ;i< splitTax.size(); i++){
      int taxId = atoi(splitTax[i].c_str());
      //cout <<" taxId "  <<taxId <<endl;
      bool taxFound(false);
      for(TAXMAP::iterator it(taxInformation.begin()) ;it != taxInformation.end();it++){
	if((*it).first == taxId){  
	  //	  cout << (i-1) << " "<< countForTaxLevel[i-1] <<endl;
	  taxFound = true;
	  break;
	}
      }
      if(taxFound)
	taxInformation[taxId].files_.push_back(fileNum);
      else{
	TaxInformation tax;
	if(i>1){
	  tax.parentId_ = atoi(splitTax[i-1].c_str());
	  if(tax.parentId_ == 0) 
	    tax.parentId_ =  atoi(splitTax[i-2].c_str());
	  //	  cerr << taxId << " parent " << tax.parentId_ <<endl;
	}
	else //set the root as a parent Id
	  tax.parentId_ = 1;
	
	tax.files_.push_back(fileNum);
	tax.taxLevel_ = i-1;
	tax.name_ = idToScientificName[taxId];
	tax.seqLengthSum_ = 0;
	countForTaxLevel[i-1] +=1;
	//	cout <<countForTaxLevel[i-1] << endl;
	//initialise all word counts with zero
	//cycle 0 starts with two chars
	tax.normalisedCount_ = 0;
	tax.wordCountPerSize_= new double[cycleSize+2];
       
	for(unsigned int i(0) ; i < cycleSize+1;i++)
	  tax.wordCountPerSize_[i] = 0;
	taxInformation[taxId] = tax;
      }
    }
  }
  //set the root
  taxInformation[1].taxLevel_ = -1;
  taxInformation[1].wordCountPerSize_= new double[cycleSize+2];
  taxInformation[1].name_ = "root";
  for(unsigned int i(0) ; i < cycleSize+1;i++)
    taxInformation[1].wordCountPerSize_[i] = 0;
  dbTax.close();
  return taxInformation;
}

int getTaxonomicLevel(string s) {
  //string levelNames[] = {"superkingdom","phylum" ,"class", "order", "family", "genus", "species"};                                                    
  int level;
  if(s.compare("superkingdom") ==0)
    level = 0;
  else if(s.compare("phylum") ==0)
    level = 1;
  else if(s.compare("class") ==0)
    level= 2;
  else if(s.compare("order") == 0)
    level= 3;
  else if(s.compare("family") ==0)
    level= 4;
  else if(s.compare("genus") == 0)
    level= 5;
  else if(s.compare("species") == 0)
    level= 6;
  else
    level= 10;
  return level;
}

// To repair the referenceTEstOutput
vector<bool> intervalInSameTaxa(vector<unsigned int>& sharedTaxIds, vector<unsigned short>& fileNumbers){
  //first get the matching fileNumbers out of the file with the file numbers corresponding to bwt positions of the mergin
  vector<bool> taxSame;
  taxSame.resize(taxLevelSize);
  //look if the files have at each point of the taxonomic tree different tax level or if they are the same
  bool sameTaxa(false);
  cout << "intervalInSame " <<fileNumbers.size() <<endl;
  cout << fileNumToTaxIds.size() <<endl;
  for(unsigned int i(0) ; i< taxLevelSize;i++){
    for(unsigned int j(0); j < fileNumbers.size()-1; j++){
      
      cout << "Number " << fileNumbers[j]<<endl;
      if(fileNumToTaxIds[fileNumbers[j]][i] == fileNumToTaxIds[fileNumbers[j+1]][i]
         && fileNumToTaxIds[fileNumbers[j]][i] != 0
         && fileNumToTaxIds[fileNumbers[j+1]][i] != 0)
        {
	  cout << "true" <<endl;
          sameTaxa = true;
        }
      else
        {
          //if one of the taxa is different than the one before it is enough to set all to false                                                        
          sameTaxa =false;
          break;
        }
    }
    taxSame[i] = sameTaxa;
    if(taxSame[i])
      sharedTaxIds[i] = fileNumToTaxIds[fileNumbers[0]][i];
    else
      sharedTaxIds[i] = 0;
  }
  return taxSame;
}

void loadFileNumToTaxIds(string taxIdNames){
  ifstream taxas(taxIdNames.c_str(), ios::in);
  string line;
  while(taxas.good()){
    vector<int> taxIDs;
    getline(taxas,line);
    if(line.compare("") !=0){
      unsigned int fileNum = atoi(line.substr(0,line.find(" ")).c_str());
      
      line =line.substr(line.find(" " ) +1, line.length());
      while(line.find(" ") != string::npos){
	taxIDs.push_back(atoi(line.substr(0,line.find(" ")).c_str()));
	line = line.substr(line.find(" ") + 1, line.length());
      }
      taxIDs.push_back(atoi(line.c_str()));
      if(taxIDs.size() < taxLevelSize){
	cerr << "Tax Ids have not enough taxonomic Information. Only  "<< taxIDs.size() << " could be found " << endl
	     << "Will add unknown taxa until size is right" <<endl;
	for(unsigned int i(taxIDs.size()-1) ; i< taxLevelSize; i++)
	  taxIDs.push_back(0);
      }
      if(taxIDs.size() > taxLevelSize)
	cerr << "Tax Ids have to much taxonomic information. "
	     << "Please note, that the taxonomic information about one file should be shown as: " <<endl
	     << "FileNumber Superkingdom Phylum Order Family Genus Species" <<endl;
      fileNumToTaxIds.push_back(taxIDs);
      unsigned int test = fileNum +1;
      if(test!= fileNumToTaxIds.size())
	cout<< "Wrong filenumber " << fileNum << " " << fileNumToTaxIds.size() <<endl;
    }
  }
  cout << " fineNumToTaxIds " << fileNumToTaxIds.size() <<endl;
}

void parseForCoverage(string countWordOutput, vector<FILE*> mergeAOutput,FILEMAP &fileInfo){
  ifstream wordCount(countWordOutput.c_str(),ios::in);
  string line;
  int count(0);
  while(wordCount.good()){
    count++;
    if((count%10000) ==0){
      cerr <<count  << " " << fileInfo.size()<<endl;
    }
    int fileNumber;
    FileInformation info;
    getline(wordCount,line);
    //    cout << line<<endl;
    vector<string> splitLine = split(line," ");
    if(splitLine[0].compare("MTAXA") == 0){
      //MTAXA 2 91061 GGCTGCCAACTAA 1197996044 0:0:2:3:0:5 0:0:0:0:0:7  1121:1123:816:77:75:1460:1462:
      unsigned long BWTPosition = atol(splitLine[4].c_str());
     
      vector<string> fileNumbersStrings = (splitLine[7].compare("") == 0) ? split(splitLine[8], ":") : split(splitLine[7],":");
      vector<int> fileNumbers;
      for(uint i(0); i< fileNumbersStrings.size()-1 ; i++){
	fileNumbers.push_back(atoi(fileNumbersStrings[i].c_str()));
	//cerr << "string " << fileNumbersStrings[i] << " int " << atoi(fileNumbersStrings[i].c_str()) <<endl;
      }
      int fileCounts = fileNumbers.size();
      fseek(mergeAOutput[whichPile(splitLine[3][0])], (BWTPosition *sizeof(unsigned)) , SEEK_SET);
      
      unsigned *suffStarts = (unsigned*) malloc ((fileCounts)* ( sizeof(unsigned) ) );
      //      cout <<"got unsigned " <<endl;
      //read the start of the words out of the mergeA-output

      fread(suffStarts, sizeof(unsigned), fileCounts, mergeAOutput[whichPile(splitLine[3][0])]);

      // cout << "Read suff worked " <<endl;
      int readsSuffCount(0);
      vector<string> countsA = split(splitLine[5],":");
      for(uint i(0); i< countsA.size(); i++)
	readsSuffCount += atoi(countsA[i].c_str());
      for(uint i(0); i<fileCounts; i++){
	//if(fileNumbers[i] != fileNum[i] )
	int fileNumber = fileNumbers[i];
	//cerr<<"fileNumber " << fileNumber << " " << splitLine[7] << endl;
	
	fileInfo[fileNumber].suffixPos_.push_back(suffStarts[i]);
	fileInfo[fileNumber].suffixLengths_.push_back(splitLine[3].length());
	fileInfo[fileNumber].suffixCounts_.push_back(readsSuffCount);
	fileInfo[fileNumber].suffixChar_.push_back(splitLine[3][0]);
	fileInfo[fileNumber].bwtPositions_.push_back(BWTPosition);
	//if(fileNumber != fileNum[i])
	// cerr <<"There is something wrong with the fileNumber" <<endl;
      }
      //      cout << count <<endl;
      delete suffStarts;
      //      delete fileNum;
    }
  }
  //get only the single suffix information
  for(FILEMAP::iterator it = fileInfo.begin(); it != fileInfo.end() ; it ++){
    vector<unsigned> uniqueSuffixPos;
    vector< int > uniqueSuffixCounts;
    vector<unsigned short > uniqueSuffixLengths;
    vector<char> uniqueSuffixChar;
    vector<unsigned long > uniqueBWTs;
    if((*it).second.suffixPos_.size() >1)
    for(unsigned int i(0); i < (*it).second.suffixPos_.size() ; i++){
      bool foundSuff (false);
      int suffPosition(0);
      for(unsigned int j(0); j < uniqueSuffixPos.size() ; j++){
	if(uniqueSuffixPos[j] == (*it).second.suffixPos_[i]){
	  foundSuff = true;
	  suffPosition = j;
	  break;
	}
      }
      if(!foundSuff){
	uniqueSuffixPos.push_back((*it).second.suffixPos_[i]);
	uniqueSuffixCounts.push_back((*it).second.suffixCounts_[i]);
	uniqueSuffixLengths.push_back((*it).second.suffixLengths_[i]);
	uniqueSuffixChar.push_back((*it).second.suffixChar_[i]);
	uniqueBWTs.push_back((*it).second.bwtPositions_[i]);
      }
      else if(uniqueSuffixLengths[suffPosition] < (*it).second.suffixLengths_[i]){
	uniqueSuffixLengths[suffPosition] = (*it).second.suffixLengths_[i];
	uniqueSuffixCounts[suffPosition] = (*it).second.suffixCounts_[i];
      }
    }
    /*    if((*it).second.suffixPos_.size() >1){
      cerr << "upos " << uniqueSuffixPos.size() <<endl;
      cerr << "ulen " << uniqueSuffixLengths.size() <<endl;
      cerr << "ucha " << uniqueSuffixChar.size() <<endl;
      cerr << "ucou " << uniqueSuffixCounts.size() <<endl; 
      }*/
    (*it).second.suffixPos_ = uniqueSuffixPos;
    (*it).second.suffixCounts_ = uniqueSuffixCounts;
    (*it).second.suffixLengths_ = uniqueSuffixLengths;
    (*it).second.suffixChar_ = uniqueSuffixChar;
    (*it).second.bwtPositions_= uniqueBWTs;
  }
  //Get the sequence length of the files
  short fileNumberRead;
  unsigned sequenceLength;
  //set the sequence lengths of each file   
  for(FILEMAP::iterator it = fileInfo.begin(); it != fileInfo.end() ; it++){  
    if((*it).second.suffixPos_.size() >1){
      stringstream ss;
      cerr << "print to " << (*it).first<<endl;
      ss<< "F_"<< (*it).first;
      ofstream output(ss.str().c_str(), ios::out);
      cerr << "suffixCount " << (*it).second.suffixPos_.size() <<endl;
      for(unsigned int i(0);i< (*it).second.suffixPos_.size();i++){
	output << (*it).second.suffixPos_[i] << " ";
	output << (*it).second.suffixLengths_[i] << " " ;
	output << (*it).second.suffixCounts_[i] << " ";
	output << (*it).second.suffixChar_[i] << " ";
	output << (*it).second.bwtPositions_[i] << endl;
	
      }
      output.close();
    }
    
  }
}


/*  

vector<SequenceInformation> loadSequenceInfo(vector<string> sequenceFiles, vector<string> mergingOutput){
  FILE* mergeZero;
  vector<SequenceInformation> seqInfo;
  SequenceInformation oneSeq;
  seqInfo.resize(sequenceFiles.size());
  cout << "sequence File count "<< sequenceFiles.size() << endl;
  vector<unsigned short> posVector;
  //just to test the amout of input files
  for(int i(0); i< mergingOutput.size(); i++) {
    FILE* merg = fopen(mergingOutput[i].c_str(), "r");
    if(getFileSize(merg) % sizeof(unsigned short) != 0){
      cerr << "MergfileC " << mergingOutput[i] << " has the wrong size" << endl;
    }
    if(mergingOutput[i].find("-C00")!=string::npos){
      //      cout << "found it " << mergingOutput[i] << endl;
      mergeZero = fopen(mergingOutput[i].c_str(),"r");
    }
  }
   
  assert(mergeZero!=NULL);
  unsigned short fNum ;
  long size = getFileSize(mergeZero);
  
  while(fread(&fNum, sizeof(unsigned short),1,mergeZero) == 1){
    posVector.push_back(fNum);
  }
  cout << "positionSize " << posVector.size() << " seq length " << sequenceFiles. size() << endl; 

  assert(posVector.size() == sequenceFiles.size());

  fclose(mergeZero);
  for(int i(0); i< sequenceFiles.size(); i++){
    unsigned int sequencelength = 0;
    string fastaTag; 
    ifstream fasFile(sequenceFiles[i].c_str(), ios::in);
    string line;
    
    while(fasFile.good()){
      getline(fasFile,line);
      if(line[0]=='>')
	oneSeq.fastaTag_ = line.substr(1,line.length()-1);
      else
	sequencelength += line.length();
    }
    oneSeq.fileNum_ = (unsigned short) i;
    int found;
    found = sequenceFiles[i].find_last_of("/\\");
    
    oneSeq.fileName_ = sequenceFiles[i].substr(found+1);
    oneSeq.sequenceLength_ = sequencelength;
    seqInfo[i] = oneSe+
    
*/


int printUsage(string message){
  cerr<<endl <<message <<endl <<endl;
  cerr<< "Usage: " <<endl
      << "parse -f -b file -t file -n file -z mergedAZeroFile -w <int>" <<endl
      << "parse -p -b file -t file -n file -w <int>"<<endl
      << "parse -c -b file -m <file>  " <<endl
      << endl;
    
  cerr  << "-f\t Fast and simple version, just counts the occurences of the single taxIds and read counts." <<endl
	<< "-p \t -b should contain the bwt info which was already parsed once " <<endl
	<< "-c get the coverage information for each file"<< endl
	<< "-t taxonomic input which was used for the countWord -m call." << endl
	<< "-z A-00 file which is a result of the merging of the genomes. Can be found in the same directory as the B-00 file" <<endl
	<< "-m A-0* files, there is for each BWT-position saved where the corresponding suffix occurs" << endl
	<< "-n  ncbiNames dmp file "<< endl
	<< "-w minimal word counts. Gives the results for everything bigger than the chosen words out." <<endl;
  return 1;
}



