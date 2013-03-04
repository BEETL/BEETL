#include<iostream>
#include<fstream>
#include<string>
#include<stdlib.h>
#include<vector>
#include<sstream>
#include<string>
#include<map>

using namespace std;


struct NCBINode
{
  string name_;
  int id_;
  int parentId_;
  int taxLevel_;
  int parentLevel_;
};

struct SequenceInformation
{
  unsigned short fileNum_;
  string tag_;
  string fileName_;
  unsigned long  giNumber_;
  string giString_;
};

vector<string>  split (string s, string token){
  vector<string> vs;
  while(s.find(token) != string::npos){
    vs.push_back(s.substr(0, s.find(token)));
    s = s.substr(s.find(token)+(token.length()));
  }
  vs.push_back(s);
  return vs;
}


vector<SequenceInformation> loadSequenceInfo(string headerFile, string fileCounter);

void findNCBITaxa(vector<SequenceInformation> &seqs, string namesDMP, string nodesDMP, string mergedDMP, string fileNumFile);

map<int, int> getGiTotaxId(string gitoTaxId, vector<int> giNumbers);

int getTaxonomicLevel(string s);

string getAncestors(map<int,NCBINode > idToNode, int id);

//string levelNames[] = {"superkingdom","phylum" ,"class", "order", "family", "genus", "species"};
 
map<int,int> giToTaxId;

int main(int argc, char** argv){
  string namesFile;
  string nodesFile;
  string mergedFile;
  string headerFile;
  string output;
  string fileCountFile;
  string giToTaxIdFile;

  vector<string> inputGenomeFiles;
  if(argc < 5)
    cout << "-nA namesfile, -nO nodesFile, -nM mergedFile, -nG giToTaxIdFile, -h headerFile, -f fileCountFile " <<endl;
  for(unsigned int i(0) ; i < argc ; i++){
    if(strcmp(argv[i], "-nA") ==0)
      namesFile = argv[i+1];
    else if(strcmp(argv[i] , "-nO") ==0)
      nodesFile = argv[i+1];
    else if(strcmp(argv[i] , "-nM") == 0)
      mergedFile = argv[i+1];
    else if(strcmp(argv[i], "-h") ==0)
      headerFile = argv[i+1];
    else if(strcmp(argv[i], "-f") ==0)
      fileCountFile = argv[i+1];
    else if(strcmp(argv[i], "-nG") ==0)
      giToTaxIdFile = argv[i+1];
  }
  vector <SequenceInformation > seqs = loadSequenceInfo(headerFile, fileCountFile);
  vector<int> gis;
  gis.push_back(seqs[0].giNumber_);
  for(uint i(1); i <seqs.size(); i++){
    bool found =false;
    for(uint j(0); j<gis.size(); j++){
      if(gis[j] == seqs[i].giNumber_){
	found =true;
	break;
      }
    }
    if(! found){
      gis.push_back(seqs[i].giNumber_);
      cerr << seqs[i].giNumber_ <<"|" ;
    }
  }
  cout << endl;
  cerr << "Gi Id size " <<gis.size() <<endl;
  giToTaxId = getGiTotaxId(giToTaxIdFile, gis);
  for(uint i (0); i < seqs.size() ; i++){
    bool foundTaxa(false);
    for(map<int,int>::iterator it = giToTaxId.begin() ; it != giToTaxId.end(); it++){
      if((*it).first ==seqs[i].giNumber_)
	foundTaxa =true;
    }
    if(!foundTaxa)
      cout << "NO taxId for " << seqs[i].giNumber_ <<endl;
    
  }
  
  cerr << "got GI IDs " << giToTaxId.size() <<endl;
  findNCBITaxa(seqs, namesFile, nodesFile, mergedFile, output);
  return 0;
}

vector<SequenceInformation> loadSequenceInfo(string headerFile, string fileCounter){
  vector<SequenceInformation> seqsInfo;
  ifstream fileCount(fileCounter.c_str(),ios::in);
  SequenceInformation oneSeq;
  string line;
  string fName;
  while  (fileCount.good()){
    getline(fileCount, line);
    unsigned short fileCNum = (unsigned short) strtoul(line.substr(0,line.find(",")).c_str(),NULL,0);
    oneSeq.fileNum_ = fileCNum;
    fName = line.substr(line.find(",")+1 );
    oneSeq.fileName_ = fName;
    seqsInfo.push_back(oneSeq);
  }
  ifstream head(headerFile.c_str() ,ios::in);
  int tagCount(0);

  while(head.good()){
    getline(head,line);
    string fileCount = line.substr(0,line.find(","));
    string tag = line.substr(line.find(",") + 1);
    string f = "G_" + fileCount;
    string fr = "G_"+ fileCount + "_rev";
    int revCount(0);
    bool foundFile(false);    
    
    for(unsigned int i (0); i<seqsInfo.size(); i++) {
      //      cout << seqsInfo[i].fileName_ << " f " << f << " fr "  <<fr <<endl;
      if(seqsInfo[i].fileName_.compare(f) ==0 || seqsInfo[i].fileName_.compare(fr) == 0){
        revCount++;
	//	cout << seqsInfo[i].fileName_ << " f " << f << " fr "  <<fr <<endl;
	if(i != seqsInfo[i].fileNum_) 
	  cerr << "wrong file number " << endl;
	//1,gi|15604717|ref|NC_000117.1| Chlamydia trachomatis D/UW-3/CX, complete genome
        seqsInfo[i].tag_ = tag;
	vector<string> tagSplit = split(line,"|");
	seqsInfo[i].giNumber_ =  atol(tagSplit[1].c_str());
	
       	seqsInfo[i].giString_ = tagSplit[1];
	tagCount++;
	foundFile = true;
      }
    }
    if(!foundFile){
      cerr << "Found no File for " << line <<endl;
      cerr << "file should be " << f << " or " << fr <<endl;
    }
  }
  cerr << "TagCount " << tagCount <<endl;
  cerr << "Got Sequenceinformation " << seqsInfo.size() <<endl;
  return seqsInfo;
}

map<int,int> getGiTotaxId(string giToTaxIdFName, vector<int> giNumbers){
  ifstream giToTaxFile(giToTaxIdFName.c_str(),ios::in);
  string line;
  map<int,int> giToTaxId;
  int smallestGiNumber = giNumbers[0];
  cout << "all GI numbers " <<giNumbers.size() <<endl;
  for(uint i(0); i<giNumbers.size();i++)
    smallestGiNumber = (smallestGiNumber>giNumbers[i] && giNumbers[i] !=0) ? giNumbers[i] : smallestGiNumber;
  cout << "smallest GI " << smallestGiNumber <<endl; 
  int lineNumber(0);
  while(giToTaxFile.good()){
    getline(giToTaxFile,line);
    //    cout <<line <<endl;
    lineNumber++;
    if(line.length() > 5){
      vector<string> lineVector = split(line,"\t");
      //cerr << "split 1 " << lineVector[1] <<endl;
      //    cerr << "split 0 " << lineVector[0] <<endl;
      int giNumber = atoi(lineVector[0].c_str());
      
      //cerr <<giNumber<<endl;
      for(uint i(0); i< giNumbers.size(); i++){
	if(giNumber == giNumbers[i]){
	  //	cerr<< "foudn gi tax |" <<lineVector[0] <<"|" <<endl;
	  int taxId = atoi(lineVector[1].c_str());
	  giToTaxId[giNumber] = taxId;
	  //cerr << giNumber << " " << taxId <<endl;
	  break;
	}
      }
    }
  }

  cout <<" got all GIs " << giToTaxId.size()<<endl;

  giToTaxFile.close();
  return giToTaxId;
}

void findNCBITaxa(vector<SequenceInformation> &seqs, string namesDMP, string nodesDMP, string mergedDMP, string outputInfo){
  ifstream ncbiNames;
  ncbiNames.open(namesDMP.c_str(), ios::in);
  vector<string> lineVector;
  string line;
  map<string,int> nameToId;
  map<int,string> idToScientificName;
  map<int,NCBINode> idToNode;
  vector<  map<int,vector<unsigned short> > > taxLevelToFileCount;
  taxLevelToFileCount.resize(7);
  
  while(ncbiNames.good()){
    getline(ncbiNames,line);
    lineVector = split(line , "\t|\t");
    string name = lineVector[1];
    //    cout << "ncbiName >"<< name <<"<" <<endl;
    int id = atoi(lineVector[0].c_str());
    nameToId[name] = id;
    if(line.find("scientific") != -1)
      idToScientificName[id] = name;
  }
  ifstream ncbiNodes;
  ncbiNodes.open(nodesDMP.c_str(), ios::in);
  NCBINode node;
  while(ncbiNodes.good()){
    getline(ncbiNodes,line);
    lineVector = split(line,"\t|\t");
    
    int id = atoi(lineVector[0].c_str());
    int parentId = atoi(lineVector[1].c_str());
    node.name_ = idToScientificName[id];
    node.id_ = id;
    if(lineVector.size() > 2){
      int taxLevel =  getTaxonomicLevel(lineVector[2]);
      node.taxLevel_=taxLevel;
    }
    else
      cerr <<line <<endl;
    node.parentId_ = parentId;
    if(idToNode.find(id) != idToNode.end())
      cerr << "already has id " << id << " " << idToScientificName[id] <<endl;
    idToNode[id] = node;
  }
  cerr << "got ncbi names " <<endl;
  //get parent levels, and change the taxLevel under species to strain
  for(map< int, NCBINode >::iterator it = idToNode.begin(); it != idToNode.end(); it++){
    NCBINode parent = idToNode[(*it).second.parentId_];
    if(parent.taxLevel_ == getTaxonomicLevel("species"))
      (*it).second.taxLevel_ = getTaxonomicLevel("strain");
    (*it).second.parentLevel_ = parent.taxLevel_;
  }
  cerr << "node to id " << idToNode.size() <<endl;
  
  /*  ifstream mergedStream( mergedDMP.c_str(), ios::in);
  while(mergedStream.good()) {
    getline(mergedStream, line);
    lineVector = split(line,"\t|\t");

    int oldId = atoi(lineVector[0].c_str());
    int newId = atoi(lineVector[1].c_str());
    cerr << oldId << " " << newId <<endl;

    NCBINode oldNode = idToNode[oldId]; 
    cerr << "name " << oldNode.name_ <<endl;
    idToNode[newId] = oldNode;
    }*/
  
  cerr << "after merged out " <<idToNode.size() <<endl;
  cerr <<"seq size " << seqs.size() <<endl;
  for (uint j (0) ; j< seqs.size(); j++) {
    unsigned short fileNum = seqs[j].fileNum_;
    cerr <<fileNum <<endl;
    //unsigned short fileNum = (unsigned short) strtoul(lineVec[1].c_str(), NULL, 0);
    //    string line = seqs[j].tag_;
    //1,gi|15604717|ref|NC_000117.1| Chlamydia trachomatis D/UW-3/CX, complete genome
    //vector<string> lineVector = split(line,"|");*/
      int gi = seqs[j].giNumber_;
      cerr <<gi << " " << giToTaxId[gi] << endl; 
      if(giToTaxId[gi] != 0){
	cout <<  fileNum;      
	cout << getAncestors( idToNode, giToTaxId[gi]);
      }
  }
}


string getAncestors(map<int,NCBINode> idToNode, int id){
  stringstream tree;
  cerr << "start id " << id <<endl;
  stringstream level(" ");
  int taxIds[7];
  tree<< endl;
  int i(7), j(0);
  while(id != 1){
    NCBINode n = idToNode[id];
    //cout <<"Level 10 reached" <<endl;
    if(n.taxLevel_ !=10){
      level << " " << n.taxLevel_; 
      //  cerr << " i " << i << " l " << n.taxLevel_ << endl;
      while(n.taxLevel_ < i){
	taxIds[j] = 0;
	j++;
	//	tree << " " << 0;
	i--;
      }
      taxIds[j] = id;
      //tree << " " << id;
      i--;j++;
    }
    cerr << "id " <<id << " parent " << n.parentId_ <<endl;
    id = n.parentId_;
  }
  for (int k(7) ; k>-1 ; k--)
    cout << " " << taxIds[k];
  //  cout <<endl;
    if(j!=8){
      //      cerr<< id << " something wrong " <<j << " i " << level.str() << endl;
  //    cerr << tree.str()<<endl;
    }
  return tree.str();
}


int getTaxonomicLevel(string s){
  string levelNames[] = {"superkingdom","phylum" ,"class", "order", "family", "genus", "species", "strain"};
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
  else if(s.compare("strain") ==0)
    level =7;
  else 
    level= 10;
  return level;
}

vector<string>  splitString (string s, string token){
  vector<string> vs;
  while(s.find(token) != string::npos){
    vs.push_back(s.substr(0, s.find(token)));
    s = s.substr(s.find(token)+(token.length()));
  }
  vs.push_back(s);
  return vs;
}
