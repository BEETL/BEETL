#include<iostream>;
#include<string>;

void printUsage(string message);

int main(int argc, char** argv){
  if(argc < 4){
    printUssage();
    return 1;
  }
  vector<string> inputFiles;
  string savingDir;
  string ncbiTaxonomy;
  
  

}




void printUsage(string message){
  cout << "This is a wrapper to create or update a BWT database of different fasta files." << endl
       << "The database will include files for the BWTs (B01 - B06),"
       << "files with unsigned shorts for the corresponding BWT character (C01-C06) "
       << "and files with unsinged long for the corresponding suffix positions of the BWTs (A01-A06)."<<endl;
  cout << "For this instance this wrapper will only work, if there is the possibility to use qsub" <<endl << endl;
  cout << "Usage: " <<endl
       << "wrapper create -f <fastaFiles> -n directory - s savingDirectory"<<endl;
       
  cout << "To create a BWT database please make sure the fasta files have a header with their taxonomic origin encoded in them." <<endl
       << "Example: >GL397224 Prevotella marshii DSM 16973 genomic scaffold SCAFFOLD11, whole genome shotgun sequence. [Prevotella marshii DSM 16973]"
       << "Please make sure the name corresponding to the ncbi taxonomy. " 
       << "-n directory where the ncbi_taxonomy files are stored (names.dmp, nodes.dmp, merged.dmp)";
}
