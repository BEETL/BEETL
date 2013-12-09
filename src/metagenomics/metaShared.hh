#include <map>

using namespace std;


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
    double *wordCountPerSizeOfChildren_;
    string name_;
    vector <uint64_t> bwtPositions_;
    vector <double>  wordCounts_;
    vector <unsigned short> wordLengths_;
    vector <unsigned short> pileNumbers_;
    int parentId_;
    double normalisedCount_;
    uint64_t seqLengthSum_;
};

struct FileInformation
{
    vector<unsigned> suffixPos_;
    vector< int > suffixCounts_;
    vector<unsigned short > suffixLengths_;
    vector<uint64_t > bwtPositions_;
    int sequenceLength_;
    vector<char> suffixChar_;
};

int countForTaxLevel[7];


typedef map<int, TaxInformation> TAXMAP;
typedef map<int, FileInformation> FILEMAP;

typedef map <uint64_t, BWTInformation > BWTMAP;

typedef map <unsigned int, string> READMAP;
