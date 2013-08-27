#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;


class CharPair: public pair<char, char>
{
public:
    CharPair( const char c1, const char c2 )
    {
        first = c1;
        second = c2;
    }

    bool operator<( const CharPair &rhs ) const
    {
        return first < rhs.first;
    }
};


int main( int argc, char **argv )
{
    if ( argc != 4 )
    {
        cout << "Usage: " << argv[0] << " <main file> <other file> <outfile>" << endl;
        exit( -1 );
    }

    char *filename1 = argv[1];
    char *filename2 = argv[2];
    char *filenameOut = argv[3];

    // Check that both files have the same size
    struct stat stat1, stat2;
    if ( stat( filename1, &stat1 ) )
    {
        perror( "Error reading file 1: " );
        exit( -1 );
    }
    if ( stat( filename2, &stat2 ) )
    {
        perror( "Error reading file 2: " );
        exit( -1 );
    }
    if ( stat1.st_size != stat2.st_size )
    {
        cerr << "Error: the 2 files don't have the same size" << endl;
        exit ( -1 );
    }

    size_t fileSize = stat1.st_size;
    vector<char> buf1( fileSize );
    vector<char> buf2( fileSize );

    ifstream is1( filename1 );
    ifstream is2( filename2 );

    is1.read( &buf1[0], fileSize );
    is2.read( &buf2[0], fileSize );

    vector<CharPair> buf3;
    buf3.reserve( fileSize );

    for ( int i = 0; i < fileSize; ++i )
        buf3.push_back( CharPair( buf1[i], buf2[i] ) );

    std::stable_sort( buf3.begin(), buf3.end() );

    for ( int i = 0; i < fileSize; ++i )
        buf2[i] = buf3[i].second;

    ofstream os( filenameOut );
    os.write( &buf2[0], fileSize );

    return 0;
}
