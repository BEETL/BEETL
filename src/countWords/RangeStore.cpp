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

#include "RangeStore.hh"

#include <cstring>
#include <inttypes.h>
#include <unistd.h>

using namespace std;


//
// RangeStoreRAM - hold BWT intervals in guess where
//

RangeStoreRAM::RangeStoreRAM() : pThis( &r1 ), pNext( &r2 ) {}

RangeStoreRAM::~RangeStoreRAM() {}

void RangeStoreRAM::swap( )
{
    pTemp = pNext;
    pNext = pThis;
    pThis = pTemp;
} // ~swap


void RangeStoreRAM::setPortion( int pileNum, int portionNum )
{
    i_ = ( *pThis )[pileNum][portionNum].begin();
    end_ = ( *pThis )[pileNum][portionNum].end();
} // ~setPortion

bool RangeStoreRAM::getRange( Range &thisRange )
{
    if ( i_ != end_ )
    {
        thisRange = *( i_++ );
        return true;
    }
    else return false;
}  // ~getRange


void RangeStoreRAM::addRange( int pileNum, int portionNum, const string &seq,
                              LetterCountType pos, LetterCountType num, const string &subset )
{
    assert( subset.empty() && "todo" );
    ( *pNext )[pileNum][portionNum].push_back( Range( seq, pos, num ) );
#ifdef DEBUG
    cout << "add range: " << alphabet[pileNum] << " " << alphabet[portionNum]
         << " " << seq << " " << pos << " " << num << endl;
#endif
} // ~addRange

// clear range store ready for next iter
void RangeStoreRAM::clear( void )
{
    for ( int i( 1 ); i < alphabetSize; ++i )
    {
        for ( int j( 1 ); j < alphabetSize; ++j )
        {
            ( *pThis )[i][j].clear();
        } // ~for j
    } // ~for i
} // ~clear

//
// RangeState: this is a proxy class that sits in RangeStoreExternal
// and handles the conversion between compressed representations on file
// and the structs the code expects
//

RangeState::RangeState() : pFile_( NULL )
{
    clear();
}

void RangeState::clear( void )
{
    //    wordLast_.clear();
    // next line makes sure string comparison fails at first char
    //  wordLast_+=notInAlphabet;
    wordLast_[0] = notInAlphabet;
    wordLast_[1] = '\0';
    posLast_ = 0;
    if ( pFile_ != NULL ) fclose( pFile_ );
    pFile_ = NULL;
} // ~clear

void RangeState::addSeq( const string &seq )
{
    //#define COMPRESS_SEQ
#ifdef COMPRESS_SEQ
    if ( wordLast_[0] != notInAlphabet )
    {
        const char *pSeq( seq.c_str() );
        assert( seq.size() < 65536 );
        uint16_t pSeqLen = seq.size();
        char *pLast( wordLast_ );
        while ( *pSeq == *pLast && pSeqLen > 0 )
        {
            ++pSeq;
            --pSeqLen;
            ++pLast;
        }
#ifdef DEBUG
        cout << "FF " << wordLast_ << " " << seq << " " << pSeq << endl;
#endif
        //        fprintf( pFile_, "%s\n", pSeq );
        assert( pSeqLen < 256 );
        fwrite( &pSeqLen, sizeof( uint16_t ), 1, pFile_ );
        fwrite( pSeq, pSeqLen, 1, pFile_ );

        strcpy( pLast, pSeq );
    }
    else
    {
#endif //ifdef COMPRESS_SEQ
        //        fprintf( pFile_, "%s\n", seq.c_str() );
        assert( seq.size() < 65536 );
        uint16_t seqLen = seq.size();
        fwrite( &seqLen, sizeof( uint16_t ), 1, pFile_ );
        fwrite( seq.c_str(), seqLen, 1, pFile_ );

#ifdef COMPRESS_SEQ
        strcpy( wordLast_, seq.c_str() );
    }
#endif //ifdef COMPRESS_SEQ
}

void RangeState::getSeq( string &word )
{
#ifdef COMPRESS_SEQ
    if ( wordLast_[0] != notInAlphabet )
    {
#ifdef DEBUG
        cout << "FF " << word;
#endif
        word = wordLast_;
        unsigned int wordLen = strlen( wordLast_ ); //word.size();

        uint16_t seqLen;
        fread( &seqLen, sizeof( uint16_t ), 1, pFile_ );
        assert( seqLen < 256 );
        assert( seqLen <= wordLen );
        if ( seqLen == 0 || fread( wordLast_ + ( wordLen - seqLen ), seqLen, 1, pFile_ ) != 1 )
        {
            cerr << "Could not get valid data from file. Aborting." << endl;
            exit( -1 );
        }
        //        wordLast_[ seqLen ] = 0;

        unsigned int wordLastLen = seqLen; //strlen( wordLast_ );
        //wordLast_[ wordLastLen - 1] = '\0';
        //--wordLastLen;
        //        for ( unsigned int i( 0 ); i < wordLastLen; i++ )
        //            word[wordLen - wordLastLen + i] = wordLast_[i];
        word = wordLast_;
#ifdef DEBUG
        cout << " -> " << word << endl;
#endif
    }
    else
    {
#endif //ifdef COMPRESS_SEQ
        uint16_t seqLen;
        fread( &seqLen, sizeof( uint16_t ), 1, pFile_ );
        assert( seqLen < 256 );
        if ( fread( wordLast_, seqLen, 1, pFile_ ) != 1 )
        {
            cerr << "Could not get valid data from file. Aborting." << endl;
            exit( -1 );
        }
        wordLast_[ seqLen ] = 0;

        word = wordLast_;
#ifdef COMPRESS_SEQ
    }
#endif //ifdef COMPRESS_SEQ
}

void RangeState::addNum( LetterCountType num )
{
#ifdef DEBUG
    cout << "AN: send " << num << endl;
#endif
    if ( ( num >> 60 ) != 0 )
    {
        cerr << "Overflow in RangeState::addNum. Aborting." << endl;
        exit( -1 );
    }
    num <<= 4;
    LetterCountType num2 = num >> 8;
    unsigned char extraByteCount = 0;
    while ( num2 != 0 )
    {
        ++extraByteCount;
        num2 >>= 8;
    }
    if ( extraByteCount > 15 )
    {
        cerr << "Overflow(2) in RangeState::addNum. Aborting." << endl;
        exit( -1 );
    }
    num |= extraByteCount;

    if ( fwrite( &num, 1 + extraByteCount, 1, pFile_ ) != 1 )
    {
        cerr << "Could not write " << extraByteCount
             << " chars to file. Aborting." << endl;
        exit( -1 );
    }
}

bool RangeState::getNum( LetterCountType &num )
{
    num = 0;
    if ( fread( &num, 1, 1, pFile_ ) != 1 )
        return false;
#ifdef DEBUG
    cout << "AN: got " << num << endl;
#endif
    int count = num & 0xF;
    if ( count > 0 && fread( ( ( char * )&num ) + 1, count, 1, pFile_ ) != 1 )
    {
        cerr << "Could not read " << count << " chars from file " << fileno( pFile_ ) << " pos " << ftell( pFile_ ) << " . Aborting." << endl;
        cerr << "Press any key to exit" << endl;
        for ( ;; )
        {
            sleep( 1 );
            if ( getchar() == '\n' )
                break;
        }
        exit( -1 );
    }
    num >>= 4;

#ifdef DEBUG
    cout << "AN: decoded " << num << endl;
#endif
    return true;
}


//
// RangeStoreExternal
//

RangeStoreExternal::RangeStoreExternal( const string fileStemIn,
                                        const string fileStemOut ) :
    fileStemIn_( fileStemIn ), fileStemOut_( fileStemOut ), stateIn_()
{
    string fileName;
    for ( int i( 1 ); i < alphabetSize; ++i )
    {
        for ( int j( 1 ); j < alphabetSize; ++j )
        {
            stateOut_[i][j].clear();
            // stateOut_[i][j].pFile_=NULL;
            // (*pThis)[i][j].clear();
            getFileName( fileStemIn_, i, j, fileName );
            if ( remove( fileName.c_str() ) == 0 )
            {
#ifdef DEBUG
                cerr << "Removed " << fileName << endl;
#endif
            }
            getFileName( fileStemOut_, i, j, fileName );
            if ( remove( fileName.c_str() ) == 0 )
            {
#ifdef DEBUG
                cerr << "Removed " << fileName << endl;
#endif
            }
        } // ~for j
    } // ~for i
} // ~ctor

RangeStoreExternal::~RangeStoreExternal() {}

void RangeStoreExternal::swap( void )
{
#ifdef DEBUG
    cout << "swap" << endl;
    cout << fileStemIn_ << " " << fileStemOut_ << endl;
#endif
    string temp = fileStemIn_;
    fileStemIn_ = fileStemOut_;
    fileStemOut_ = temp;
#ifdef DEBUG
    cout << fileStemIn_ << " " << fileStemOut_ << endl;
#endif
} // ~swap

void RangeStoreExternal::setPortion( int pileNum, int portionNum )
{
#ifdef DEBUG
    cout << "set portion " << alphabet[pileNum] << alphabet[portionNum] << endl;
#endif
    stateIn_.clear();
    //    if (stateIn_.pFile_!=NULL) fclose(stateIn_.pFile_);
    string fileName;
    getFileName( fileStemIn_, pileNum, portionNum, fileName );
#ifdef DEBUG
    cout << "Made input file name " << fileName << endl;
#endif
    stateIn_.pFile_ = TemporaryRamFile::fopen( fileName.c_str(), "rb" );
#ifdef DEBUG
    if ( stateIn_.pFile_ == NULL )
    {
        cerr << "Warning: no file " << fileName
             << " found, presuming no ranges of interest in this region"
             << endl;
    } // ~if
#endif
}

bool RangeStoreExternal::getRange( Range &thisRange )
{
#ifdef DEBUG
    cout << "get range: " << fileStemIn_ << endl;
#endif
    if ( stateIn_.pFile_ == NULL )
        return false;
    //    else if (fread(&thisRange.pos_,sizeof(LetterCountType),1,
    //     stateIn_.pFile_)!=1)
    else if ( stateIn_.getNum( thisRange.pos_ ) == false )
    {
        fclose( stateIn_.pFile_ );
        stateIn_.pFile_ = NULL;
        return false;

    }
    else
    {
        thisRange.pos_ = ( ( thisRange.pos_ & 1 ) << 63 ) | ( thisRange.pos_ >> 1 );

        //      assert(fread(&thisRange.num_,sizeof(LetterCountType),1,
        //     stateIn_.pFile_)==1);

        if ( !stateIn_.getNum( thisRange.num_ ) )
        {
            cerr << "getNum did not return true. Aborting." << endl;
        }

#ifdef PROPAGATE_PREFIX
        stateIn_.getSeq( thisRange.word_ );
#endif
        //      assert(fgets(buf_,256,stateIn_.pFile_)!=NULL);
        //  buf_[strlen(buf_)-1]='\0';
        //  thisRange.word_=buf_;


#ifdef DEBUG
        cout << "got range: " << fileStemIn_ << " " << thisRange.word_ << " " << thisRange.pos_
             << " " << thisRange.num_ << endl;
#endif
        return true;
    }
} // ~getRange

void RangeStoreExternal::addRange( int pileNum, int portionNum, const string &seq,
                                   LetterCountType pos, LetterCountType num, const string &subset )
{
    if ( !subset.empty() )
    {
        if ( seq.size() >= subset.size() )
        {
            if ( !hasSuffix( seq, subset ) )
                return;
        }
        else
        {
            if ( !hasSuffix( subset, seq ) )
                return;
        }
    }
#ifdef DEBUG
    cout << "set range: " << fileStemIn_ << " " << alphabet[pileNum] << " " << alphabet[portionNum]
         << " " << seq << " " << pos << " " << num << endl;
#endif
    if ( stateOut_[pileNum][portionNum].pFile_ == NULL )
    {
        string fileName;
        getFileName( fileStemOut_, pileNum, portionNum, fileName );
#ifdef DEBUG
        cout << "Made output file name " << fileName << endl;
#endif
        stateOut_[pileNum][portionNum].pFile_ = TemporaryRamFile::fopen( fileName.c_str(), "wb", subset.empty() ? 0 : ( 32 * 1024 * 1024 ) ); // 32 MB * 64 files = 2GB, good for our cluster jobs
    } // ~if
    //    assert(fwrite(&pos,sizeof(LetterCountType),1,
    //   stateOut_[pileNum][portionNum].pFile_)==1);
    //    assert(fwrite(&num,sizeof(LetterCountType),1,
    //   stateOut_[pileNum][portionNum].pFile_)==1);
    stateOut_[pileNum][portionNum].addNum( ( pos << 1 ) | ( ( pos >> 63 ) & 1 ) ); // Move the last bit flag to bit 0 for better downstream compression
    stateOut_[pileNum][portionNum].addNum( num );


#ifdef DEBUG
    cout << "add seq " << seq.c_str() << " " << strlen( seq.c_str() ) << endl;
#endif
    //    fputs(seq.c_str(),fileOut_[pileNum][portionNum]);
    //    fprintf(stateOut_[pileNum][portionNum].pFile_,"%s\n",seq.c_str());
#ifdef PROPAGATE_PREFIX
    stateOut_[pileNum][portionNum].addSeq( seq );
#endif
}

void RangeStoreExternal::clear( void )
{

    string fileName;

    for ( int i( 1 ); i < alphabetSize; ++i )
    {
        for ( int j( 1 ); j < alphabetSize; ++j )
        {
            // if (stateOut_[i][j].pFile_!=NULL)
            //  fclose( stateOut_[i][j].pFile_ );
            // stateOut_[i][j].pFile_=NULL;
            stateOut_[i][j].clear();

            getFileName( fileStemIn_, i, j, fileName );
            if ( remove( fileName.c_str() ) != 0 )
            {
#ifdef DEBUG
                cerr << "Could not remove " << fileName << endl;
#endif
            }

            // (*pThis)[i][j].clear();
        } // ~for j
    } // ~for i
} // ~clear

void RangeStoreExternal::getFileName( const string &stem, const int pile, const int portion,
                                      string &fileName )
{
    fileName = stem;
    fileName += '-';
    fileName += '0';
    fileName += ( char )( 48 + pile );
    fileName += '-';
    fileName += '0';
    fileName += ( char )( 48 + portion );
    if ( pile > 9 || portion > 9 )
    {
        cerr << "Alphabet seems to be larger than 9 chars. Aborting." << endl;
        exit( -1 );
    }

}

