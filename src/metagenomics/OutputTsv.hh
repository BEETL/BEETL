void printTsvHeader( ofstream &output )
{
    output << "#TaxId\tTaxLevel\tCount\tCountIncludingChildren" << endl;
}

void printTsvChildren( TAXMAP::iterator &iter, ofstream &output, int level, TAXMAP &taxInfo, unsigned int wordMinSizeCount )
{
    bool dataAvailable = false;
    for ( unsigned int s ( 0 ); s < wordMinSizeCount; s++ )
    {
        int magnitude = iter->second.wordCountPerSize_[s] + iter->second.wordCountPerSizeOfChildren_[s];
        if ( magnitude != 0 ) dataAvailable = true;
    }
    if ( !dataAvailable ) return;

    int id = iter->first;
    unsigned int taxLevel = iter->second.taxLevel_;

    output << id << '\t' << taxLevel;
    for ( unsigned int s ( 0 ); s < wordMinSizeCount; s++ )
    {
        uint64_t magnitude = iter->second.wordCountPerSize_[s] + iter->second.wordCountPerSizeOfChildren_[s];
        output << '\t' << iter->second.wordCountPerSize_[s] << '\t' << magnitude;
    }
    output << endl;

    for ( TAXMAP::iterator iter = taxInfo.begin() ; iter != taxInfo.end(); ++iter )
    {
        if ( iter->second.parentId_ == id && id != 0 )
        {
            printTsvChildren( iter, output, level + 2, taxInfo, wordMinSizeCount );
        }
    }
}


