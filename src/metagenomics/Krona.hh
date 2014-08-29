void printKronaHeader( ofstream &output )
{
    output << "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">" << endl;
    output << "<!-- adapted from http://krona.sourceforge.net/examples/mg-rast.krona.html -->" << endl;
    output << "<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\" lang=\"en\">" << endl;
    output << "<head>" << endl;
    output << "  <meta charset=\"utf-8\"/>" << endl;
    output << "  <base href=\"https://s3.amazonaws.com/metabeetl/krona/\" target=\"_blank\"/>" << endl;
    output << "  <link rel=\"shortcut icon\" href=\"favicon.ico\"/>" << endl;
    output << "  <script id=\"notfound\">window.onload=function(){document.body.innerHTML=\"Could not get resources from \\\"https://s3.amazonaws.com/metabeetl/krona\\\".\"}</script>" << endl;
    output << "  <script src=\"krona-2.0.js\"></script>" << endl;
    output << " <title>Krona - all</title></head>" << endl;
    output << " <body style=\"padding:0;position:relative\">" << endl;
    output << "  <img id=\"hiddenImage\" src=\"hidden.png\" style=\"display:none\">" << endl;
    output << "  <noscript>Javascript must be enabled to view this page.</noscript>" << endl;
    output << "  <div style=\"display:none\">" << endl;
    output << "  <krona collapse=\"false\" key=\"true\">" << endl;
    output << "   <attributes magnitude=\"magnitude\">" << endl;
    output << "    <attribute display=\"Abundance\">magnitude</attribute>" << endl;
    output << "    <attribute display=\"Rank\">rank</attribute>" << endl;
    output << "    <attribute display=\"Tax id\">taxid</attribute>" << endl;
    output << "   </attributes>" << endl;
}

void printKronaDatasets( ofstream &output, vector<int> &wordMinSize )
{
    output << "   <datasets>" << endl;
    for ( unsigned int s ( 0 ); s < wordMinSize.size(); s++ )
        output << "     <dataset>" << wordMinSize[s] << "</dataset>" << endl;
    output << "   </datasets>" << endl;
}

void printKronaFooter( ofstream &output )
{
    output << "</krona>" << endl;
    output << "</div>" << endl;
    output << "</body></html>" << endl;
}

void printKronaChildren( TAXMAP::iterator &iter, ofstream &output, int level, TAXMAP &taxInfo, unsigned int wordMinSizeCount )
{
    bool dataAvailable = false;
    for ( unsigned int s ( 0 ); s < wordMinSizeCount; s++ )
    {
        int magnitude = iter->second.wordCountPerSize_[s] + iter->second.wordCountPerSizeOfChildren_[s];
        if ( magnitude != 0 ) dataAvailable = true;
    }
    if ( !dataAvailable ) return;

    int id = iter->first;
    string indent( level, ' ' );
    output << indent << "<node name=\"" << iter->second.name_ /* << "(" << id << ")" */ << "\">" << endl;
    output << indent << "<magnitude>";
    for ( unsigned int s ( 0 ); s < wordMinSizeCount; s++ )
    {
        uint64_t magnitude = iter->second.wordCountPerSize_[s] + iter->second.wordCountPerSizeOfChildren_[s];
        output << "<val>" << magnitude << "</val>";
    }
    output << "</magnitude>" << endl;

    unsigned int taxLevel = iter->second.taxLevel_;
    if ( taxLevel < taxLevelSize )
        output << indent << "<rank><val>" << taxLevelNames[taxLevel] << "</val></rank>" << endl;

    output << indent << "<taxid><val>" << id << "</val></taxid>" << endl;

    //    cerr << "krona id " << id << " " << level << endl;

    for ( TAXMAP::iterator iter = taxInfo.begin() ; iter != taxInfo.end(); ++iter )
    {
        if ( iter->second.parentId_ == id && id != 0 )
        {
            printKronaChildren( iter, output, level + 2, taxInfo, wordMinSizeCount );
        }
    }
    output << indent << "</node>" << endl;
}


