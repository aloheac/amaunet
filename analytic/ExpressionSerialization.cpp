//
// Created by loheac on 2/28/17.
//

#include <fstream>
#include "ExpressionSerialization.h"

using namespace std;

int saveSumToFile( Sum &expr, string filename ) {
    ofstream ofs;

    ofs.open( filename.c_str() );

    if ( not ofs.is_open() ) {
        cout << "***ERROR: Failed to open file '" << filename << "' for writing." << endl;
        return -1;
    }

    boost::archive::text_oarchive oa{ ofs };
    oa << expr;

    cout << "Expression written to file '" << filename << "'." << endl;
    return 0;
}

Sum loadSumFromFile( string filename ) {
    ifstream ifs;

    ifs.open( filename.c_str() );

    if ( not ifs.is_open() ) {
        cout << "***ERROR: Failed to open file '" << filename << "' for reading." << endl;
        return Sum();
    }

    boost::archive::text_iarchive ia{ ifs };
    Sum expr;
    ia >> expr;

    cout << "Expression loaded from file '" << filename << "'." << endl;
    return expr;
}


