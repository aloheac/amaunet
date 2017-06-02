//
// Created by loheac on 2/28/17.
//

#include <fstream>
#include <sstream>
#include "ExpressionSerialization.h"
#include "Multithreading.h"

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
    ofs.close();

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
    ifs.close();

    cout << "Expression loaded from file '" << filename << "'." << endl;
    return expr;
}

int splitSumToFiles( Sum &expr, int blockSize, string saveDir ) {
    cout << ">> Expression contains " << expr.getNumberOfTerms() << " terms to write. " << expr.getNumberOfTerms() / blockSize << " files required." << endl;

    // Check for case where the number of terms in expr is less than or equal blockSize. If this is the case, save the
    // entire sum to one file.
    if ( expr.getNumberOfTerms() <= blockSize ) {
        stringstream ssfilename;
        ssfilename << saveDir << "/EX0.out";
        cout << "***NOTE: Length of expression (" << expr.getNumberOfTerms() << " terms) less than block size. Saving expression to single file." << endl;
        saveSumToFile( expr, ssfilename.str() );
        return 1;
    }

    Sum subsum;
    int fileNo = 0;
    for ( vector<SymbolicTermPtr>::iterator iter = expr.getIteratorBegin(); iter != expr.getIteratorEnd(); ++iter ) {
        if ( subsum.getNumberOfTerms() == blockSize ) {
            // Concatenate filename for the next expression to serialize and output to file.
            stringstream ssfilename;
            ssfilename << saveDir << "/EX" << fileNo << ".out";

            // Save the file.
            int result = saveSumToFile( subsum, ssfilename.str() );

            // Check for an error from the above call.
            if ( result != 0 ) {
                cout << "***ERROR: Failed to save a partial sum." << endl;
                exit( -1 );  // Critical failure -- must terminate calculation.
            }

            // Set subsum to a new instance of Sum.
            subsum = Sum();

            fileNo++;
        }

        subsum.addTerm( (*iter)->copy() );
    }

    return fileNo;
}

Sum loadAndEvaluateSumFromFiles( string saveDir, int numberOfFiles, int EXPANSION_ORDER_IN_A, int POOL_SIZE ) {

    Sum nextPartialSum;
    SumPtr evaluatedPartialSum;
    Sum completeSum;
    for ( int fileNo = 0; fileNo < numberOfFiles; fileNo++ ) {
        stringstream ssfilename;
        ssfilename << saveDir << "/EX" << fileNo << ".out";
        cout << ">> Loading expression from file '" << ssfilename.str() << "'..." << endl;
        nextPartialSum = loadSumFromFile( ssfilename.str() );

        cout << ">> Evaluating partial sum..." << endl;
        evaluatedPartialSum = fullyEvaluateExpressionByParts( static_pointer_cast<Sum>( nextPartialSum.copy() ), EXPANSION_ORDER_IN_A, POOL_SIZE );
        completeSum.addTerm( evaluatedPartialSum->copy() );
        completeSum.reduceTree();

        cout << ">> Combining like terms..." << endl;
        completeSum = combineLikeTerms( completeSum, POOL_SIZE );
    }

    return completeSum;
}