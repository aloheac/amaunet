//
// Created by loheac on 2/28/17.
//

#include <fstream>
#include <sstream>
#include "omp.h"
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

int splitDualExpansionByPartsToFiles( SumPtr exprA, SumPtr exprB, int blockSize, string saveDir ) {
    exprA->reduceTree();
    exprB->reduceTree();

    Sum expandedExpression;
    SymbolicTermPtr exprBCopy = exprB->copy();

    cout << ">> Approximately " << pow( exprA->getNumberOfTerms(), 2 ) / blockSize << " files required to save expansion to disk." << endl;

    int fileNo = 0;
    int termID = 1;
    for ( vector<SymbolicTermPtr>::iterator term = exprA->getIteratorBegin(); term != exprA->getIteratorEnd(); ++term ) {
        cout << ">> Performing dual expression expansion for term " << termID << " of " << exprB->getNumberOfTerms() << "..." << endl;

        Product nextExpansion;
        nextExpansion.addTerm( (*term)->copy() );
        nextExpansion.addTerm( exprBCopy );

        SumPtr expanded = static_pointer_cast<Sum>( nextExpansion.getExpandedExpr().copy() );
        expanded->reduceTree();

        expandedExpression.addTerm( static_pointer_cast<SymbolicTerm>( expanded ) );
        nextExpansion.clear();
        expandedExpression.reduceTree();  // Required so an accurate count of the total number of terms is obtained.

        termID++;

        if ( expandedExpression.getNumberOfTerms() >= blockSize ) {
            cout << ">> Dumping expanded expression to file..." << endl;

            // Concatenate filename for the next expression to serialize and output to file.
            stringstream ssfilename;
            ssfilename << saveDir << "/EX" << fileNo << ".out";

            // Save the file.
            int result = saveSumToFile( expandedExpression, ssfilename.str() );

            // Check for an error from the above call.
            if ( result != 0 ) {
                cout << "***ERROR: Failed to save a partial sum." << endl;
                exit( -1 );  // Critical failure -- must terminate calculation.
            }

            // Set subsum to a new instance of Sum.
            expandedExpression = Sum();

            fileNo++;
        }
    }

    // Check for case where the fileNo is still zero, indicating that the expression was too small to be split across
    // multiple files. If so, save the entire expression to one file now.
    if ( fileNo == 0 ) {
        stringstream ssfilename;
        ssfilename << saveDir << "/EX0.out";
        cout << "***NOTE: Length of expression (" << expandedExpression.getNumberOfTerms() << " terms) less than block size. Saving expression to single file." << endl;
        saveSumToFile( expandedExpression, ssfilename.str() );
        return 1;
    }

    cout << ">> Dual expansion complete. " << fileNo << " files saved." << endl;
    return fileNo;
}

int multithreaded_splitDualExpansionByPartsToFiles( SumPtr exprA, SumPtr exprB, int blockSize, string saveDir, int NUM_THREADS ) {
    exprA->reduceTree();
    exprB->reduceTree();

    Sum expandedExpression;
    SymbolicTermPtr exprBCopy = exprB->copy();
    int fileNo = 0;
    const int numOfBlocks = (int)ceil( (float)exprA -> getNumberOfTerms() / (float)blockSize );

    omp_set_num_threads( NUM_THREADS );

    cout << ">> Note " << numOfBlocks << " files required to save expansion to disk." << endl;



    for ( int block = 0; block < numOfBlocks; block++ ) {
        cout << ">> Expanding block " << block + 1 << " of " << numOfBlocks << "..." << endl;

        int numTermsComplete = 0;
        SumPtr parallelParts[ blockSize ];

#pragma omp parallel for shared( exprA, exprBCopy, parallelParts )
        for (int term = 0; term < blockSize; term++) {
#pragma omp critical(printcout)
            {
                cout << ">> Performing expression expansion for term " << term << " of "
                     << blockSize << " (" << numTermsComplete << " terms complete) in block " << block + 1 << " of "
                     << numOfBlocks << "..." << endl;
                numTermsComplete++;
            }

            // Since we are working in the units of blocks, verify that the total term index (block * blockSize + term)
            // does not exceed the number of terms in the expression.
            if ( block * blockSize + term < exprA->getNumberOfTerms() ) {
                Product nextExpansion;
                nextExpansion.addTerm(exprA->getTerm(block * blockSize + term)->copy());
                nextExpansion.addTerm(exprBCopy);

                SumPtr expanded = static_pointer_cast<Sum>(nextExpansion.getExpandedExpr().copy());
                expanded->reduceTree();

                parallelParts[term] = expanded;
                nextExpansion.clear();
            }
        }

        cout << ">> Evaluation of block complete. Dumping expanded expression to file..." << endl;

        // Reduce parallel results.
        Sum reducedExpression;
        for ( int term = 0; term < blockSize; term++ ) {
            // Again, verify that we are not going past the end of the total expression.
            if ( block * blockSize + term < exprA->getNumberOfTerms() ) {
                reducedExpression.addTerm( parallelParts[ term ] );
            }
        }

        // Concatenate filename for the next expression to serialize and output to file.
        stringstream ssfilename;
        ssfilename << saveDir << "/EX" << fileNo << ".out";

        // Save the file.
        int result = saveSumToFile(reducedExpression, ssfilename.str());

        // Check for an error from the above call.
        if (result != 0) {
            cout << "***ERROR: Failed to save a partial sum." << endl;
            exit(-1);  // Critical failure -- must terminate calculation.
        }

        // Set subsum to a new instance of Sum.
        expandedExpression = Sum();

        fileNo++;
    }

    // Check for case where the fileNo is still zero, indicating that the expression was too small to be split across
    // multiple files. If so, save the entire expression to one file now.
    if ( fileNo == 0 ) {
        stringstream ssfilename;
        ssfilename << saveDir << "/EX0.out";
        cout << "***NOTE: Length of expression (" << expandedExpression.getNumberOfTerms() << " terms) less than block size. Saving expression to single file." << endl;
        saveSumToFile( expandedExpression, ssfilename.str() );
        return 1;
    }

    cout << ">> Dual expansion complete. " << fileNo << " files saved." << endl;
    return fileNo;
}