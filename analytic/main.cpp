/*
 * Amaunet - High-order Lattice Perturbation Theory
 * 			 for Non-relativistic Quantum Matter
 *
 * High-order Perturbation Theory Analytics
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Executable Entry Point
 *
 * v. 0.1 (alpha)		1 Feb 2016		no-commit
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at Chapel Hill
 */

#include <iostream>
#include <cstdlib>
#include <memory>
#include "PTSymbolicObjects.h"
#include "PathIntegration.h"
#include "Multithreading.h"
#include "ExpressionSerialization.h"

using namespace std;

int main( int argc, char** argv ) {
	// Output introductory information.

	const string WELCOME_HEADER =
			"**********************************************************************\n"
					"    amaunet\n"
					"    --------------------------------------------------------\n"
					"    high-order lattice perturbation theory\n"
					"    for non-relativistic quantum matter\n\n"
					"    High-order Perturbation Theory Analytics\n"
					"**********************************************************************";

	const string PHYSICAL_SYSTEM =
			"Weak-coupling Expansion for Fermionic Contact Interactions";

	const string VERSION_STRING = "v. 0.1 (alpha)";

	const string BUILD_DATE = "1 Feb 2016";

	const string COMMIT_ID = "no-commit";

	cout << WELCOME_HEADER << endl << endl;
	cout << PHYSICAL_SYSTEM << endl;
	cout << VERSION_STRING << "\t\t" << BUILD_DATE << "\t\t" << COMMIT_ID << endl << endl;

	// Load global system parameters.
	int EXPANSION_ORDER_IN_A = 6;
	int SPLIT_SUMS_BY_LINE = 1;
	int EVALUATION_METHOD = 2;
    int POOL_SIZE = 5000;
    int NUM_THREADS = 10;

	cout << "Loaded parameters:" << endl;
	cout << "\tExpansion order in A:\t\t" << EXPANSION_ORDER_IN_A << endl;
	cout << "\tSplit sums by line:\t\t" << SPLIT_SUMS_BY_LINE << endl;
    cout << "\tEvaluation method:\t\t" << EVALUATION_METHOD << endl;
    cout << "\tTerm pool size:\t\t" << POOL_SIZE << endl;
    cout << "\tNumber of threads:\t\t" << NUM_THREADS << endl;
	cout << endl;

	if ( EXPANSION_ORDER_IN_A > 10 ) {
		cout << "***ERROR: The highest perturbation theory order currently implemented is 10th order in A." << endl;
		exit( -1 );
	}

	cout << endl << "Initializing..." << endl << endl;
	initializeStaticReferences();

	Sum Z, Zup, Zdn;
	cout << "Generating series for fermion determinant..." << endl;
	Zup = generateDeterminantExpansion( EXPANSION_ORDER_IN_A, "", true);
    Zdn = generateDeterminantExpansion( EXPANSION_ORDER_IN_A, "", true);

	cout << "Expanding spin-up fermion determinant..." << endl;
	Zup = Zup.getExpandedExpr();

    cout << "Expanding spin-down fermion determinant..." << endl;
    Zdn = Zdn.getExpandedExpr();

	cout << "Reducing expression tree and mathematically simplifying expansion..." << endl;
	Zup.reduceTree();
	Zup.simplify();
    Zdn.reduceTree();
    Zdn.simplify();

    if ( EVALUATION_METHOD == 0 ) {
        cout << "Evaluation method is STANDARD." << endl;
        cout << "Generating product of fermion determinants..." << endl;
        Product fermionDeterminants;
        fermionDeterminants.addTerm( Zup.copy() );
        fermionDeterminants.addTerm( Zdn.copy() );
        Z = Sum( fermionDeterminants.copy() );

        cout << "Expanding product of fermion determinants..." << endl;
        Z = Z.getExpandedExpr();

        cout << "Reducing expression tree..." << endl;
        Z.reduceTree();
        Z.simplify();

        cout << "Truncating high-order terms in expansion..." << endl;
        Z = truncateAOrder( Z.copy(), EXPANSION_ORDER_IN_A );

        cout << "Truncating odd order terms in expansion..." << endl;
        Z = truncateOddOrders( Z.copy() );

        cout << "Sorting traces by order..." << endl;
        Z = sortTracesByOrder( Z );

        SymbolicTermPtr ZPtr = Z.copy();
        cout << "Indexing trace arguments..." << endl;
        indexExpression( ZPtr );

        cout << "Computing path integral of expression..." << endl;
        ZPtr->reduceTree();
        Z = pathIntegrateExpression( ZPtr );

        cout << "Expanding integrated expression..." << endl;
        Z = Z.getExpandedExpr();

        cout << "Reducing expression tree..." << endl;
        Z.reduceTree();

        cout << "Computing symbolic Fourier transform..." << endl;
        Z = fourierTransformExpression( Z.copy() );

        cout << "Reducing dummy indices of Fourier transform..." << endl;
        Z.reduceFourierSumIndices();

        splitSumToFiles( Z, 1000, "./" );

        cout << "Combining like terms..." << endl;
        Z = combineLikeTerms( Z, 1000 );

        cout << Z << endl;

    } else if ( EVALUATION_METHOD == 1 ) {
        cout << "Evaluation method is BY PARTS WITH MULTITHREADING." << endl;
        SumPtr ZPtr = multithreaded_expandAndEvaluateExpressionByParts( static_pointer_cast<Sum>( Zup.copy() ), static_pointer_cast<Sum>( Zdn.copy() ), EXPANSION_ORDER_IN_A, POOL_SIZE, NUM_THREADS );

        cout << ZPtr->to_string() << endl;

    } else if ( EVALUATION_METHOD == 2 ) {
        cout << "Evaluation method is BY PARTS WRITTEN TO FILE WITH MULTITHREADING." << endl;

        SumPtr ZPtr = multithreaded_getDualExpansionByParts( static_pointer_cast<Sum>( Zup.copy() ), static_pointer_cast<Sum>( Zdn.copy() ), NUM_THREADS );

        cout << "Writing split sums to file..." << endl;
        int numFiles;
        numFiles = splitSumToFiles( *ZPtr, 20000, "." );

        Z = loadAndEvaluateSumFromFiles( ".", numFiles, EXPANSION_ORDER_IN_A, POOL_SIZE );
        cout << Z << endl;

    } else {
        cout << "***ERROR: Invalid evaluation method identifier." << endl;
    }
}