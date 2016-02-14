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
#include "PTSymbolicObjects.h"

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
	int EXPANSION_ORDER_IN_A = 8;
	int SPLIT_SUMS_BY_LINE = 1;

	cout << "Loaded parameters:" << endl;
	cout << "\tExpansion order in A:\t\t" << EXPANSION_ORDER_IN_A << endl;
	cout << "\tSplit sums by line:\t\t\t" << SPLIT_SUMS_BY_LINE << endl;
	cout << endl;

	if ( EXPANSION_ORDER_IN_A > 8 ) {
		cout << "***ERROR: The highest perturbation theory order currently implemented is 8th order in A." << endl;
		exit( -1 );
	}
	Sum D1, D2, D3, D4, D5, D6, D7, D8;
	DetM detUp("up");

	if (EXPANSION_ORDER_IN_A >= 2) {
		cout << "Computing first derivative of the determinant..." << endl;
		D1 = detUp.getDerivative();
		D1 = D1.getExpandedExpr();
		D1.reduceTree();
		cout << "Computing second derivative of the determinant..." << endl;
		D2 = D1.getDerivative();
		D2 = D2.getExpandedExpr();
		D2.reduceTree();
	}

	if ( EXPANSION_ORDER_IN_A >= 4 ) {
		cout << "Computing third derivative of the determinant..." << endl;
		D3 = D2.getDerivative();
		D3 = D3.getExpandedExpr();
		D3.reduceTree();
		cout << "Computing fourth derivative of the determinant..." << endl;
		D4 = D3.getDerivative();
		D4 = D4.getExpandedExpr();
		D4.reduceTree();
	}

	if ( EXPANSION_ORDER_IN_A >= 6 ) {
		cout << "Computing fifth derivative of the determinant..." << endl;
		D5 = D4.getDerivative();
		D5 = D5.getExpandedExpr();
		D5.reduceTree();
		cout << "Computing sixth derivative of the determinant..." << endl;
		D6 = D5.getDerivative();
		D6 = D6.getExpandedExpr();
		D6.reduceTree();
	}

	if ( EXPANSION_ORDER_IN_A >= 8 ) {
		cout << "Computing seventh derivative of the determinant..." << endl;
		D7 = D6.getDerivative();
		D7 = D7.getExpandedExpr();
		D7.reduceTree();
		cout << "Computing eighth derivative of the determinant..." << endl;
		D8 = D7.getDerivative();
		D8 = D8.getExpandedExpr();
		D8.reduceTree();
	}
}
