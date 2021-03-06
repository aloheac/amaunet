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
	int EXPANSION_ORDER_IN_A = 4;
	int SPLIT_SUMS_BY_LINE = 1;
	int EVAL_BY_PARTS = 0;

	cout << "Loaded parameters:" << endl;
	cout << "\tExpansion order in A:\t\t" << EXPANSION_ORDER_IN_A << endl;
	cout << "\tSplit sums by line:\t\t\t" << SPLIT_SUMS_BY_LINE << endl;
	cout << endl;

	if ( EXPANSION_ORDER_IN_A > 10 ) {
		cout << "***ERROR: The highest perturbation theory order currently implemented is 10th order in A." << endl;
		exit( -1 );
	}
	Sum D1, D2, D3, D4, D5, D6, D7, D8, D9, D10;
	DetM detUp("up");

	initializeStaticReferences();

	if (EXPANSION_ORDER_IN_A >= 2) {
		cout << "Computing first derivative of the determinant..." << endl;
		D1 = detUp.getDerivative();
		D1 = D1.getExpandedExpr();
		D1.reduceTree();
		D1 = distributeAllTraces( D1.copy() );
		cout << "Computing second derivative of the determinant..." << endl;
		D2 = D1.getDerivative();
		D2 = D2.getExpandedExpr();
		D2.reduceTree();
		D2 = distributeAllTraces( D2.copy() );
	}

	if ( EXPANSION_ORDER_IN_A >= 4 ) {
		cout << "Computing third derivative of the determinant..." << endl;
		D3 = D2.getDerivative();
		D3 = D3.getExpandedExpr();
		D3.reduceTree();
		D3 = distributeAllTraces( SymbolicTermPtr( D3.copy() ) );
		cout << "Computing fourth derivative of the determinant..." << endl;
		D4 = D3.getDerivative();
		D4 = D4.getExpandedExpr();
		D4.reduceTree();
		D4 = distributeAllTraces( SymbolicTermPtr( D4.copy() ) );
	}

	if ( EXPANSION_ORDER_IN_A >= 6 ) {
		cout << "Computing fifth derivative of the determinant..." << endl;
		D5 = D4.getDerivative();
		D5 = D5.getExpandedExpr();
		D5.reduceTree();
		D5 = distributeAllTraces( SymbolicTermPtr( D5.copy() ) );
		cout << "Computing sixth derivative of the determinant..." << endl;
		D6 = D5.getDerivative();
		D6 = D6.getExpandedExpr();
		D6.reduceTree();
		D6 = distributeAllTraces( SymbolicTermPtr( D6.copy() ) );
	}

	if ( EXPANSION_ORDER_IN_A >= 8 ) {
		cout << "Computing seventh derivative of the determinant..." << endl;
		D7 = D6.getDerivative();
		D7 = D7.getExpandedExpr();
		D7.reduceTree();
		D7 = distributeAllTraces( SymbolicTermPtr( D7.copy() ) );
		cout << "Computing eighth derivative of the determinant..." << endl;
		D8 = D7.getDerivative();
		D8 = D8.getExpandedExpr();
		D8.reduceTree();
		D8 = distributeAllTraces( SymbolicTermPtr( D8.copy() ) );
	}

	if ( EXPANSION_ORDER_IN_A >= 10 ) {
		cout << "Computing ninth derivative of the determinant..." << endl;
		D9 = D8.getDerivative();
		D9 = D9.getExpandedExpr();
		D9.reduceTree();
		D9 = distributeAllTraces( SymbolicTermPtr( D9.copy() ) );
		cout << "Computing tenth derivative of the determinant..." << endl;
		D10 = D9.getDerivative();
		D10 = D10.getExpandedExpr();
		D10.reduceTree();
		D10 = distributeAllTraces( SymbolicTermPtr( D10.copy() ) );
	}

	cout << "Computing expansion of the partition function..." << endl;

	SumPtr Zup( new Sum() );
	if ( EXPANSION_ORDER_IN_A >= 2 ) {
		Zup->addTerm( detUp.copy() );

		Product P1;
		P1.addTerm( SymbolicTermPtr( new TermA() ) );
		P1.addTerm( D1.copy() );
		Zup->addTerm( P1.copy() );

		Product P2;
		P2.addTerm( SymbolicTermPtr( new CoefficientFraction( 1, 2 ) ) );
		P2.addTerm( SymbolicTermPtr( new TermA() ) );
		P2.addTerm( SymbolicTermPtr( new TermA() ) );
		P2.addTerm( D2.copy() );
		Zup->addTerm( P2.copy() );
	}

	if ( EXPANSION_ORDER_IN_A >= 4 ) {
		Product P3;
		P3.addTerm( SymbolicTermPtr( new CoefficientFraction( 1, 6 ) ) );
		P3.addTerm( SymbolicTermPtr( new TermA() ) );
		P3.addTerm( SymbolicTermPtr( new TermA() ) );
		P3.addTerm( SymbolicTermPtr( new TermA() ) );
		P3.addTerm( D3.copy() );
		Zup->addTerm( P3.copy() );

		Product P4;
		P4.addTerm( SymbolicTermPtr( new CoefficientFraction( 1, 24 ) ) );
		P4.addTerm( SymbolicTermPtr( new TermA() ) );
		P4.addTerm( SymbolicTermPtr( new TermA() ) );
		P4.addTerm( SymbolicTermPtr( new TermA() ) );
		P4.addTerm( SymbolicTermPtr( new TermA() ) );
		P4.addTerm( D4.copy() );
		Zup->addTerm( P4.copy() );
	}

	if ( EXPANSION_ORDER_IN_A >= 6 ) {
		Product P5;
		P5.addTerm( SymbolicTermPtr( new CoefficientFraction( 1, 120 ) ) );
		P5.addTerm( SymbolicTermPtr( new TermA() ) );
		P5.addTerm( SymbolicTermPtr( new TermA() ) );
		P5.addTerm( SymbolicTermPtr( new TermA() ) );
		P5.addTerm( SymbolicTermPtr( new TermA() ) );
		P5.addTerm( SymbolicTermPtr( new TermA() ) );
		P5.addTerm( D5.copy() );
		Zup->addTerm( P5.copy() );

		Product P6;
		P6.addTerm( SymbolicTermPtr( new CoefficientFraction( 1, 270 ) ) );
		P6.addTerm( SymbolicTermPtr( new TermA() ) );
		P6.addTerm( SymbolicTermPtr( new TermA() ) );
		P6.addTerm( SymbolicTermPtr( new TermA() ) );
		P6.addTerm( SymbolicTermPtr( new TermA() ) );
		P6.addTerm( SymbolicTermPtr( new TermA() ) );
		P6.addTerm( SymbolicTermPtr( new TermA() ) );
		P6.addTerm( D6.copy() );
		Zup->addTerm( P6.copy() );
	}

	if ( EXPANSION_ORDER_IN_A >= 8 ) {
		Product P7;
		P7.addTerm( SymbolicTermPtr( new CoefficientFraction( 1, 5040 ) ) );
		P7.addTerm( SymbolicTermPtr( new TermA() ) );
		P7.addTerm( SymbolicTermPtr( new TermA() ) );
		P7.addTerm( SymbolicTermPtr( new TermA() ) );
		P7.addTerm( SymbolicTermPtr( new TermA() ) );
		P7.addTerm( SymbolicTermPtr( new TermA() ) );
		P7.addTerm( SymbolicTermPtr( new TermA() ) );
		P7.addTerm( SymbolicTermPtr( new TermA() ) );
		P7.addTerm( D7.copy() );
		Zup->addTerm( P7.copy() );

		Product P8;
		P8.addTerm( SymbolicTermPtr( new CoefficientFraction( 1, 40320 ) ) );
		P8.addTerm( SymbolicTermPtr( new TermA() ) );
		P8.addTerm( SymbolicTermPtr( new TermA() ) );
		P8.addTerm( SymbolicTermPtr( new TermA() ) );
		P8.addTerm( SymbolicTermPtr( new TermA() ) );
		P8.addTerm( SymbolicTermPtr( new TermA() ) );
		P8.addTerm( SymbolicTermPtr( new TermA() ) );
		P8.addTerm( SymbolicTermPtr( new TermA() ) );
		P8.addTerm( SymbolicTermPtr( new TermA() ) );
		P8.addTerm( D8.copy() );
		Zup->addTerm( P8.copy() );
	}

	if ( EXPANSION_ORDER_IN_A >= 10 ) {
		Product P9;
		P9.addTerm( SymbolicTermPtr( new CoefficientFraction( 1, 362880 ) ) );
		P9.addTerm( SymbolicTermPtr( new TermA() ) );
		P9.addTerm( SymbolicTermPtr( new TermA() ) );
		P9.addTerm( SymbolicTermPtr( new TermA() ) );
		P9.addTerm( SymbolicTermPtr( new TermA() ) );
		P9.addTerm( SymbolicTermPtr( new TermA() ) );
		P9.addTerm( SymbolicTermPtr( new TermA() ) );
		P9.addTerm( SymbolicTermPtr( new TermA() ) );
		P9.addTerm( SymbolicTermPtr( new TermA() ) );
		P9.addTerm( SymbolicTermPtr( new TermA() ) );
		P9.addTerm( D9.copy() );
		Zup->addTerm( P9.copy() );

		Product P10;
		P10.addTerm( SymbolicTermPtr( new CoefficientFraction( 1, 3628800 ) ) );
		P10.addTerm( SymbolicTermPtr( new TermA() ) );
		P10.addTerm( SymbolicTermPtr( new TermA() ) );
		P10.addTerm( SymbolicTermPtr( new TermA() ) );
		P10.addTerm( SymbolicTermPtr( new TermA() ) );
		P10.addTerm( SymbolicTermPtr( new TermA() ) );
		P10.addTerm( SymbolicTermPtr( new TermA() ) );
		P10.addTerm( SymbolicTermPtr( new TermA() ) );
		P10.addTerm( SymbolicTermPtr( new TermA() ) );
		P10.addTerm( SymbolicTermPtr( new TermA() ) );
		P10.addTerm( SymbolicTermPtr( new TermA() ) );
		P10.addTerm( D10.copy() );
		Zup->addTerm( P10.copy() );
	}

	cout << "Expanding flavor determinants..." << endl;
	Zup = static_pointer_cast<Sum>( Zup->getExpandedExpr().copy() );

	SumPtr pathIntegral;
	if ( not EVAL_BY_PARTS ) {
		cout << "Expanding full determinant..." << endl;
		SumPtr Z;
		Z = static_pointer_cast<Sum>(getDualExpansionByParts(static_pointer_cast<Sum>(Zup->copy()),
															 static_pointer_cast<Sum>(Zup->copy())).copy());

		cout << "Reducing expression tree..." << endl;
		Z->reduceTree();

		cout << "Truncating high orders in A of expansion..." << endl;
		Z = static_pointer_cast<Sum>(truncateAOrder(Z, EXPANSION_ORDER_IN_A).copy());

		cout << "Truncating odd orders in A of expansion..." << endl;
		Z = static_pointer_cast<Sum>(truncateOddOrders(Z).copy());

		cout << "Setting all objects in expansion as non-interacting..." << endl;
		Z->setAsNonInteracting();

		cout << "Rewriting expansion in terms of KS formalism..." << endl;
		rewriteSumInKSFormalism(Z);

		cout << "Indexing terms in expansion..." << endl;
		indexExpression(Z);

		cout << "Reducing expression tree..." << endl;
		Z->reduceTree();

		cout << "Computing path integral of expression..." << endl;
		Sum interimPathIntegral;
		interimPathIntegral = pathIntegrateExpression(Z);
		Z.reset();

		cout << "Reducing expression tree..." << endl;
		interimPathIntegral.reduceTree();

		cout << "Performing trivial mathematical simplification..." << endl;
		interimPathIntegral.simplify();

		cout << "Expanding integrated expression..." << endl;
		interimPathIntegral = interimPathIntegral.getExpandedExpr();

		cout << "Reducing expression tree..." << endl;
		interimPathIntegral.reduceTree();

		cout << "Computing analytic Fourier transform of expression..." << endl;
		interimPathIntegral = fourierTransformExpression(interimPathIntegral.copy());

		cout << "Reducing dummy indices of Fourier transformation..." << endl;
		interimPathIntegral.reduceFourierSumIndices();

		cout << "Combining like terms..." << endl;
		interimPathIntegral = combineLikeTerms(interimPathIntegral, 500);

		cout << "Performing trivial mathematical simplification..." << endl;
		interimPathIntegral.simplify();

		pathIntegral = static_pointer_cast<Sum>( interimPathIntegral.copy() );

	} else {
		pathIntegral = expandAndEvaluateExpressionByParts( Zup, Zup, EXPANSION_ORDER_IN_A, 1000 );
	}

	cout << *pathIntegral << endl;
}