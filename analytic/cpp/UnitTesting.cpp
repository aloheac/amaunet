/* ***********************************************************************
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * Numerical CUDA Perturbation Theory Fourier Transform Evaluation
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Primary Unit Testing
 *
 * v. 0.1		09 Feb 2016
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at Chapel Hill
 * ***********************************************************************
 */

#include <iostream>
#include <sstream>
#include <string>
#include "PTSymbolicObjects.h"
#include "PathIntegration.h"

using namespace std;

class UnitTest {
public:
	UnitTest( string name, string (*tst)(), string result ) {
		cout << "Running unit test '" << name << "'..." << endl;
		string ret = tst();

		if ( ret == result ) {
			passedTests++;
			cout << ">> Unit test PASSED." << endl;
		} else {
			failedTests++;
			cout << ">> Unit test FAILED. Result:" << endl;
			cout << ret << endl;
			cout << "   Correct result: " << endl;
			cout << result << endl;
 		}
	}

	static int passedTests;

	static int failedTests;
};

int UnitTest::passedTests = 0;
int UnitTest::failedTests = 0;

/*
 * Unit test functions.
 */

string i01() {
	stringstream ss;
	GenericTestTerm A = GenericTestTerm( 0, 0 );
	ss << A;
	return ss.str();
}

string i02() {
	stringstream ss;
	GenericTestTerm A = GenericTestTerm( 0, 0 );
	ss << A.getDerivative();
	return ss.str();
}

string i03() {
	stringstream ss;
	GenericTestTerm A = GenericTestTerm( 1, 0 );
	GenericTestTermPtr B( static_pointer_cast<GenericTestTerm>( A.copy() ) );
	ss << *B << " " << A;
	//delete B;
	return ss.str();
}

string i04() {
	stringstream ss;
	GenericTestTerm A = GenericTestTerm( 1, 0 );
	GenericTestTermPtr B( static_pointer_cast<GenericTestTerm>( A.copy() ) );
	Sum C = A.getDerivative();
	ss << *B << " " << C;
	//delete B;
	return ss.str();
}

string i05() {
	stringstream ss;
	GenericTestTerm A = GenericTestTerm( 0, 0 );
	ss << A.getTermID();
	return ss.str();
}

string i06() {
	stringstream ss;
	initializeStaticReferences();
	ss << "1: " << Amaunet::SINE_PATH_INTEGRALS[ 1 ] << "    2: " << Amaunet::SINE_PATH_INTEGRALS[ 2 ]
		<< "    3: " << Amaunet::SINE_PATH_INTEGRALS[ 3 ] << "    4: " << Amaunet::SINE_PATH_INTEGRALS[ 4 ];
	return ss.str();
}

string A01() {
	stringstream ss;
	SymbolicTerm A = SymbolicTerm();
	ss << A;
	return ss.str();
}

string A02() {
	stringstream ss;
	SymbolicTerm A = SymbolicTerm();
	ss << A.getDerivative();
	return ss.str();
}

string A03() {
	stringstream ss;
	SymbolicTerm A = SymbolicTerm();
	ss << A.getTermID();
	return ss.str();
}

string A04() {
	stringstream ss;
	SymbolicTerm A;
	SymbolicTermPtr B = A.copy();
	A.setAsNonInteracting();
	ss << *B << "    " << B->isTermInteracting();
	return ss.str();
}

string B01() {
	stringstream ss;
	TermA A = TermA();
	ss << A;
	return ss.str();
}

string B02() {
	stringstream ss;
	TermA A = TermA();
	TermA B = TermA();
	ss << (A == B);
	return ss.str();
}

string B03() {
	stringstream ss;
	TermA A = TermA();
	TermAPtr B( static_pointer_cast<TermA>( A.copy() ) );
	ss << *B;
	//delete B;
	return ss.str();
}

string B04() {
	stringstream ss;
	TermA A = TermA();
	ss << A.getTermID();
	return ss.str();
}

string B05() {
	stringstream ss;
	TermA A = TermA();
	TermAPtr B( static_pointer_cast<TermA>( A.copy() ) );
	ss << B->getTermID();
	//delete B;
	return ss.str();
}

string C01() {
	stringstream ss;
	MatrixM A = MatrixM();
	ss << A;
	return ss.str();
}

string C02() {
	stringstream ss;
	MatrixM A = MatrixM();
	A.setAsNonInteracting();
	ss << A;
	return ss.str();
}

string C03() {
	stringstream ss;
	MatrixM A = MatrixM( "up" );
	ss << A;
	return ss.str();
}

string C04() {
	stringstream ss;
	MatrixM A = MatrixM( "up" );
	A.setAsNonInteracting();
	ss << A;
	return ss.str();
}

string C05() {
	stringstream ss;
	MatrixM A = MatrixM( "up" );
	ss << A.getDerivative();
	return ss.str();
}

string C06() {
	stringstream ss;
	MatrixM A = MatrixM( "up" );
	ss << A.getDerivative().getDerivative();
	return ss.str();
}

string C07() {
	stringstream ss;
	MatrixM A = MatrixM( "up" );
	A.setAsNonInteracting();
	ss << A.getDerivative();
	return ss.str();
}

string C08() {
	stringstream ss;
	MatrixM A = MatrixM( "up" );
	MatrixMPtr B( static_pointer_cast<MatrixM>( A.copy() ) );
	B->setAsNonInteracting();
	ss << A << " " << *B << " " << A.isTermInteracting() << " " << B->isTermInteracting();
	//delete B;
	return ss.str();
}

string C09() {
	stringstream ss;
	MatrixM A = MatrixM( "up" );
	MatrixMPtr B( static_pointer_cast<MatrixM>( A.copy() ) );
	Sum C = B->getDerivative();
	ss << C;
	//delete B;
	return ss.str();
}

string D01() {
	stringstream ss;
	CoefficientFloat A = CoefficientFloat( 0.0 );
	ss << A;
	return ss.str();
}

string D02() {
	stringstream ss;
	CoefficientFloat A = CoefficientFloat( 3.0 );
	ss << A << "    " << A.getDerivative();
	return ss.str();
}

string E01() {
	stringstream ss;
	MatrixB A = MatrixB();
	ss << A;
	return ss.str();
}

string E02() {
	stringstream ss;
	MatrixB A = MatrixB();
	A.setAsNonInteracting();
	ss << A;
	return ss.str();
}

string E03() {
	stringstream ss;
	MatrixB A = MatrixB( "up" );
	ss << A;
	return ss.str();
}

string E04() {
	stringstream ss;
	MatrixB A = MatrixB( "up" );
	A.setAsNonInteracting();
	ss << A;
	return ss.str();
}

string E05() {
	stringstream ss;
	MatrixB A = MatrixB( "up" );
	ss << A.getDerivative();
	return ss.str();
}

string E06() {
	stringstream ss;
	MatrixB A = MatrixB( "up" );
	Sum B = A.getDerivative().getDerivative();
	B.reduceTree();
	ss << B;
	return ss.str();
}

string E07() {
	stringstream ss;
	MatrixB A = MatrixB( "up" );
	Sum B = A.getDerivative().getDerivative().getDerivative();
	B.reduceTree();
	ss << B;
	return ss.str();
}

string F01() {
	stringstream ss;
	MatrixK A;
	ss << A;
	return ss.str();
}

string F02() {
	stringstream ss;
	MatrixK A( "up" );
	ss << A;
	return ss.str();
}

string F03() {
	stringstream ss;
	MatrixK A( "up" );
	A.setAsNonInteracting();
	ss << A;
	return ss.str();
}

string F04() {
	stringstream ss;
	MatrixK A( "up" );
	A.fourierTransform();
	ss << A;
	return ss.str();
}

string G01() {
	stringstream ss;
	MatrixS A;
	ss << A;
	return ss.str();
}

string H01() {
	stringstream ss;
	MatrixM m( "up" );
	DetM det( m );
	ss << det;
	return ss.str();
}

string H02() {
	stringstream ss;
	MatrixM m( "up" );
	DetM det( m );
	SymbolicTermPtr cpy( det.copy() );
	det.setAsNonInteracting();
	ss << det << "    " << *cpy;
	//delete cpy;
	return ss.str();
}

string H03() {
	stringstream ss;
	MatrixM m( "up" );
	DetM det( m );
	Sum derivative = det.getDerivative();
	ss << derivative;
	return ss.str();
}

string H04() {
	stringstream ss;
	MatrixM m( "up" );
	DetM det( m );
	Sum derivative = det.getDerivative().getDerivative();
	derivative.reduceTree();
	ss << derivative;
	return ss.str();
}

string H05() {
	stringstream ss;
	MatrixM m( "up" );
	DetM det( m );
	Sum derivative = det.getDerivative().getDerivative().getDerivative();
	derivative.reduceTree();
	ss << derivative;
	return ss.str();
}

string I01() {
	stringstream ss;
	CoefficientFraction A( 3, 7 );
	ss << A;
	return ss.str();
}

string I02() {
	stringstream ss;
	CoefficientFraction A( 5, 9 );
	ss << A << "    " << A.getDerivative();
	return ss.str();
}

string I03() {
	stringstream ss;
	Sum A;
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( CoefficientFractionPtr( new CoefficientFraction( 3, 7 ) ) );
	ss << A << "    ";
	Sum B = A.getDerivative();
	ss << B << "    ";
	B.simplify();
	ss << B;
	return ss.str();
}

string I04() {
	stringstream ss;
	CoefficientFraction A( 4, 8 );
	ss << A << "    ";
	A.reduce();
	ss << A;
	return ss.str();
}

string I05() {
	stringstream ss;
	CoefficientFraction A( 84, 174 );
	ss << A << "    ";
	A.reduce();
	ss << A;
	return ss.str();
}

string I06() {
	stringstream ss;
	CoefficientFraction A( 3, 7 );
	ss << A << "    ";
	A.reduce();
	ss << A;
	return ss.str();
}

string I07() {
	stringstream ss;
	CoefficientFraction A( 3, 7 );
	CoefficientFraction B( 4, 9 );
	ss << A << "    " << B << "    " << A * B;
	return ss.str();
}

string I08() {
	stringstream ss;
	CoefficientFraction A( 1, 2 );
	CoefficientFraction B( 1, 4 );
	ss << A << "    " << B << "    " << A + B;
	return ss.str();
}

string I09() {
	stringstream ss;
	CoefficientFraction A( 5, 3 );
	CoefficientFraction B( 8, 3 );
	ss << A << "    " << B << "    " << A + B;
	return ss.str();
}

string I10() {
	stringstream ss;
	CoefficientFraction A( 7, 2 );
	CoefficientFraction B( 8, 4 );
	ss << A << "    " << B << "    " << A + B;
	return ss.str();
}

string I11() {
	stringstream ss;
	CoefficientFraction A( 5.4, 7.9 );
	ss << A << "    ";
	A.reduce();
	ss << A;
	return ss.str();
}

string I12() {
	stringstream ss;
	CoefficientFraction A( 1, 3 );
	CoefficientFloat B( 5 );
	ss << A << "    " << B << "    " << A * B;
	return ss.str();
}

string I13() {
	stringstream ss;
	CoefficientFraction A( 5, 8 );
	CoefficientFloat B( 2 );
	CoefficientFraction C = A * B;
	C.reduce();
	ss << A << "    " << B << "    " << C;
	return ss.str();
}

string I14() {
	stringstream ss;
	CoefficientFraction A( 5, 8 );
	CoefficientFloat B( -2 );
	CoefficientFraction C = A * B;
	C.reduce();
	ss << A << "    " << B << "    " << C;
	return ss.str();
}

string I15() {
	stringstream ss;
	CoefficientFraction A( 1, 2 );
	CoefficientFloat B( -1 );
	CoefficientFraction C = A * B;
	C.reduce();
	ss << A << "    " << B << "    " << C << "    " << C * B;
	return ss.str();
}

string I16() {
	stringstream ss;
	CoefficientFraction A( 1, 2 );
	CoefficientFloat B( 1 );
	CoefficientFraction C = A + B;
	C.reduce();
	ss << A << "    " << B << "    " << C;
	return ss.str();
}

string I17() {
	stringstream ss;
	CoefficientFraction A( 1, 2 );
	CoefficientFloat B( -1 );
	CoefficientFraction C = A + B;
	C.reduce();
	ss << A << "    " << B << "    " << C;
	return ss.str();
}

string I18() {
	stringstream ss;
	CoefficientFraction A( 2, 8 );
	CoefficientFloat B( -1 );
	CoefficientFraction C = A + B;
	C.reduce();
	ss << A << "    " << B << "    " << C;
	return ss.str();
}


string J01() {
	stringstream ss;
	Sum A = Sum( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	ss << A;
	return ss.str();
}

string J02() {
	stringstream ss;
	Sum A = Sum();
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	ss << A;
	return ss.str();
}

string J03() {
	stringstream ss;
	Sum A = Sum();
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	A = A.getDerivative();
	ss << A;
	return ss.str();
}

string J04() {
	stringstream ss;
	Sum A = Sum();
	A.addTerm( MatrixMPtr( new MatrixM( "a" ) ) );
	A.addTerm( MatrixMPtr( new MatrixM( "b" ) ) );
	A.addTerm( MatrixMPtr( new MatrixM( "c" ) ) );
	SumPtr B( dynamic_pointer_cast<Sum>( A.copy() ) );
	B->setAsNonInteracting();
	ss << A << " " << *B;
	//delete B;
	return ss.str();
}

string J05() {
	stringstream ss;
	Sum A = Sum();
	A.addTerm( MatrixMPtr( new MatrixM( "a" ) ) );
	A.addTerm( MatrixMPtr( new MatrixM( "b" ) ) );
	A.addTerm( CoefficientFloatPtr( new CoefficientFloat( 0.0 ) ) );
	A.addTerm( MatrixMPtr( new MatrixM( "c" ) ) );
	ss << A << " ";
	A.simplify();
	ss << A;
	return ss.str();
}

string J06() {
	stringstream ss;
	Sum A = Sum();
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	Sum B = Sum();
	B.addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	B.addTerm( GenericTestTermPtr( new GenericTestTerm(3,0) ) );
	Sum C = Sum();
	C.addTerm( SymbolicTermPtr( A.copy() ) );
	C.addTerm( SymbolicTermPtr( B.copy() ) );
	ss << C << "   " << C.getNumberOfTerms() << "    ";
	C.reduceTree();
	ss << C << "   " << C.getNumberOfTerms();
	return ss.str();
}

string J07() {
	stringstream ss;
	Sum A;
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	ProductPtr B( new Product() );
	B->addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	ProductPtr C( new Product() );
	C->addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	C->addTerm( GenericTestTermPtr( new GenericTestTerm(3,0) ) );
	B->addTerm( C );
	A.addTerm( B );
	ss << A << "    " << A.getNumberOfTerms() << "    ";
	A.reduceTree();
	ss << A << "    " << A.getNumberOfTerms();
	return ss.str();
}

string J08() {
	stringstream ss;
	Sum A;
	Product B;
	B.addTerm( SymbolicTermPtr( new CoefficientFloat( -1.0 ) ) );
	Product C;
	C.addTerm( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	C.addTerm( SymbolicTermPtr( new GenericTestTerm(1,0) ) );
	B.addTerm( C.copy() );
	B.addTerm( SymbolicTermPtr( new GenericTestTerm(2,0) ) );
	A.addTerm( B.copy() );
	ss << A << "    " << A.getNumberOfTerms() << "    ";
	A.reduceTree();
	ss << A << "    " << A.getNumberOfTerms();
	return ss.str();
}

string J09() {
	stringstream ss;
	Sum A;
	Product B;
	Product C;
	C.addTerm( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	C.addTerm( SymbolicTermPtr( new GenericTestTerm(1,0) ) );
	B.addTerm( C.copy() );
	A.addTerm( B.copy() );
	ss << A << "    " << A.getNumberOfTerms() << "    ";
	A.reduceTree();
	ss << A << "    " << A.getNumberOfTerms();
	return ss.str();
}

string J10() {
	stringstream ss;
	Sum A;
	Product B;
	Product C;
	Product D;
	D.addTerm( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	D.addTerm( SymbolicTermPtr( new GenericTestTerm(1,0) ) );
	C.addTerm( D.copy() );
	B.addTerm( C.copy() );
	B.addTerm( SymbolicTermPtr( new GenericTestTerm(2,0) ) );
	A.addTerm( B.copy() );
	ss << A << "    " << A.getNumberOfTerms() << "    ";
	A.reduceTree();
	ss << A << "    " << A.getNumberOfTerms();
	return ss.str();
}

string J11() {
	stringstream ss;
	Sum A;
	Product B;
	B.addTerm( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	B.addTerm( SymbolicTermPtr( new GenericTestTerm(1,0) ) );
	Product C;
	C.addTerm( SymbolicTermPtr( new CoefficientFloat( 3.0 ) ) );
	C.addTerm( B.copy() );
	Product D;
	D.addTerm( SymbolicTermPtr( new CoefficientFloat( 5.0 ) ) );
	D.addTerm( C.copy() );
	Product E;
	E.addTerm( D.copy() );
	A.addTerm( E.copy() );
	ss << A << "    " << A.getNumberOfTerms() << "    ";
	A.reduceTree();
	ss << A << "    " << A.getNumberOfTerms();
	return ss.str();
}

string J12() {
	stringstream ss;
	Sum A;
	A.addTerm( CoefficientFloatPtr( new CoefficientFloat( 0.0 ) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( CoefficientFloatPtr( new CoefficientFloat( 0.0 ) ) );
	A.simplify();
	ss << A;
	return ss.str();
}

string K01() {
	stringstream ss;
	Product A = Product( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	ss << A;
	return ss.str();
}

string K02() {
	stringstream ss;
	Product A = Product();
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(3,0) ) );
	ss << A;
	return ss.str();
}

string K03() {
	stringstream ss;
	Product A = Product();
	A.addTerm( MatrixMPtr( new MatrixM() ) );
	A.addTerm( MatrixMPtr( new MatrixM() ) );
	A.addTerm( MatrixBPtr( new MatrixB() ) );
	A.setAsNonInteracting();
	ss << A;
	return ss.str();
}

string K04() {
	stringstream ss;
	Product A = Product();
	A.addTerm( MatrixMPtr( new MatrixM() ) );
	A.addTerm( SumPtr( new Sum() ) );
	A.addTerm( MatrixBPtr( new MatrixB() ) );
	ss << A.containsSum();
	return ss.str();
}

string K05() {
	stringstream ss;
	Product A = Product();
	A.addTerm( MatrixMPtr( new MatrixM() ) );
	A.addTerm( TermAPtr( new TermA() ) );
	A.addTerm( MatrixBPtr( new MatrixB() ) );
	ss << A.containsSum();
	return ss.str();
}

string K06() {
	stringstream ss;
	Product A = Product();
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	Sum B = A.getDerivative();
	ss << B;
	return ss.str();
}

string K07() {
	stringstream ss;
	Product A = Product();
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	A.addTerm( CoefficientFloatPtr( new CoefficientFloat( 0.0 ) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	A.simplify();
	ss << A;
	return ss.str();
}

string K08() {
	stringstream ss;
	Product A = Product();
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	A.addTerm( CoefficientFloatPtr( new CoefficientFloat( 1.0 ) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	A.simplify();
	ss << A;
	return ss.str();
}

string K09() {
	stringstream ss;
	Product A = Product();
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	A.addTerm( CoefficientFloatPtr( new CoefficientFloat( 1.0 ) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	A.addTerm( CoefficientFloatPtr( new CoefficientFloat( 0.0 ) ) );
	A.simplify();
	ss << A;
	return ss.str();
}

string K10() {
	stringstream ss;
	Sum A;
	A.addTerm( ProductPtr( new Product( GenericTestTermPtr( new GenericTestTerm(0,0) ) ) ) );
	ss << A;
	return ss.str();
}

string K11() {
	stringstream ss;
	Product A;
	A.addTerm( GenericTestTermPtr(  new GenericTestTerm(0,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	ProductPtr B = ProductPtr( new Product() );
	B->addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	B->addTerm( GenericTestTermPtr( new GenericTestTerm(3,0) ) );
	A.addTerm( B );
	ss << A << "    " << A.getNumberOfTerms() << "    ";
	A.reduceTree();
	ss << A << "    " << A.getNumberOfTerms();
	return ss.str();
}

string K12() {
	stringstream ss;
	Product A;
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	ss << A << "    ";
	Sum B = A.getExpandedExpr();
	ss << B;
	return ss.str();
}

string K13() {
	stringstream ss;
	Sum A;
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	Sum B;
	B.addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	B.addTerm( GenericTestTermPtr( new GenericTestTerm(3,0) ) );
	Product C;
	C.addTerm( A.copy() );
	C.addTerm( B.copy() );
	ss << C << "    ";
	Sum D = C.getExpandedExpr();
	ss << D;
	return ss.str();
}

string K14() {
	stringstream ss;
	Sum A;
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	Sum B;
	B.addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	B.addTerm( GenericTestTermPtr( new GenericTestTerm(3,0) ) );
	B.addTerm( GenericTestTermPtr( new GenericTestTerm(4,0) ) );
	Product C;
	C.addTerm( GenericTestTermPtr( new GenericTestTerm(5,0) ) );
	C.addTerm( A.copy() );
	C.addTerm( B.copy() );
	ss << C << "    ";
	Sum D = C.getExpandedExpr();
	ss << D;
	return ss.str();
}

string K15() {
	stringstream ss;
	Sum A;
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	Sum B;
	B.addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	B.addTerm( GenericTestTermPtr( new GenericTestTerm(3,0) ) );
	B.addTerm( GenericTestTermPtr( new GenericTestTerm(4,0) ) );
	Product C;
	C.addTerm( A.copy() );
	C.addTerm( GenericTestTermPtr( new GenericTestTerm(5,0) ) );
	C.addTerm( B.copy() );
	ss << C << "    ";
	Sum D = C.getExpandedExpr();
	D.reduceTree();
	ss << D;
	return ss.str();
}

string K16() {
	stringstream ss;
	Sum A;
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(1,0) ) );
	Sum B;
	B.addTerm( GenericTestTermPtr( new GenericTestTerm(2,0) ) );
	B.addTerm( GenericTestTermPtr( new GenericTestTerm(3,0) ) );
	Product C;
	C.addTerm( GenericTestTermPtr( new GenericTestTerm(4,0) ) );
	C.addTerm( A.copy() );
	Product D;
	D.addTerm( GenericTestTermPtr( new GenericTestTerm(5,0) ) );
	D.addTerm( B.copy() );
	Product E;
	E.addTerm( C.copy() );
	E.addTerm( D.copy() );
	ss << E << "    ";
	E.reduceTree();
	ss << E << "    ";
	Sum F = E.getExpandedExpr();
	F.reduceTree();
	ss << F;
	return ss.str();
}

string K17() {
	stringstream ss;
	Product A;
	A.addTerm( CoefficientFloatPtr( new CoefficientFloat( 1.0 ) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( CoefficientFloatPtr( new CoefficientFloat( 1.0 ) ) );
	A.simplify();
	ss << A;
	return ss.str();
}

string K18() {
	stringstream ss;
	Product A;
	A.addTerm( GenericTestTermPtr( new GenericTestTerm( 0, 0 ) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm( 1, 0 ) ) );
	A.addTerm( GenericTestTermPtr( new GenericTestTerm( 2, 0 ) ) );
	ss << A << "    ";
	A.zero();
	ss << A;
	return ss.str();
}

string L01() {
	stringstream ss;
	SymbolicTermPtr term = GenericTestTermPtr( new GenericTestTerm(0,0) );
	Trace A( term );
	ss << A;
	return ss.str();
}

string L02() {
	stringstream ss;
	SymbolicTermPtr term = GenericTestTermPtr( new GenericTestTerm(0,0) );
	Trace A( term );
	Sum B = A.getDerivative();
	ss << B;
	return ss.str();
}

string L03() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new MatrixB( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	SymbolicTermPtr B = A.copy();
	Trace C( B );
	ss << C << "    ";
	C.rewriteInKSFormalism();
	ss << C;
	return ss.str();
}

string L04() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new MatrixB( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixB( "dn" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "dn" ) ) );
	SymbolicTermPtr B = A.copy();
	Trace C( B );
	ss << C << "    ";
	C.rewriteInKSFormalism();
	ss << C;
	return ss.str();
}

string M01() {
	stringstream ss;
	Delta d( 0, 1 );
	ss << d;
	return ss.str();
}

string M02() {
	stringstream ss;
	Delta d( 0, 1, true );
	ss << d;
	return ss.str();
}

string N01() {
	stringstream ss;
	vector<IndexContraction> A;
	A.push_back( IndexContraction( 0, 1 ) );
	A.push_back( IndexContraction( 2, 3 ) );
	A.push_back( IndexContraction( 4, 5 ) );
	FourierSum B( A, 3 );
	ss << B;
	return ss.str();
}

string N02() {
	stringstream ss;
	vector<IndexContraction> A;
	A.push_back( IndexContraction( 0, 1 ) );
	A.push_back( IndexContraction( 2, 3 ) );
	A.push_back( IndexContraction( 4, 5 ) );
	FourierSum B( A, 3 );
	vector<IndexContraction> C;
	C.push_back( IndexContraction( 0, 1 ) );
	C.push_back( IndexContraction( 2, 3 ) );
	C.push_back( IndexContraction( 4, 5 ) );
	FourierSum D( C, 3 );
	ss << B << "    " << D << "    " << ( B == D );
	return ss.str();
}

string N03() {
	stringstream ss;
	vector<IndexContraction> A;
	A.push_back( IndexContraction( 0, 1 ) );
	A.push_back( IndexContraction( 2, 3 ) );
	A.push_back( IndexContraction( 4, 5 ) );
	FourierSum B( A, 3 );
	vector<IndexContraction> C;
	C.push_back( IndexContraction( 0, 1 ) );
	C.push_back( IndexContraction( 2, 3 ) );
	C.push_back( IndexContraction( 7, 9 ) );
	FourierSum D( C, 3 );
	ss << B << "    " << D << "    " << ( B == D );
	return ss.str();
}

string N04() {
	stringstream ss;
	vector<IndexContraction> A;
	A.push_back( IndexContraction( 0, 1 ) );
	A.push_back( IndexContraction( 2, 3 ) );
	A.push_back( IndexContraction( 4, 5 ) );
	FourierSum B( A, 3 );
	vector<IndexContraction> C;
	C.push_back( IndexContraction( 0, 1 ) );
	C.push_back( IndexContraction( 4, 5 ) );
	C.push_back( IndexContraction( 2, 3 ) );
	FourierSum D( C, 3 );
	ss << B << "    " << D << "    " << ( B == D );
	return ss.str();
}

string N05() {
	stringstream ss;
	vector<IndexContraction> A;
	A.push_back( IndexContraction( 0, 1 ) );
	A.push_back( IndexContraction( 2, 3 ) );
	A.push_back( IndexContraction( 4, 5 ) );
	FourierSum B( A, 3 );
	vector<IndexContraction> C;
	C.push_back( IndexContraction( 0, 1 ) );
	C.push_back( IndexContraction( 4, 5 ) );
	C.push_back( IndexContraction( 3, 2 ) );
	FourierSum D( C, 3 );
	ss << B << "    " << D << "    " << ( B == D );
	return ss.str();
}

string N06() {
	stringstream ss;
	vector<IndexContraction> A;
	A.push_back( IndexContraction( 0, 1 ) );
	A.push_back( IndexContraction( 2, 3 ) );
	A.push_back( IndexContraction( 4, 5 ) );
	FourierSum B( A, 3 );
	vector<IndexContraction> C;
	C.push_back( IndexContraction( 0, 1 ) );
	C.push_back( IndexContraction( 2, 3 ) );
	C.push_back( IndexContraction( 4, 5 ) );
	C.push_back( IndexContraction( 6, 7 ) );
	FourierSum D( C, 4 );
	ss << B << "    " << D << "    " << ( B == D );
	return ss.str();
}

string N07() {
	stringstream ss;
	vector<IndexContraction> A;
	A.push_back( IndexContraction( 0, 1 ) );
	A.push_back( IndexContraction( 2, 3 ) );
	A.push_back( IndexContraction( 4, 5 ) );
	SymbolicTermPtr B = FourierSumPtr( new FourierSum( A, 3 ) );
	ss << *B;
	return ss.str();
}

string N08() {
	stringstream ss;
	vector<IndexContraction> A;
	A.push_back( IndexContraction( 0, 1 ) );
	A.push_back( IndexContraction( 3, 3 ) );
	A.push_back( IndexContraction( 0, 0 ) );
	FourierSum B = FourierSum( A, 3 );
	ss << B << "    ";
	B.reduceDummyIndices();
	ss << B;
	return ss.str();
}

string O01() {
	stringstream ss;
	Product A = Product();
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	SymbolicTermPtr B( A.copy() );
	unpackTrivialExpression( B );
	ss << *B << " " << B->getTermID();
	//delete B;
	return ss.str();
}

string O02() {
	stringstream ss;
	Product A = Product();
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(1,0) ) );
	SymbolicTermPtr B( A.copy() );
	unpackTrivialExpression( B );
	ss << *B << " " << B->getTermID();
	//delete B;
	return ss.str();
}

string O03() {
	stringstream ss;
	Sum A = Sum();
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	SymbolicTermPtr B( A.copy() );
	unpackTrivialExpression( B );
	ss << *B << " " << B->getTermID();
	//delete B;
	return ss.str();
}

string O04() {
	stringstream ss;
	Sum A;
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(1,0) ) );
	SymbolicTermPtr B( A.copy() );
	unpackTrivialExpression( B );
	ss << *B << " " << B->getTermID();
	//delete B;
	return ss.str();
}

string O05() {
	stringstream ss;
	Sum A;
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(1,0) ) );
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(2,0) ) );
	Product B;
	B.addTerm( A.copy() );
	ss << B << "    " << B.getTermID() << "    ";
	SymbolicTermPtr C( B.copy() );
	unpackTrivialExpression( C );
	ss << *C << "    " << C->getTermID();
	return ss.str();
}

string O06() {
	stringstream ss;
	Sum A;
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	Product B;
	B.addTerm( A.copy() );
	Sum C;
	C.addTerm( B.copy() );
	ss << C << "    ";
	SymbolicTermPtr D( C.copy() );
	unpackTrivialExpression( D );
	ss << *D;
	//delete D;
	return ss.str();
}

string P01() {
	stringstream ss;
	Sum A;
	Trace B( A.copy() );
	ss << isZeroTrace( B.copy() );
	return ss.str();
}

string P02() {
	stringstream ss;
	Sum A;
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	Trace B( A.copy() );
	ss << isZeroTrace( B.copy() );
	return ss.str();
}

string P03() {
	stringstream ss;
	Product A;
	Trace B( A.copy() );
	ss << isZeroTrace( B.copy() );
	return ss.str();
}

string P04() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	Trace B( A.copy() );
	ss << isZeroTrace( B.copy() );
	return ss.str();
}

string P05() {
	stringstream ss;
	GenericTestTermPtr A = GenericTestTermPtr( new GenericTestTerm(0,0) );
	ss << isZeroTrace( A );
	//delete A;
	return ss.str();
}

string Q01() {
	stringstream ss;
	SymbolicTermPtr A( new Trace( SymbolicTermPtr( new MatrixM( "up" ) ) ) );
	ss << distributeTrace( A );
	return ss.str();
}

string Q02() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new CoefficientFloat( -1.0 ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	SymbolicTermPtr B( new Trace( A.copy() ) );
	ss << distributeTrace( B );
	return ss.str();
}

string Q03() {
	stringstream ss;
	Sum A;
	A.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "dn" ) ) );
	Product B;
	B.addTerm( SymbolicTermPtr( new CoefficientFloat( -1.0 ) ) );
	B.addTerm( A.copy() );
	SymbolicTermPtr C( new Trace( B.copy() ) );
	ss << distributeTrace( C );
	return ss.str();
}

string Q04() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new CoefficientFloat( 3.0 ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	Product B;
	B.addTerm( SymbolicTermPtr( new CoefficientFloat( 5.0 ) ) );
	B.addTerm( SymbolicTermPtr( new MatrixM( "dn" ) ) );
	Sum C;
	C.addTerm( A.copy() );
	C.addTerm( B.copy() );
	Product D;
	D.addTerm( SymbolicTermPtr( new CoefficientFloat( -1.0 ) ) );
	D.addTerm( C.copy() );
	SymbolicTermPtr E( new Trace( D.copy() ) );
	ss << *E << "    " << distributeTrace( E );
	return ss.str();
}

string R01() {
	stringstream ss;
	Sum A;
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(0,0) ) );
	A.addTerm( SymbolicTermPtr( new GenericTestTerm(1,0) ) );
	Product B;
	B.addTerm( SymbolicTermPtr( new TermA() ) );
	B.addTerm( SymbolicTermPtr( new Trace( A.copy() ) ) );
	SymbolicTermPtr C = B.copy();
	ss << *C << "    " << distributeAllTraces( C );
	return ss.str();
}

string R02() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new CoefficientFloat( 3.0 ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	Product B;
	B.addTerm( SymbolicTermPtr( new CoefficientFloat( 5.0 ) ) );
	B.addTerm( SymbolicTermPtr( new MatrixM( "dn" ) ) );
	Trace C( A.copy() );
	Trace D( B.copy() );
	Product E;
	E.addTerm( SymbolicTermPtr( new CoefficientFloat( -1.0 ) ) );
	E.addTerm( C.copy() );
	E.addTerm( D.copy() );
	SymbolicTermPtr F = E.copy();
	Sum G = distributeAllTraces( F );
	ss << *F << "    " << G;
	return ss.str();
}

string R03() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new CoefficientFloat( 3.0 ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	Product B;
	B.addTerm( SymbolicTermPtr( new CoefficientFloat( 5.0 ) ) );
	B.addTerm( SymbolicTermPtr( new MatrixM( "dn" ) ) );
	Trace C( A.copy() );
	Trace D( B.copy() );
	Sum E;
	E.addTerm( C.copy() );
	E.addTerm( D.copy() );
	SymbolicTermPtr F = E.copy();
	Sum G = distributeAllTraces( F );
	ss << *F << "    " << G;
	return ss.str();
}

string R04() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new CoefficientFloat( 3.0 ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	Product B;
	B.addTerm( SymbolicTermPtr( new CoefficientFloat( 5.0 ) ) );
	B.addTerm( SymbolicTermPtr( new MatrixM( "dn" ) ) );
	Trace C( A.copy() );
	Trace D( B.copy() );
	Product E;
	E.addTerm( SymbolicTermPtr( new TermA() ) );
	E.addTerm( C.copy() );
	Product F;
	F.addTerm( SymbolicTermPtr( new TermA() ) );
	F.addTerm( D.copy() );
	Sum G;
	G.addTerm( E.copy() );
	G.addTerm( F.copy() );
	SymbolicTermPtr H = G.copy();
	Sum I = distributeAllTraces( H );
	ss << *H << "    " << I;
	return ss.str();
}

string R05() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new TermA() ) );
	Product B;
	B.addTerm( SymbolicTermPtr( new CoefficientFloat( 3.0 ) ) );
	B.addTerm( SymbolicTermPtr( new MatrixB( "up" ) ) );
	B.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new Trace( B.copy() ) ) );
	Product C;
	C.addTerm( SymbolicTermPtr( new TermA() ) );
	Product D;
	D.addTerm( SymbolicTermPtr( new CoefficientFloat( 5.0 ) ) );
	D.addTerm( SymbolicTermPtr( new MatrixB( "md" ) ) );
	D.addTerm( SymbolicTermPtr( new MatrixM( "md" ) ) );
	C.addTerm( SymbolicTermPtr( new Trace( D.copy() ) ) );
	Product E;
	E.addTerm( SymbolicTermPtr( new TermA() ) );
	Product F;
	F.addTerm( SymbolicTermPtr( new CoefficientFloat( 7.0 ) ) );
	F.addTerm( SymbolicTermPtr( new MatrixB( "dn" ) ) );
	F.addTerm( SymbolicTermPtr( new MatrixM( "dn" ) ) );
	E.addTerm( SymbolicTermPtr( new Trace( F.copy() ) ) );
	Sum G;
	G.addTerm( A.copy() );
	G.addTerm( C.copy() );
	G.addTerm( E.copy() );
	SymbolicTermPtr H = G.copy();
	ss << *H << "    " << distributeAllTraces( H );
	return ss.str();
}

string S01() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new MatrixB( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixB( "dn" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "dn" ) ) );
	SymbolicTermPtr B = A.copy();
	Trace C( B );
	Sum D;
	D.addTerm( C.copy() );
	ss << D << "    ";
	SymbolicTermPtr E = D.copy();
	rewriteSumInKSFormalism( E );
	ss << *E;
	return ss.str();
}

string S02() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new MatrixB( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixB( "dn" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "dn" ) ) );
	SymbolicTermPtr B = A.copy();
	Trace C( B );
	Sum D;
	D.addTerm( C.copy() );
	Product E;
	E.addTerm( SymbolicTermPtr( new MatrixB( "a" ) ) );
	E.addTerm( SymbolicTermPtr( new MatrixM( "a" ) ) );
	E.addTerm( SymbolicTermPtr( new MatrixB( "a" ) ) );
	E.addTerm( SymbolicTermPtr( new MatrixM( "a" ) ) );
	SymbolicTermPtr F = E.copy();
	Trace G( F );
	D.addTerm( G.copy() );
	Product H;
	H.addTerm( SymbolicTermPtr( new MatrixB( "c" ) ) );
	H.addTerm( SymbolicTermPtr( new MatrixM( "d" ) ) );
	H.addTerm( SymbolicTermPtr( new MatrixB( "e" ) ) );
	H.addTerm( SymbolicTermPtr( new MatrixM( "f" ) ) );
	H.addTerm( SymbolicTermPtr( new MatrixB( "g" ) ) );
	H.addTerm( SymbolicTermPtr( new MatrixM( "h" ) ) );
	SymbolicTermPtr I = H.copy();
	Trace J( I );
	D.addTerm( J.copy() );
	ss << D << "    ";
	SymbolicTermPtr K = D.copy();
	rewriteSumInKSFormalism( K );
	ss << *K;
	return ss.str();
}

string T01() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new MatrixB( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixB( "dn" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "dn" ) ) );
	Trace B( A.copy() );
	Product E;
	E.addTerm( B.copy() );
	Sum C( E.copy() );
	ss << C << "    ";
	SymbolicTermPtr D = C.copy();
	rewriteSumInKSFormalism( D );
	indexExpression( D );
	ss << *D;
	return ss.str();
}

string T02() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new MatrixB( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixB( "dn" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "dn" ) ) );
	Trace B( A.copy() );
	Product E;
	E.addTerm( B.copy() );

	Product F;
	F.addTerm( SymbolicTermPtr( new MatrixB( "up" ) ) );
	F.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	F.addTerm( SymbolicTermPtr( new MatrixB( "dn" ) ) );
	F.addTerm( SymbolicTermPtr( new MatrixM( "dn" ) ) );
	F.addTerm( SymbolicTermPtr( new MatrixB( "md" ) ) );
	F.addTerm( SymbolicTermPtr( new MatrixM( "md" ) ) );
	Trace G( F.copy() );
	Product H;
	H.addTerm( G.copy() );

	Sum C;
	C.addTerm( E.copy() );
	C.addTerm( H.copy() );

	ss << C << "    ";
	SymbolicTermPtr D = C.copy();
	rewriteSumInKSFormalism( D );
	indexExpression( D );
	ss << *D;
	return ss.str();
}

string T03() {
	stringstream ss;
	Product A;
	A.addTerm( SymbolicTermPtr( new MatrixB( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixB( "dn" ) ) );
	A.addTerm( SymbolicTermPtr( new MatrixM( "dn" ) ) );
	Trace B( A.copy() );
	Product C;
	C.addTerm( B.copy() );
	Product D;
	D.addTerm( SymbolicTermPtr( new MatrixB( "up" ) ) );
	D.addTerm( SymbolicTermPtr( new MatrixM( "up" ) ) );
	Trace E( D.copy() );
	C.addTerm( E.copy() );
	Sum F;
	F.addTerm( C.copy() );
	ss << F << "    ";
	SymbolicTermPtr G = F.copy();
	rewriteSumInKSFormalism( G );
	indexExpression( G );
	ss << *G;
	return ss.str();
}

string U01() {
	stringstream ss;
	IndexContraction A;
	A.i = 3;
	A.j = 7;
	ss << A.i << "    " << A.j;
	return ss.str();
}

string V01() {
	stringstream ss;
	DeltaContractionSet A;
	IndexContraction B( 1, 2 );
	IndexContraction C( 3, 5 );
	IndexContraction D( 7, 11 );
	A.addContraction( B );
	A.addContraction( C );
	A.addContraction( D );
	ss << A.to_string();
	return ss.str();
}

string V02() {
	stringstream ss;
	DeltaContractionSet A;
	IndexContraction B( 1, 2 );
	IndexContraction C( 3, 5 );
	IndexContraction D( 7, 11 );
	A.addContraction( B );
	A.addContraction( C );
	A.addContraction( D );
	ss << A;
	return ss.str();
}

string V03() {
	stringstream ss;
	DeltaContractionSet A;
	IndexContraction B( 1, 2 );
	IndexContraction C( 5, 3 );
	IndexContraction D( 7, 11 );
	IndexContraction E( 5, 5 );
	IndexContraction F( 6, 5 );
	A.addContraction( B );
	A.addContraction( C );
	A.addContraction( D );
	A.addContraction( E );
	A.addContraction( F );
	ss << A << "    ";
	A.orderContractionIndices();
	ss << A;
	return ss.str();
}

string V04() {
	stringstream ss;
	DeltaContractionSet A;
	IndexContraction B( 1, 2 );
	IndexContraction C( 5, 3 );
	IndexContraction D( 4, 7 );
	IndexContraction E( 9, 12 );
	IndexContraction F( 8, 8 );
	A.addContraction( B );
	A.addContraction( C );
	A.addContraction( D );
	A.addContraction( E );
	A.addContraction( F );
	ss << A << "    ";
	A.sortContractions();
	ss << A;
	return ss.str();
}


string V05() {
	stringstream ss;
	DeltaContractionSet A;
	IndexContraction B( 1, 2 );
	IndexContraction C( 5, 3 );
	IndexContraction D( 1, 0 );
	IndexContraction E( 8, 12 );
	IndexContraction F( 8, 8 );
	A.addContraction( B );
	A.addContraction( C );
	A.addContraction( D );
	A.addContraction( E );
	A.addContraction( F );
	ss << A << "    ";
	A.sortContractions();
	ss << A;
	return ss.str();
}

string X01() {
	stringstream ss;
	DeltaSignature A;
	DeltaContractionSet B;
	IndexContraction C( 1, 2 );
	IndexContraction D( 3, 5 );
	IndexContraction E( 7, 11 );
	B.addContraction( C );
	B.addContraction( D );
	B.addContraction( E );
	A.addContractionSet( B );
	DeltaContractionSet F;
	IndexContraction G( 13, 17 );
	IndexContraction H( 19, 31 );
	F.addContraction( G );
	F.addContraction( H );
	A.addContractionSet( F );
	DeltaContractionSet I;
	A.addContractionSet( I );
	ss << A.to_string();
	return ss.str();
}

string X02() {
	stringstream ss;
	DeltaSignature A;
	DeltaContractionSet B;
	IndexContraction C( 1, 2 );
	IndexContraction D( 3, 5 );
	IndexContraction E( 7, 11 );
	B.addContraction( C );
	B.addContraction( D );
	B.addContraction( E );
	A.addContractionSet( B );
	DeltaContractionSet F;
	IndexContraction G( 13, 17 );
	IndexContraction H( 19, 31 );
	F.addContraction( G );
	F.addContraction( H );
	A.addContractionSet( F );
	DeltaContractionSet I;
	A.addContractionSet( I );
	ss << A;
	return ss.str();
}

string Y01() {
	stringstream ss;
	vector<int> A;
	A.push_back( 2 );
	TotalSignature B = getDeltaSignature( A );
	ss << B.deltas << "    " << B.deltaBars;
	return ss.str();
}

string Y02() {
	stringstream ss;
	vector<int> A;
	A.push_back( 4 );
	TotalSignature B = getDeltaSignature( A );
	ss << B.deltas << "    " << B.deltaBars;
	return ss.str();
}

string Y03() {
	stringstream ss;
	vector<int> A;
	A.push_back( 2 );
	A.push_back( 2 );
	TotalSignature B = getDeltaSignature( A );
	ss << B.deltas << "    " << B.deltaBars;
	return ss.str();
}

string Y04() {
	stringstream ss;
	vector<int> A;
	A.push_back( 6 );
	TotalSignature B = getDeltaSignature( A );
	ss << B.deltas << "    " << B.deltaBars;
	return ss.str();
}

string Y05() {
	stringstream ss;
	vector<int> A;
	A.push_back( 4 );
	A.push_back( 2 );
	TotalSignature B = getDeltaSignature( A );
	ss << B.deltas << "    " << B.deltaBars;
	return ss.str();
}

string Y06() {
	stringstream ss;
	vector<int> A;
	A.push_back( 2 );
	A.push_back( 2 );
	A.push_back( 2 );
	TotalSignature B = getDeltaSignature( A );
	ss << B.deltas << "    " << B.deltaBars;
	return ss.str();
}

string Y07() {
	stringstream ss;
	vector<int> A;
	A.push_back( 8 );
	TotalSignature B = getDeltaSignature( A );
	ss << B.deltas << "    " << B.deltaBars;
	return ss.str();
}

string Y08() {
	stringstream ss;
	vector<int> A;
	A.push_back( 6 );
	A.push_back( 2 );
	TotalSignature B = getDeltaSignature( A );
	ss << B.deltas << "    " << B.deltaBars;
	return ss.str();
}

string Y09() {
	stringstream ss;
	vector<int> A;
	A.push_back( 4 );
	A.push_back( 4 );
	TotalSignature B = getDeltaSignature( A );
	ss << B.deltas << "    " << B.deltaBars;
	return ss.str();
}

string Y10() {
	stringstream ss;
	vector<int> A;
	A.push_back( 4 );
	A.push_back( 2 );
	A.push_back( 2 );
	TotalSignature B = getDeltaSignature( A );
	ss << B.deltas << "    " << B.deltaBars;
	return ss.str();
}

string Y11() {
	stringstream ss;
	vector<int> A;
	A.push_back( 2 );
	A.push_back( 2 );
	A.push_back( 2 );
	A.push_back( 2 );
	TotalSignature B = getDeltaSignature( A );
	ss << B.deltas << "    " << B.deltaBars;
	return ss.str();
}

string Z01() {
	stringstream ss;
	vector<int> A;
	A.push_back( 1 );
	ss << combinations( A, 1 );
	return ss.str();
}

string Z02() {
	stringstream ss;
	vector<int> A;
	A.push_back( 1 );
	A.push_back( 2 );
	A.push_back( 3 );
	ss << combinations( A, 1 );
	return ss.str();
}

string Z03() {
	stringstream ss;
	vector<int> A;
	A.push_back( 1 );
	A.push_back( 2 );
	A.push_back( 3 );
	A.push_back( 4 );
	A.push_back( 5 );
	ss << combinations( A, 2 ) << "    " << combinations( A, 2 ).size();
	return ss.str();
}

string Z04() {
	stringstream ss;
	vector<int> A;
	A.push_back( 1 );
	A.push_back( 2 );
	A.push_back( 3 );
	A.push_back( 4 );
	A.push_back( 5 );
	ss << combinations( A, 3 ) << "    " << combinations( A, 3 ).size();
	return ss.str();
}

string Z05() {
	stringstream ss;
	vector<int> A;
	A.push_back( 1 );
	A.push_back( 2 );
	A.push_back( 3 );
	A.push_back( 4 );
	A.push_back( 5 );
	ss << combinations( A, 4 ) << "    " << combinations( A, 4 ).size();
	return ss.str();
}

string Z06() {
	stringstream ss;
	vector<int> A;
	A.push_back( 1 );
	A.push_back( 2 );
	A.push_back( 3 );
	A.push_back( 4 );
	A.push_back( 5 );
	ss << combinations( A, 5 ) << "    " << combinations( A, 5 ).size();
	return ss.str();
}

string Z07() {
	stringstream ss;
	vector<int> A;
	A.push_back( 1 );
	A.push_back( 2 );
	A.push_back( 3 );
	A.push_back( 4 );
	A.push_back( 5 );
	A.push_back( 6 );
	A.push_back( 7 );
	ss << combinations( A, 1 ) << "    " << combinations( A, 1 ).size();
	return ss.str();
}

string Z08() {
	stringstream ss;
	vector<int> A;
	A.push_back( 1 );
	A.push_back( 2 );
	A.push_back( 3 );
	A.push_back( 4 );
	A.push_back( 5 );
	A.push_back( 6 );
	A.push_back( 7 );
	ss << combinations( A, 2 ) << "    " << combinations( A, 2 ).size();
	return ss.str();
}

string Z09() {
	stringstream ss;
	vector<int> A;
	A.push_back( 1 );
	A.push_back( 2 );
	A.push_back( 3 );
	A.push_back( 4 );
	A.push_back( 5 );
	A.push_back( 6 );
	A.push_back( 7 );
	ss << combinations( A, 3 ) << "    " << combinations( A, 3 ).size();
	return ss.str();
}

string Z10() {
	stringstream ss;
	vector<int> A;
	A.push_back( 1 );
	A.push_back( 2 );
	A.push_back( 3 );
	A.push_back( 4 );
	A.push_back( 5 );
	A.push_back( 6 );
	A.push_back( 7 );
	ss << combinations( A, 4 ) << "    " << combinations( A, 4 ).size();
	return ss.str();
}

string Z11() {
	stringstream ss;
	vector<int> A;
	A.push_back( 1 );
	A.push_back( 2 );
	A.push_back( 3 );
	A.push_back( 4 );
	A.push_back( 5 );
	A.push_back( 6 );
	A.push_back( 7 );
	ss << combinations( A, 5 ) << "    " << combinations( A, 5 ).size();
	return ss.str();
}

string Z12() {
	stringstream ss;
	vector<int> A;
	ss << combinations( A, 0 );
	return ss.str();
}

string Z13() {
	stringstream ss;
	vector<int> A;
	A.push_back( 1 );
	A.push_back( 2 );
	A.push_back( 3 );
	ss << combinations( A, 0 );
	return ss.str();
}

string AA01() {
	stringstream ss;
	vector<int> A;
	A.push_back( 2 );
	ss << getIndexPermutations( A );
	return ss.str();
}

string AA02() {
	stringstream ss;
	vector<int> A;
	A.push_back( 4 );
	ss << getIndexPermutations( A );
	return ss.str();
}

string AA03() {
	stringstream ss;
	vector<int> A;
	A.push_back( 2 );
	A.push_back( 2 );
	ss << getIndexPermutations( A );
	return ss.str();
}

string AA04() {
	stringstream ss;
	vector<int> A;
	A.push_back( 6 );
	ss << getIndexPermutations( A );
	return ss.str();
}

string AA05() {
	stringstream ss;
	vector<int> A;
	A.push_back( 4 );
	A.push_back( 2 );
	ss << getIndexPermutations( A );
	return ss.str();
}

string AA06() {
	stringstream ss;
	vector<int> A;
	A.push_back( 2 );
	A.push_back( 2 );
	A.push_back( 2 );
	ss << getIndexPermutations( A );
	return ss.str();
}

string AB01() {
	stringstream ss;
	vector<int> A;
	A.push_back( 2 );
	TotalSignature B = getDeltaSignature( A );
	vector< vector<int> > C = getIndexPermutations( A );
	ss << generateSignaturePermutations( C, B );
	return ss.str();
}

string AB02() {
	stringstream ss;
	vector<int> A;
	A.push_back( 4 );
	TotalSignature B = getDeltaSignature( A );
	vector< vector<int> > C = getIndexPermutations( A );
	ss << generateSignaturePermutations( C, B );
	return ss.str();
}

string AB03() {
	stringstream ss;
	vector<int> A;
	A.push_back( 2 );
	A.push_back( 2 );
	TotalSignature B = getDeltaSignature( A );
	vector< vector<int> > C = getIndexPermutations( A );
	ss << generateSignaturePermutations( C, B );
	return ss.str();
}

string AB04() {
	stringstream ss;
	vector<int> A;
	A.push_back( 6 );
	TotalSignature B = getDeltaSignature( A );
	vector< vector<int> > C = getIndexPermutations( A );
	ss << generateSignaturePermutations( C, B );
	return ss.str();
}

string AB05() {
	stringstream ss;
	vector<int> A;
	A.push_back( 4 );
	A.push_back( 2 );
	TotalSignature B = getDeltaSignature( A );
	vector< vector<int> > C = getIndexPermutations( A );
	ss << generateSignaturePermutations( C, B );
	return ss.str();
}

string AB06() {
	stringstream ss;
	vector<int> A;
	A.push_back( 2 );
	A.push_back( 2 );
	A.push_back( 2 );
	TotalSignature B = getDeltaSignature( A );
	vector< vector<int> > C = getIndexPermutations( A );
	ss << generateSignaturePermutations( C, B );
	return ss.str();
}

string AC01() {
	stringstream ss;
	vector<int> A;
	A.push_back( 2 );
	TotalSignature B = getDeltaSignature( A );
	TotalSignature C;
	DeltaContractionSet D;
	D.addContraction( IndexContraction( 2, 3 ) );
	C.deltas = D;
	DeltaContractionSet E;
	C.deltaBars = E;
	ss << B.deltas << "    " << C.deltas << "    " << B.areSignaturesDegenerate( C );
	return ss.str();
}

string AC02() {
	stringstream ss;
	vector<int> A;
	A.push_back( 2 );
	TotalSignature B = getDeltaSignature( A );
	TotalSignature C;
	DeltaContractionSet D;
	D.addContraction( IndexContraction( 2, 1 ) );
	C.deltas = D;
	DeltaContractionSet E;
	C.deltaBars = E;
	ss << B.deltas << "    " << C.deltas << "    " << B.areSignaturesDegenerate( C );
	return ss.str();
}

string AC03() {
	stringstream ss;
	TotalSignature A;
	TotalSignature B;
	DeltaContractionSet C;
	DeltaContractionSet D;
	DeltaContractionSet E;
	DeltaContractionSet F;
	C.addContraction( IndexContraction( 2, 3 ) );
	C.addContraction( IndexContraction( 1, 5 ) );
	D.addContraction( IndexContraction( 1, 2 ) );
	D.addContraction( IndexContraction( 5, 4 ) );
	A.deltas = C;
	B.deltas = D;
	A.deltaBars = E;
	A.deltaBars = F;
	ss << A.deltas << "    " << B.deltas << "    " << A.areSignaturesDegenerate( B );
	return ss.str();
}

string AC04() {
	stringstream ss;
	TotalSignature A;
	TotalSignature B;
	DeltaContractionSet C;
	DeltaContractionSet D;
	DeltaContractionSet E;
	DeltaContractionSet F;
	C.addContraction( IndexContraction( 2, 3 ) );
	C.addContraction( IndexContraction( 1, 5 ) );
	D.addContraction( IndexContraction( 1, 2 ) );
	D.addContraction( IndexContraction( 5, 4 ) );
	A.deltas = C;
	B.deltas = D;
	A.deltaBars = E;
	A.deltaBars = F;
	ss << A.deltas << "    " << B.deltas << "    " << B.areSignaturesDegenerate( A );
	return ss.str();
}

string AC05() {
	stringstream ss;
	TotalSignature A;
	TotalSignature B;
	DeltaContractionSet C;
	DeltaContractionSet D;
	DeltaContractionSet E;
	DeltaContractionSet F;
	C.addContraction( IndexContraction( 1, 2 ) );
	C.addContraction( IndexContraction( 3, 4 ) );
	D.addContraction( IndexContraction( 3, 4 ) );
	D.addContraction( IndexContraction( 1, 2 ) );
	A.deltas = C;
	B.deltas = D;
	A.deltaBars = E;
	A.deltaBars = F;
	ss << A.deltas << "    " << B.deltas << "    " << B.areSignaturesDegenerate( A );
	return ss.str();
}

string AD01() {
	stringstream ss;
	ss << calculateAllContractions( 2 );
	return ss.str();
}

string AD02() {
	stringstream ss;
	ss << calculateAllContractions( 4 );
	return ss.str();
}

string AD03() {
	stringstream ss;
	ss << calculateAllContractions( 6 );
	return ss.str();
}

string AD04() {
	stringstream ss;
	ss << calculateAllContractions( 8 );
	return ss.str();
}

string AD05() {
	stringstream ss;
	ss << calculateAllContractions( 10 );
	return ss.str();
}

string AE01() {
	stringstream ss;
	ss << generateCoordinateSpacePathIntegral( 2 );
	return ss.str();
}

string AE02() {
	stringstream ss;
	ss << generateCoordinateSpacePathIntegral( 4 );
	return ss.str();
}

string AE03() {
	stringstream ss;
	ss << generateCoordinateSpacePathIntegral( 6 );
	return ss.str();
}

string AF01() {
	stringstream ss;
	Sum A;
	Product B;
	MatrixK C;
	C.setIndices( 0, 1 );
	MatrixS D;
	D.setIndices( 1, 0 );
	B.addTerm( C.copy() );
	B.addTerm( D.copy() );
	C.setIndices( 2, 3 );
	D.setIndices( 3, 2 );
	B.addTerm( C.copy() );
	B.addTerm( D.copy() );
	A.addTerm( B.copy() );
	ss << A << "    " << pathIntegrateExpression( A.copy() );
	return ss.str();
}

string AF02() {
	stringstream ss;
	Sum A;
	Product B;
	MatrixK C;
	C.setIndices( 0, 1 );
	MatrixS D;
	D.setIndices( 1, 0 );
	B.addTerm( C.copy() );
	B.addTerm( D.copy() );
	C.setIndices( 2, 3 );
	D.setIndices( 3, 2 );
	B.addTerm( C.copy() );
	B.addTerm( D.copy() );
	C.setIndices( 4, 5 );
	D.setIndices( 5, 4 );
	B.addTerm( C.copy() );
	B.addTerm( D.copy() );
	C.setIndices( 6, 7 );
	D.setIndices( 7, 6 );
	B.addTerm( C.copy() );
	B.addTerm( D.copy() );
	A.addTerm( B.copy() );
	ss << A << "    " << pathIntegrateExpression( A.copy() );
	return ss.str();
}

string AG01() {
	stringstream ss;
	Sum A;
	Product B;
	B.addTerm( SymbolicTermPtr( new TermA() ) );
	A.addTerm( B.copy() );
	B.addTerm( SymbolicTermPtr( new TermA() ) );
	A.addTerm( B.copy() );
	ss << A << "    " << truncateAOrder( A.copy(), 1 );
	return ss.str();
}

string AH01() {
	stringstream ss;
	Sum A;
	Product B;
	B.addTerm( SymbolicTermPtr( new TermA() ) );
	A.addTerm( B.copy() );
	B.addTerm( SymbolicTermPtr( new TermA() ) );
	A.addTerm( B.copy() );
	B.addTerm( SymbolicTermPtr( new TermA() ) );
	A.addTerm( B.copy() );
	B.addTerm( SymbolicTermPtr( new TermA() ) );
	A.addTerm( B.copy() );
	B.addTerm( SymbolicTermPtr( new TermA() ) );
	A.addTerm( B.copy() );
	B.addTerm( SymbolicTermPtr( new TermA() ) );
	A.addTerm( B.copy() );
	ss << A << "    " << truncateOddOrders( A.copy() );
	return ss.str();
}

string AI01() {
	stringstream ss;
	map<int, int> A;
	A[ 3 ] = 2;
	A[ 2 ] = 1;
	A[ 1 ] = 0;
	A[ 5 ] = 4;
	ss << "1: " << getTerminatedContraction( A, 1 ) << "    2: " << getTerminatedContraction( A, 2 ) << "    3: "
		<< getTerminatedContraction( A, 3 ) << "    4: " << getTerminatedContraction( A, 4 ) << "    5: "
		<< getTerminatedContraction( A, 5 );
	return ss.str();
}

string AI02() {
	stringstream ss;
	DeltaContractionSet A;
	A.addContraction( IndexContraction( 0, 1 ) );
	A.addContraction( IndexContraction( 1, 2 ) );
	A.addContraction( IndexContraction( 2, 3 ) );
	ss << A << "    " << constructContractionDictionary( A );
	return ss.str();
}

string AI03() {
	stringstream ss;
	DeltaContractionSet A;
	A.addContraction( IndexContraction( 0, 1 ) );
	A.addContraction( IndexContraction( 2, 3 ) );
	A.addContraction( IndexContraction( 0, 2 ) );
	ss << A << "    " << constructContractionDictionary( A );
	return ss.str();
}

string AI04() {
	stringstream ss;
	DeltaContractionSet A;
	A.addContraction( IndexContraction( 0, 1 ) );
	A.addContraction( IndexContraction( 2, 3 ) );
	A.addContraction( IndexContraction( 3, 4 ) );
	ss << A << "    " << constructContractionDictionary( A );
	return ss.str();
}

string AI05() {
	stringstream ss;
	DeltaContractionSet A;
	A.addContraction( IndexContraction( 0, 2 ) );
	A.addContraction( IndexContraction( 2, 3 ) );
	A.addContraction( IndexContraction( 0, 1 ) );
	ss << A << "    " << constructContractionDictionary( A );
	return ss.str();
}

string AI06() {
	stringstream ss;
	DeltaContractionSet A;
	A.addContraction( IndexContraction( 1, 2 ) );
	A.addContraction( IndexContraction( 0, 1 ) );
	A.addContraction( IndexContraction( 0, 2 ) );
	ss << A << "    " << constructContractionDictionary( A );
	return ss.str();
}

string AI07() {
	stringstream ss;
	DeltaContractionSet A;
	A.addContraction( IndexContraction( 0, 2 ) );
	A.addContraction( IndexContraction( 0, 3 ) );
	A.addContraction( IndexContraction( 3, 2 ) );
	ss << A << "    " << constructContractionDictionary( A );
	return ss.str();
}

string AI08() {
	stringstream ss;
	DeltaContractionSet A;
	A.addContraction( IndexContraction( 0, 1 ) );
	A.addContraction( IndexContraction( 2, 3 ) );
	A.addContraction( IndexContraction( 4, 5 ) );
	A.addContraction( IndexContraction( 6, 7 ) );
	A.addContraction( IndexContraction( 0, 2 ) );
	A.addContraction( IndexContraction( 2, 4 ) );
	A.addContraction( IndexContraction( 4, 6 ) );
	ss << A << "    " << constructContractionDictionary( A );
	return ss.str();
}

string AI09() {
	stringstream ss;
	DeltaContractionSet A;
	A.addContraction( IndexContraction( 1, 2 ) );
	A.addContraction( IndexContraction( 3, 4 ) );
	A.addContraction( IndexContraction( 0, 5 ) );
	A.addContraction( IndexContraction( 0, 2 ) );
	A.addContraction( IndexContraction( 4, 6 ) );
	ss << A << "    " << constructContractionDictionary( A );
	return ss.str();
}

string AI10() {
	stringstream ss;
	DeltaContractionSet A;
	A.addContraction( IndexContraction( 1, 2 ) );
	A.addContraction( IndexContraction( 3, 4 ) );
	A.addContraction( IndexContraction( 5, 6 ) );
	A.addContraction( IndexContraction( 0, 7 ) );
	A.addContraction( IndexContraction( 0, 4 ) );
	A.addContraction( IndexContraction( 2, 6 ) );
	ss << A << "    " << constructContractionDictionary( A );
	return ss.str();
}

string AI11() {
	stringstream ss;
	DeltaContractionSet A;
	A.addContraction( IndexContraction( 1, 2 ) );
	A.addContraction( IndexContraction( 3, 4 ) );
	A.addContraction( IndexContraction( 2, 6 ) );
	A.addContraction( IndexContraction( 5, 6 ) );
	A.addContraction( IndexContraction( 0, 7 ) );
	A.addContraction( IndexContraction( 0, 4 ) );
	ss << A << "    " << constructContractionDictionary( A );
	return ss.str();
}

string AI12() {
	stringstream ss;
	DeltaContractionSet A;
	A.addContraction( IndexContraction( 1, 2 ) );
	A.addContraction( IndexContraction( 3, 4 ) );
	A.addContraction( IndexContraction( 5, 6 ) );
	A.addContraction( IndexContraction( 7, 0 ) );
	A.addContraction( IndexContraction( 2, 4 ) );
	A.addContraction( IndexContraction( 6, 0 ) );
	A.addContraction( IndexContraction( 4, 6 ) );
	ss << A << "    " << constructContractionDictionary( A );
	return ss.str();
}

string AJ01() {
	stringstream ss;
	Sum A;
	Product B;
	MatrixK C;
	C.setIndices( 0, 1 );
	B.addTerm( SymbolicTermPtr( C.copy() ) );
	C.setIndices( 1, 0 );
	B.addTerm( SymbolicTermPtr( C.copy() ) );
	B.addTerm( DeltaPtr( new Delta( 0, 1 ) ) );
	A.addTerm( B.copy() );
	ss << A << "    " << fourierTransformExpression( A.copy() );
	return ss.str();
}

string AK01() {
	stringstream ss;
	Product A;
	A.addTerm( TermAPtr( new TermA() ) );
	A.addTerm( TermAPtr( new TermA() ) );
	vector<IndexContraction> B;
	B.push_back( IndexContraction( 0, 1 ) );
	B.push_back( IndexContraction( 2, 3 ) );
	B.push_back( IndexContraction( 4, 5 ) );
	A.addTerm( FourierSumPtr( new FourierSum( B, 3 ) ) );
	Product C;
	C.addTerm( TermAPtr( new TermA() ) );
	C.addTerm( TermAPtr( new TermA() ) );
	C.addTerm( FourierSumPtr( new FourierSum( B, 3 ) ) );
	ss << A << "    " << C << "    " << areTermsCommon( A.copy(), C.copy() );
	return ss.str();
}

string AK02() {
	stringstream ss;
	Product A;
	A.addTerm( TermAPtr( new TermA() ) );
	A.addTerm( TermAPtr( new TermA() ) );
	A.addTerm( TermAPtr( new TermA() ) );
	vector<IndexContraction> B;
	B.push_back( IndexContraction( 0, 1 ) );
	B.push_back( IndexContraction( 2, 3 ) );
	B.push_back( IndexContraction( 4, 5 ) );
	A.addTerm( FourierSumPtr( new FourierSum( B, 3 ) ) );
	Product C;
	C.addTerm( TermAPtr( new TermA() ) );
	C.addTerm( TermAPtr( new TermA() ) );
	C.addTerm( FourierSumPtr( new FourierSum( B, 3 ) ) );
	ss << A << "    " << C << "    " << areTermsCommon( A.copy(), C.copy() );
	return ss.str();
}

string AK03() {
	stringstream ss;
	Product A;
	A.addTerm( TermAPtr( new TermA() ) );
	A.addTerm( TermAPtr( new TermA() ) );
	A.addTerm( TermAPtr( new TermA() ) );
	vector<IndexContraction> B;
	B.push_back( IndexContraction( 0, 1 ) );
	B.push_back( IndexContraction( 2, 3 ) );
	B.push_back( IndexContraction( 4, 5 ) );
	A.addTerm( FourierSumPtr( new FourierSum( B, 3 ) ) );
	Product C;
	vector<IndexContraction> D;
	D.push_back( IndexContraction( 0, 1 ) );
	D.push_back( IndexContraction( 2, 3 ) );
	D.push_back( IndexContraction( 6, 7 ) );
	C.addTerm( TermAPtr( new TermA() ) );
	C.addTerm( TermAPtr( new TermA() ) );
	C.addTerm( FourierSumPtr( new FourierSum( D, 3 ) ) );
	ss << A << "    " << C << "    " << areTermsCommon( A.copy(), C.copy() );
	return ss.str();
}

string AL01() {
	stringstream ss;
	Product A;
	A.addTerm( TermAPtr( new TermA() ) );
	A.addTerm( TermAPtr( new TermA() ) );
	vector<IndexContraction> B;
	B.push_back( IndexContraction( 0, 1 ) );
	B.push_back( IndexContraction( 2, 3 ) );
	B.push_back( IndexContraction( 4, 5 ) );
	A.addTerm( FourierSumPtr( new FourierSum( B, 3 ) ) );
	Product C;
	vector<IndexContraction> D;
	D.push_back( IndexContraction( 0, 1 ) );
	D.push_back( IndexContraction( 2, 3 ) );
	D.push_back( IndexContraction( 6, 7 ) );
	C.addTerm( TermAPtr( new TermA() ) );
	C.addTerm( TermAPtr( new TermA() ) );
	C.addTerm( FourierSumPtr( new FourierSum( D, 3 ) ) );
	Sum E;
	E.addTerm( A.copy() );
	E.addTerm( C.copy() );
	ss << E << "    ";
	ss << combineLikeTerms( E );
	return ss.str();
}

string AL02() {
	stringstream ss;
	Product A;
	A.addTerm( TermAPtr( new TermA() ) );
	A.addTerm( TermAPtr( new TermA() ) );
	vector<IndexContraction> B;
	B.push_back( IndexContraction( 0, 1 ) );
	B.push_back( IndexContraction( 2, 3 ) );
	B.push_back( IndexContraction( 4, 5 ) );
	A.addTerm( FourierSumPtr( new FourierSum( B, 3 ) ) );
	Product C;
	vector<IndexContraction> D;
	D.push_back( IndexContraction( 0, 1 ) );
	D.push_back( IndexContraction( 2, 3 ) );
	D.push_back( IndexContraction( 6, 7 ) );
	C.addTerm( TermAPtr( new TermA() ) );
	C.addTerm( TermAPtr( new TermA() ) );
	C.addTerm( FourierSumPtr( new FourierSum( D, 3 ) ) );
	Sum E;
	E.addTerm( A.copy() );
	E.addTerm( C.copy() );
	ss << E << "    " << combineLikeTerms( E );
	return ss.str();
}

string AM01() {
	stringstream ss;
	ss << gcd( 6, 4 );
	return ss.str();
}

string AM02() {
	stringstream ss;
	ss << gcd( 3, 12 );
	return ss.str();
}

string AM03() {
	stringstream ss;
	ss << gcd( 3, 7 );
	return ss.str();
}

string AM04() {
	stringstream ss;
	ss << gcd( 164, 88 );
	return ss.str();
}

int main( int argc, char** argv ) {
	cout << "**********************************************************************" << endl;
	cout << "  Amaunet Primary Unit Testing" << endl;
	cout << "**********************************************************************" << endl << endl;

	/*
	 * ********************************************************************
	 * CLASS UNIT TESTS
	 * ********************************************************************
	 */

	/*
	 * Unit test objects unit tests.
	 */

	UnitTest( "i01: GenericTestTerm, Constructor", &i01, "GT_0^0" );

	UnitTest( "i02: GenericTestTerm, getDerivative()", &i02, "GT_0^1" );

	UnitTest( "i03: GenericTestTerm, copy() I", &i03, "GT_1^0 GT_1^0" );

	UnitTest( "i04: GenericTestTerm, copy() II", &i04, "GT_1^0 GT_1^1" );

	UnitTest( "i05: GenericTestTerm, getTermID()", &i05, "g" );

	UnitTest( "i06: initializeStaticReferences(), SINE_PATH_INTEGRALS", &i06, "1: 0 / 1    2: 1 / 2    3: 0 / 1    4: 3 / 8" );

	/*
	 * SymbolicTerm
	 */

	UnitTest( "A01: SymbolicTerm, Default Constructor", &A01, "<invalid_term>" );

	UnitTest( "A02: SymbolicTerm, getDerivative()", &A02, "0" );

	UnitTest( "A03: SymbolicTerm, getTermID()", &A03, "0" );

	UnitTest( "A04: SymbolicTerm, copy()", &A04, "<invalid_term>    1" );

	/*
	 * B: TermA
	 */

	UnitTest( "B01: TermA, Default Constructor", &B01, "A" );

	UnitTest( "B02: TermA, == Overload", &B02, "1" );

	UnitTest( "B03: TermA, copy()", &B03, "A" );

	UnitTest( "B04: TermA, getTermID()", &B04, "A" );

	UnitTest( "B05: TermA, getTermID(), copy()", &B05, "A" );

	/*
	 * C: MatrixM
	 */

	UnitTest( "C01: MatrixM, Default Constructor", &C01, "M" );

	UnitTest( "C02: MatrixM, setAsNonInteracting()", &C02, "M0" );

	UnitTest( "C03: MatrixM, Constructor, flavorLabel", &C03, "M_up" );

	UnitTest( "C04: MatrixM, to_string(), isInteracting = 0", &C04, "M0_up" );

	UnitTest( "C05: MatrixM, getDerivative()", &C05, "dM_up / dA" );

	UnitTest( "C06: MatrixM, Second derivative", &C06, "0" );

	UnitTest( "C07: MatrixM, getDerivative(), isInteracting = 0", &C07, "0" );

	UnitTest( "C08: MatrixM, copy() I", &C08, "M_up M0_up 1 0" );

	UnitTest( "C09: MatrixM, copy() II", &C09, "dM_up / dA" );

	/*
	 * D: CoefficientFloat
	 */

	UnitTest( "D01: CoefficientFloat, to_string(), Constructor", &D01, "0" );

	UnitTest( "D02: CoefficientFloat, getDerivative()", &D02, "3    0" );

	/*
	 * E: MatrixB
	 */

	UnitTest( "E01: MatrixB, Default Constructor", &E01, "B" );

	UnitTest( "E02: MatrixB, setAsNonInteracting()", &E02, "B0" );

	UnitTest( "E03: MatrixB, Constructor, flavorLabel", &E03, "B_up" );

	UnitTest( "E04: MatrixB, to_string(), isInteracting = 0", &E04, "B0_up" );

	UnitTest( "E05: MatrixB, getDerivative()", &E05, " {-1} {B_up} {dM_up / dA} {B_up} " );

	UnitTest( "E06: MatrixB, Second derivative", &E06, " {-1} {-1} {B_up} {dM_up / dA} {B_up} {dM_up / dA} {B_up}  +  {-1} {B_up} {dM_up / dA} {-1} {B_up} {dM_up / dA} {B_up} " );

	UnitTest( "E07: MatrixB, Third derivative", &E07, " {-1} {-1} {-1} {B_up} {dM_up / dA} {B_up} {dM_up / dA} {B_up} {dM_up / dA} {B_up}  +  {-1} {-1} {B_up} {dM_up / dA} {-1} {B_up} {dM_up / dA} {B_up} {dM_up / dA} {B_up}  +  {-1} {-1} {B_up} {dM_up / dA} {B_up} {dM_up / dA} {-1} {B_up} {dM_up / dA} {B_up}  +  {-1} {-1} {B_up} {dM_up / dA} {B_up} {dM_up / dA} {-1} {B_up} {dM_up / dA} {B_up}  +  {-1} {B_up} {dM_up / dA} {-1} {-1} {B_up} {dM_up / dA} {B_up} {dM_up / dA} {B_up}  +  {-1} {B_up} {dM_up / dA} {-1} {B_up} {dM_up / dA} {-1} {B_up} {dM_up / dA} {B_up} " );

	/*
	 * F: MatrixK
	 */

	UnitTest( "F01: MatrixK, Default Constructor", &F01, "K__( 0, 0 )" );

	UnitTest( "F02: MatrixK, Constructor, flavorLabel", &F02, "K_up_( 0, 0 )" );

	UnitTest( "F03: MatrixK, setAsNonInteracting()", &F03, "" );

	UnitTest( "F04: MatrixK, fourierTransform()", &F04, "D_up_( 0, 0 )" );

	/*
	 * G: MatrixS
	 */

	UnitTest( "G01: MatrixS, Default Constructor", &G01, "S_(0, 0)" );

	/*
	 * H: DetM
	 */

	UnitTest( "H01: DetM, Constructor", &H01, "Det[ M_up ]" );

	UnitTest( "H02: DetM, copy(), setAsNonInteracting()", &H02, "Det[ M0_up ]    Det[ M_up ]" );

	UnitTest( "H03: DetM, getDerivative()", &H03, " {Det[ M_up ]} {Trace[  {B_up} {dM_up / dA}  ]} " );

	UnitTest( "H04: DetM, Second Derivative", &H04, " {Det[ M_up ]} {Trace[  {B_up} {dM_up / dA}  ]} {Trace[  {B_up} {dM_up / dA}  ]}  +  {Det[ M_up ]} {Trace[  {-1} {B_up} {dM_up / dA} {B_up} {dM_up / dA}  ]} " );

	UnitTest( "H05: DetM, Third Derivative", &H05, " {Det[ M_up ]} {Trace[  {B_up} {dM_up / dA}  ]} {Trace[  {B_up} {dM_up / dA}  ]} {Trace[  {B_up} {dM_up / dA}  ]}  +  {Det[ M_up ]} {Trace[  {-1} {B_up} {dM_up / dA} {B_up} {dM_up / dA}  ]} {Trace[  {B_up} {dM_up / dA}  ]}  +  {Det[ M_up ]} {Trace[  {B_up} {dM_up / dA}  ]} {Trace[  {-1} {B_up} {dM_up / dA} {B_up} {dM_up / dA}  ]}  +  {Det[ M_up ]} {Trace[  {B_up} {dM_up / dA}  ]} {Trace[  {-1} {B_up} {dM_up / dA} {B_up} {dM_up / dA}  ]}  +  {Det[ M_up ]} {Trace[  {-1} {-1} {B_up} {dM_up / dA} {B_up} {dM_up / dA} {B_up} {dM_up / dA}  +  {-1} {B_up} {dM_up / dA} {-1} {B_up} {dM_up / dA} {B_up} {dM_up / dA}  ]} " );

	/*
	 * I: CoefficientFraction
	 */

	UnitTest( "I01: CoefficientFraction, Constructor, operator<< Overload", &I01, "3 / 7" );

	UnitTest( "I02: CoefficientFraction, getDerivative()", &I02, "5 / 9    0" );

	UnitTest( "I03: CoefficientFraction, Sum, getDerivative(), simplify()", &I03, "GT_0^0 + 3 / 7    GT_0^1    GT_0^1" );

	UnitTest( "I04: CoefficientFraction, reduce() I", &I04, "4 / 8    1 / 2" );

	UnitTest( "I05: CoefficientFraction, reduce() II", &I05, "84 / 174    14 / 29" );

	UnitTest( "I06: CoefficientFraction, reduce() III", &I06, "3 / 7    3 / 7" );

	UnitTest( "I07: CoefficientFraction, operator* Overload", &I07, "3 / 7    4 / 9    4 / 21" );

	UnitTest( "I08: CoefficientFraction, operator+ Overload I", &I08, "1 / 2    1 / 4    3 / 4" );

	UnitTest( "I09: CoefficientFraction, operator+ Overload II", &I09, "5 / 3    8 / 3    13 / 3" );

	UnitTest( "I10: CoefficientFraction, operator+ Overload III", &I10, "7 / 2    8 / 4    11 / 2" );

	UnitTest( "I11: CoefficientFraction, reduce() IV, double", &I11, "5.4 / 7.9    5.4 / 7.9" );

	UnitTest( "I12: CoefficientFraction, operator*( CoefficientFloat ) Overload I", &I12, "1 / 3    5    5 / 3" );

	UnitTest( "I13: CoefficientFraction, operator*( CoefficientFloat ) Overload II", &I13, "5 / 8    2    5 / 4" );

	UnitTest( "I14: CoefficientFraction, operator*( CoefficientFloat ) Overload III", &I14, "5 / 8    -2    -5 / 4" );

	UnitTest( "I15: CoefficientFraction, operator*( CoefficientFloat ) Overload IV", &I15, "1 / 2    -1    -1 / 2    1 / 2" );

	UnitTest( "I16: CoefficientFraction, operator+( CoefficientFloat ) Overload I", &I16, "1 / 2    1    3 / 2" );

	UnitTest( "I17: CoefficientFraction, operator+( CoefficientFloat ) Overload II", &I17, "1 / 2    -1    -1 / 2" );

	UnitTest( "I18: CoefficientFraction, operator+( CoefficientFloat ) Overload III", &I18, "2 / 8    -1    -3 / 4" );

	/*
	 * J: Sum
	 */

	UnitTest( "J01: Sum, Constructor, SymbolicTerm*", &J01, "GT_0^0" );

	UnitTest( "J02: Sum, Default Constructor, addTerm()", &J02, "GT_0^0 + GT_1^0 + GT_2^0" );

	UnitTest( "J03: Sum, getDerivative() I", &J03, "GT_0^1 + GT_1^1 + GT_2^1" );

	UnitTest( "J04: Sum, copy()", &J04, "M_a + M_b + M_c M0_a + M0_b + M0_c" );

	UnitTest( "J05: Sum, simplify() I", &J05, "M_a + M_b + 0 + M_c M_a + M_b + M_c" );

	UnitTest( "J06: Sum, reduceTree() I", &J06, "GT_0^0 + GT_1^0 + GT_2^0 + GT_3^0   2    GT_0^0 + GT_1^0 + GT_2^0 + GT_3^0   4" );

	UnitTest( "J07: Sum, reduceTree() II", &J07, "GT_0^0 +  {GT_1^0} { {GT_2^0} {GT_3^0} }     2    GT_0^0 +  {GT_1^0} {GT_2^0} {GT_3^0}     2" );

	UnitTest( "J08: Sum, reduceTree() III", &J08, " {-1} { {GT_0^0} {GT_1^0} } {GT_2^0}     1     {-1} {GT_0^0} {GT_1^0} {GT_2^0}     1" );

	UnitTest( "J09: Sum, reduceTree() IV", &J09, " { {GT_0^0} {GT_1^0} }     1     {GT_0^0} {GT_1^0}     1" );

	UnitTest( "J10: Sum, reduceTree() V", &J10, " { { {GT_0^0} {GT_1^0} } } {GT_2^0}     1     {GT_0^0} {GT_1^0} {GT_2^0}     1" );

	UnitTest( "J11: Sum, reduceTree() VI", &J11, " { {5} { {3} { {GT_0^0} {GT_1^0} } } }     1     {5} {3} {GT_0^0} {GT_1^0}     1" );

	UnitTest( "J12: Sum, simplify() II, Seg Fault Check, Ending Zero", &J12, "GT_0^0" );

	/*
	 * K: Product
	 */

	UnitTest( "K01: Product, Constructor, SymbolicTerm*", &K01, " {GT_0^0} " );

	UnitTest( "K02: Product, Default Constructor, addTerm()", &K02, " {GT_1^0} {GT_2^0} {GT_3^0} " );

	UnitTest( "K03: Product, setAsNonInteracting()", &K03, " {M0} {M0} {B0} " );

	UnitTest( "K04: Product, containsSum() I", &K04, "1" );

	UnitTest( "K05: Product, containsSum() II", &K05, "0" );

	UnitTest( "K06: Product, getDerivative() I", &K06, " {GT_0^1} {GT_1^0} {GT_2^0}  +  {GT_0^0} {GT_1^1} {GT_2^0}  +  {GT_0^0} {GT_1^0} {GT_2^1} " );

	UnitTest( "K07: Product, simplify() I", &K07, " {0} " );

	UnitTest( "K08: Product, simplify() II", &K08, " {GT_0^0} {GT_1^0} {GT_2^0} " );

	UnitTest( "K09: Product, simplify() III", &K09, " {0} " );

	UnitTest( "K10: Product, Destructor, Memory Leak Check I", &K10, " {GT_0^0} " );

	UnitTest( "K11: Product, reduceTree() I", &K11, " {GT_0^0} {GT_1^0} { {GT_2^0} {GT_3^0} }     3     {GT_0^0} {GT_1^0} {GT_2^0} {GT_3^0}     4" );

	UnitTest( "K12: Product, getExpandedExpr() I", &K12, " {GT_0^0}      {GT_0^0} " );

	UnitTest( "K13: Product, getExpandedExpr() II", &K13, " {GT_0^0 + GT_1^0} {GT_2^0 + GT_3^0}      {GT_0^0} {GT_2^0}  +  {GT_0^0} {GT_3^0}  +  {GT_1^0} {GT_2^0}  +  {GT_1^0} {GT_3^0} " );

	UnitTest( "K14: Product, getExpandedExpr() III, no reduceTree()", &K14, " {GT_5^0} {GT_0^0 + GT_1^0} {GT_2^0 + GT_3^0 + GT_4^0}      { {GT_5^0} {GT_0^0} } {GT_2^0}  +  { {GT_5^0} {GT_0^0} } {GT_3^0}  +  { {GT_5^0} {GT_0^0} } {GT_4^0}  +  { {GT_5^0} {GT_1^0} } {GT_2^0}  +  { {GT_5^0} {GT_1^0} } {GT_3^0}  +  { {GT_5^0} {GT_1^0} } {GT_4^0} " );

	UnitTest( "K15: Product, getExpandedExpr() IV", &K15, " {GT_0^0 + GT_1^0} {GT_5^0} {GT_2^0 + GT_3^0 + GT_4^0}      {GT_0^0} {GT_5^0} {GT_2^0}  +  {GT_0^0} {GT_5^0} {GT_3^0}  +  {GT_0^0} {GT_5^0} {GT_4^0}  +  {GT_1^0} {GT_5^0} {GT_2^0}  +  {GT_1^0} {GT_5^0} {GT_3^0}  +  {GT_1^0} {GT_5^0} {GT_4^0} " );

	UnitTest( "K16: Product, getExpandedExpr() V, reduceTree()", &K16, " { {GT_4^0} {GT_0^0 + GT_1^0} } { {GT_5^0} {GT_2^0 + GT_3^0} }      {GT_4^0} {GT_0^0 + GT_1^0} {GT_5^0} {GT_2^0 + GT_3^0}      {GT_4^0} {GT_0^0} {GT_5^0} {GT_2^0}  +  {GT_4^0} {GT_0^0} {GT_5^0} {GT_3^0}  +  {GT_4^0} {GT_1^0} {GT_5^0} {GT_2^0}  +  {GT_4^0} {GT_1^0} {GT_5^0} {GT_3^0} " );

	UnitTest( "K17: Product, simplify() IV, Seg Fault Check, Ending One", &K17, " {GT_0^0} " );

	UnitTest( "K18: Product, zero()", &K18, " {GT_0^0} {GT_1^0} {GT_2^0}      {0} " );

	/*
	 * L: Trace
	 */

	UnitTest( "L01: Trace, Constructor", &L01, "Trace[ GT_0^0 ]" );

	UnitTest( "L02: Trace, getDerivative() I", &L02, "Trace[ GT_0^1 ]" );

	UnitTest( "L03: Trace, rewriteInKSFormalism() I", &L03, "Trace[  {B_up} {M_up}  ]    Trace[  {K_up_( 0, 0 )} {S_(0, 0)}  ]" );

	UnitTest( "L04: Trace, rewriteInKSFormalism() II", &L04, "Trace[  {B_up} {M_up} {B_dn} {M_dn}  ]    Trace[  {K_up_( 0, 0 )} {S_(0, 0)} {K_dn_( 0, 0 )} {S_(0, 0)}  ]" );

	/*
	 * M: Delta
	 */

	UnitTest( "M01: Delta, Constructor", &M01, "Delta( 0, 1 )" );

	UnitTest( "M02: Delta, Constructor, isDeltaBar()", &M02, "DeltaBar( 0, 1 )" );

	/*
	 * N: FourierSum
	 */

	UnitTest( "N01: FourierSum, Constructor, operator<< Overload", &N01, "FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]" );

	UnitTest( "N02: FourierSum, operator== Overload I", &N02, "FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]    FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]    1" );

	UnitTest( "N03: FourierSum, operator== Overload II", &N03, "FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]    FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 7, 9 ) ]    0" );

	UnitTest( "N04: FourierSum, operator== Overload III", &N04, "FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]    FourierSum[ ( 0, 1 )  ( 4, 5 )  ( 2, 3 ) ]    1" );

	UnitTest( "N05: FourierSum, operator== Overload IV", &N05, "FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]    FourierSum[ ( 0, 1 )  ( 4, 5 )  ( 3, 2 ) ]    0" );

	UnitTest( "N06: FourierSum, operator== Overload V", &N06, "FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]    FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 )  ( 6, 7 ) ]    0" );

	UnitTest( "N07: FourierSum, to_string()", N07, "FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]" );

	UnitTest( "N08: FourierSum, reduceDummyIndices() I", N08, "FourierSum[ ( 0, 1 )  ( 3, 3 )  ( 0, 0 ) ]    FourierSum[ ( 0, 1 )  ( 2, 2 )  ( 0, 0 ) ]" );

	/*
	 * ********************************************************************
	 * METHOD UNIT TESTS
	 * ********************************************************************
	 */

	/*
	 * O: unpackTrivialExpression()
	 */

	UnitTest( "O01: unpackTrivialExpression(), Product I", &O01, "GT_0^0 g" );

	UnitTest( "O02: unpackTrivialExpression(), Product II", &O02, " {GT_0^0} {GT_1^0}  P" );

	UnitTest( "O03: unpackTrivialExpression(), Sum I", &O03, "GT_0^0 g" );

	UnitTest( "O04: unpackTrivialExpression(), Sum II", &O04, "GT_0^0 + GT_1^0 S" );

	UnitTest( "O05: unpackTrivialExpression(), Single Sum in Product", &O05, " {GT_0^0 + GT_1^0 + GT_2^0}     P    GT_0^0 + GT_1^0 + GT_2^0    S" );

	UnitTest( "O06: unpackTrivialExpression(), Nested Expressions", &O06, " {GT_0^0}     GT_0^0" );

	/*
	 * P: isZeroTrace()
	 */

	UnitTest( "P01: isZeroTrace(), Empty Sum", &P01, "1" );

	UnitTest( "P02: isZeroTrace(), Non-empty Sum", &P02, "0" );

	UnitTest( "P03: isZeroTrace(), Empty Product", &P03, "1" );

	UnitTest( "P04: isZeroTrace(), Non-empty Product", &P04, "0" );

	UnitTest( "P05: isZeroTrace(), Non-Trace", &P05, "0" );

	/*
	 * Q: distributeTrace()
	 */

	UnitTest( "Q01: distributeTrace(), Trivial", &Q01, "Trace[ M_up ]" );

	UnitTest( "Q02: distributeTrace(), Simplification I", &Q02, " {-1} {Trace[ M_up ]} " );

	UnitTest( "Q03: distributeTrace(), Simplification II", &Q03, " {-1} {Trace[ M_up ]}  +  {-1} {Trace[ M_dn ]} " );

	UnitTest( "Q04: distributeTrace(), Simplification III", &Q04, "Trace[  {-1} { {3} {M_up}  +  {5} {M_dn} }  ]     {-1} {3} {Trace[ M_up ]}  +  {-1} {5} {Trace[ M_dn ]} " );

	/*
	 * R: distributeAllTraces()
	 */

	UnitTest( "R01: distributeAllTraces(), Distribution I", &R01, " {A} {Trace[ GT_0^0 + GT_1^0 ]}      {A} {Trace[ GT_0^0 ]}  +  {A} {Trace[ GT_1^0 ]} " );

	UnitTest( "R02: distributeAllTraces(), Distribution II", &R02, " {-1} {Trace[  {3} {M_up}  ]} {Trace[  {5} {M_dn}  ]}      {-1} {3} {Trace[ M_up ]} {5} {Trace[ M_dn ]} " );

	UnitTest( "R03: distributeAllTraces(), Distribution III", &R03, "Trace[  {3} {M_up}  ] + Trace[  {5} {M_dn}  ]     {3} {Trace[ M_up ]}  +  {5} {Trace[ M_dn ]} " );

	UnitTest( "R04: distributeAllTraces(), Distribution IV", &R04, " {A} {Trace[  {3} {M_up}  ]}  +  {A} {Trace[  {5} {M_dn}  ]}      {A} {3} {Trace[ M_up ]}  +  {A} {5} {Trace[ M_dn ]} " );

	UnitTest( "R05: distributeAllTraces(), Distribution V", &R05, " {A} {Trace[  {3} {B_up} {M_up}  ]}  +  {A} {Trace[  {5} {B_md} {M_md}  ]}  +  {A} {Trace[  {7} {B_dn} {M_dn}  ]}      {A} {3} {Trace[  {B_up} {M_up}  ]}  +  {A} {5} {Trace[  {B_md} {M_md}  ]}  +  {A} {7} {Trace[  {B_dn} {M_dn}  ]} " );

	/*
	 * S: rewriteSumInKSFormalism()
	 */

	UnitTest( "S01: rewriteSumInKSFormalism() I", &S01, "Trace[  {B_up} {M_up} {B_dn} {M_dn}  ]    Trace[  {K_up_( 0, 0 )} {S_(0, 0)} {K_dn_( 0, 0 )} {S_(0, 0)}  ]" );

	UnitTest( "S02: rewriteSumInKSFormalism() II", &S02, "Trace[  {B_up} {M_up} {B_dn} {M_dn}  ] + Trace[  {B_a} {M_a} {B_a} {M_a}  ] + Trace[  {B_c} {M_d} {B_e} {M_f} {B_g} {M_h}  ]    Trace[  {K_up_( 0, 0 )} {S_(0, 0)} {K_dn_( 0, 0 )} {S_(0, 0)}  ] + Trace[  {K_a_( 0, 0 )} {S_(0, 0)} {K_a_( 0, 0 )} {S_(0, 0)}  ] + Trace[  {K_c_( 0, 0 )} {S_(0, 0)} {K_e_( 0, 0 )} {S_(0, 0)} {K_g_( 0, 0 )} {S_(0, 0)}  ]" );

	/*
	 * T: indexExpression()
	 */

	UnitTest( "T01: indexExpression() I", &T01, " {Trace[  {B_up} {M_up} {B_dn} {M_dn}  ]}      { {K_up_( 0, 1 )} {S_(1, 2)} {K_dn_( 2, 3 )} {S_(3, 0)} } " );

	UnitTest( "T02: indexExpression() II", &T02, " {Trace[  {B_up} {M_up} {B_dn} {M_dn}  ]}  +  {Trace[  {B_up} {M_up} {B_dn} {M_dn} {B_md} {M_md}  ]}      { {K_up_( 0, 1 )} {S_(1, 2)} {K_dn_( 2, 3 )} {S_(3, 0)} }  +  { {K_up_( 0, 1 )} {S_(1, 2)} {K_dn_( 2, 3 )} {S_(3, 4)} {K_md_( 4, 5 )} {S_(5, 0)} } " );

	UnitTest( "T03: indexExpression() III", &T03, " {Trace[  {B_up} {M_up} {B_dn} {M_dn}  ]} {Trace[  {B_up} {M_up}  ]}      { {K_up_( 0, 1 )} {S_(1, 2)} {K_dn_( 2, 3 )} {S_(3, 0)} } { {K_up_( 4, 5 )} {S_(5, 4)} } " );

	/*
	 * U: IndexContraction
	 */

	UnitTest( "U01: IndexContraction, Struct Construction", &U01, "3    7" );

	/*
	 * V: DeltaContractionSet
	 */

	UnitTest( "V01: DeltaContractionSet, Constructor, to_string()", &V01, "[ ( 1, 2 )  ( 3, 5 )  ( 7, 11 ) ]" );

	UnitTest( "V02: DeltaContractionSet, operator<< Overload", &V02, "[ ( 1, 2 )  ( 3, 5 )  ( 7, 11 ) ]" );

	UnitTest( "V03: DeltaContractionSet, orderContractionIndices()", &V03, "[ ( 1, 2 )  ( 5, 3 )  ( 7, 11 )  ( 5, 5 )  ( 6, 5 ) ]    [ ( 1, 2 )  ( 3, 5 )  ( 7, 11 )  ( 5, 5 )  ( 5, 6 ) ]" );

	UnitTest( "V04: DeltaContractionSet, sortContractions() I", &V04, "[ ( 1, 2 )  ( 5, 3 )  ( 4, 7 )  ( 9, 12 )  ( 8, 8 ) ]    [ ( 1, 2 )  ( 4, 7 )  ( 5, 3 )  ( 8, 8 )  ( 9, 12 ) ]" );

	UnitTest( "V05: DeltaContractionSet, sortContractions() II", &V05, "[ ( 1, 2 )  ( 5, 3 )  ( 1, 0 )  ( 8, 12 )  ( 8, 8 ) ]    [ ( 1, 0 )  ( 1, 2 )  ( 5, 3 )  ( 8, 8 )  ( 8, 12 ) ]" );

	/*
	 * X: DeltaSignature
	 */

	UnitTest( "X01: DeltaSignature, Constructor, to_string()", &X01, "[ [ ( 1, 2 )  ( 3, 5 )  ( 7, 11 ) ]  [ ( 13, 17 )  ( 19, 31 ) ]  [] ]" );

	UnitTest( "X02: DeltaSignature, operator<< Overload", &X02, "[ [ ( 1, 2 )  ( 3, 5 )  ( 7, 11 ) ]  [ ( 13, 17 )  ( 19, 31 ) ]  [] ]" );

	/*
	 * Y: getDeltaSignature
	 */

	UnitTest( "Y01: getDeltaSignature() I, (2)", &Y01, "[ ( 0, 1 ) ]    []" );

	UnitTest( "Y02: getDeltaSignature() II, (4)", &Y02, "[ ( 0, 1 )  ( 1, 2 )  ( 2, 3 ) ]    []" );

	UnitTest( "Y03: getDeltaSignature() III, (2,2)", &Y03, "[ ( 0, 1 )  ( 2, 3 ) ]    [ ( 1, 2 ) ]" );

	UnitTest( "Y04: getDeltaSignature() IV, (6)", &Y04, "[ ( 0, 1 )  ( 1, 2 )  ( 2, 3 )  ( 3, 4 )  ( 4, 5 ) ]    []" );

	UnitTest( "Y05: getDeltaSignature() V, (4,2)", &Y05, "[ ( 0, 1 )  ( 1, 2 )  ( 2, 3 )  ( 4, 5 ) ]    [ ( 3, 4 ) ]" );

	UnitTest( "Y06: getDeltaSignature() VI, (2,2,2)", &Y06, "[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]    [ ( 1, 2 )  ( 3, 4 ) ]" );

	UnitTest( "Y07: getDeltaSignature() VII, (8)", &Y07, "[ ( 0, 1 )  ( 1, 2 )  ( 2, 3 )  ( 3, 4 )  ( 4, 5 )  ( 5, 6 )  ( 6, 7 ) ]    []" );

	UnitTest( "Y08: getDeltaSignature() VIII, (6,2)", &Y08, "[ ( 0, 1 )  ( 1, 2 )  ( 2, 3 )  ( 3, 4 )  ( 4, 5 )  ( 6, 7 ) ]    [ ( 5, 6 ) ]" );

	UnitTest( "Y09: getDeltaSignature() IX, (4,4)", &Y09, "[ ( 0, 1 )  ( 1, 2 )  ( 2, 3 )  ( 4, 5 )  ( 5, 6 )  ( 6, 7 ) ]    [ ( 3, 4 ) ]" );

	UnitTest( "Y10: getDeltaSignature() X, (4,2,2)", &Y10, "[ ( 0, 1 )  ( 1, 2 )  ( 2, 3 )  ( 4, 5 )  ( 6, 7 ) ]    [ ( 3, 4 )  ( 5, 6 ) ]" );

	UnitTest( "Y11: getDeltaSignature() XI, (2,2,2,2)", &Y11, "[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 )  ( 6, 7 ) ]    [ ( 1, 2 )  ( 3, 4 )  ( 5, 6 ) ]" );

	/*
	 * Z: combinations()
	 */

	UnitTest( "Z01: combinations() I, operator<< Overload", &Z01, "[ [  1  ] ]" );

	UnitTest( "Z02: combinations() II", &Z02, "[ [  1  ]  [  2  ]  [  3  ] ]" );

	UnitTest( "Z03: combinations() III", &Z03, "[ [  1  2  ]  [  1  3  ]  [  1  4  ]  [  1  5  ]  [  2  3  ]  [  2  4  ]  [  2  5  ]  [  3  4  ]  [  3  5  ]  [  4  5  ] ]    10" );

	UnitTest( "Z04: combinations() IV", &Z04, "[ [  1  2  3  ]  [  1  2  4  ]  [  1  2  5  ]  [  1  3  4  ]  [  1  3  5  ]  [  1  4  5  ]  [  2  3  4  ]  [  2  3  5  ]  [  2  4  5  ]  [  3  4  5  ] ]    10" );

	UnitTest( "Z05: combinations() V", &Z05, "[ [  1  2  3  4  ]  [  1  2  3  5  ]  [  1  2  4  5  ]  [  1  3  4  5  ]  [  2  3  4  5  ] ]    5" );

	UnitTest( "Z06: combinations() VI", &Z06, "[ [  1  2  3  4  5  ] ]    1" );

	UnitTest( "Z07: combinations() VII", &Z07, "[ [  1  ]  [  2  ]  [  3  ]  [  4  ]  [  5  ]  [  6  ]  [  7  ] ]    7" );

	UnitTest( "Z08: combinations() VIII", &Z08, "[ [  1  2  ]  [  1  3  ]  [  1  4  ]  [  1  5  ]  [  1  6  ]  [  1  7  ]  [  2  3  ]  [  2  4  ]  [  2  5  ]  [  2  6  ]  [  2  7  ]  [  3  4  ]  [  3  5  ]  [  3  6  ]  [  3  7  ]  [  4  5  ]  [  4  6  ]  [  4  7  ]  [  5  6  ]  [  5  7  ]  [  6  7  ] ]    21" );

	UnitTest( "Z09: combinations() IX", &Z09, "[ [  1  2  3  ]  [  1  2  4  ]  [  1  2  5  ]  [  1  2  6  ]  [  1  2  7  ]  [  1  3  4  ]  [  1  3  5  ]  [  1  3  6  ]  [  1  3  7  ]  [  1  4  5  ]  [  1  4  6  ]  [  1  4  7  ]  [  1  5  6  ]  [  1  5  7  ]  [  1  6  7  ]  [  2  3  4  ]  [  2  3  5  ]  [  2  3  6  ]  [  2  3  7  ]  [  2  4  5  ]  [  2  4  6  ]  [  2  4  7  ]  [  2  5  6  ]  [  2  5  7  ]  [  2  6  7  ]  [  3  4  5  ]  [  3  4  6  ]  [  3  4  7  ]  [  3  5  6  ]  [  3  5  7  ]  [  3  6  7  ]  [  4  5  6  ]  [  4  5  7  ]  [  4  6  7  ]  [  5  6  7  ] ]    35" );

	UnitTest( "Z10: combinations() X", &Z10, "[ [  1  2  3  4  ]  [  1  2  3  5  ]  [  1  2  3  6  ]  [  1  2  3  7  ]  [  1  2  4  5  ]  [  1  2  4  6  ]  [  1  2  4  7  ]  [  1  2  5  6  ]  [  1  2  5  7  ]  [  1  2  6  7  ]  [  1  3  4  5  ]  [  1  3  4  6  ]  [  1  3  4  7  ]  [  1  3  5  6  ]  [  1  3  5  7  ]  [  1  3  6  7  ]  [  1  4  5  6  ]  [  1  4  5  7  ]  [  1  4  6  7  ]  [  1  5  6  7  ]  [  2  3  4  5  ]  [  2  3  4  6  ]  [  2  3  4  7  ]  [  2  3  5  6  ]  [  2  3  5  7  ]  [  2  3  6  7  ]  [  2  4  5  6  ]  [  2  4  5  7  ]  [  2  4  6  7  ]  [  2  5  6  7  ]  [  3  4  5  6  ]  [  3  4  5  7  ]  [  3  4  6  7  ]  [  3  5  6  7  ]  [  4  5  6  7  ] ]    35" );

	UnitTest( "Z11: combinations() XI", &Z11, "[ [  1  2  3  4  5  ]  [  1  2  3  4  6  ]  [  1  2  3  4  7  ]  [  1  2  3  5  6  ]  [  1  2  3  5  7  ]  [  1  2  3  6  7  ]  [  1  2  4  5  6  ]  [  1  2  4  5  7  ]  [  1  2  4  6  7  ]  [  1  2  5  6  7  ]  [  1  3  4  5  6  ]  [  1  3  4  5  7  ]  [  1  3  4  6  7  ]  [  1  3  5  6  7  ]  [  1  4  5  6  7  ]  [  2  3  4  5  6  ]  [  2  3  4  5  7  ]  [  2  3  4  6  7  ]  [  2  3  5  6  7  ]  [  2  4  5  6  7  ]  [  3  4  5  6  7  ] ]    21" );

	UnitTest( "Z12: combinations() XII, Choose Zero, Empty Vector", &Z12, "[]" );

	UnitTest( "Z13: combinations() XIII, Choose Zero, Non-empty Vector", &Z13, "[]" );

	/*
	 * getIndexPermutations()
	 */

	UnitTest( "AA01: getIndexPermutations() I, (2), operator<< Overload", &AA01, "[ [  0  1  ] ]" );

	UnitTest( "AA02: getIndexPermutations() II, (4)", &AA02, "[ [  0  1  2  3  ] ]" );

	UnitTest( "AA03: getIndexPermutations() III, (2,2)", &AA03, "[ [  0  1  2  3  ]  [  0  2  1  3  ]  [  0  3  1  2  ]  [  1  2  0  3  ]  [  1  3  0  2  ]  [  2  3  0  1  ] ]" );

	UnitTest( "AA04: getIndexPermutations() IV, (6)", &AA04, "[ [  0  1  2  3  4  5  ] ]" );

	// TODO: Verify.
	UnitTest( "AA05: getIndexPermutations() V, (4,2)", &AA05, "[ [  0  1  2  3  4  5  ]  [  0  1  2  4  3  5  ]  [  0  1  2  5  3  4  ]  [  0  1  3  4  2  5  ]  [  0  1  3  5  2  4  ]  [  0  1  4  5  2  3  ]  [  0  2  3  4  1  5  ]  [  0  2  3  5  1  4  ]  [  0  2  4  5  1  3  ]  [  0  3  4  5  1  2  ]  [  1  2  3  4  0  5  ]  [  1  2  3  5  0  4  ]  [  1  2  4  5  0  3  ]  [  1  3  4  5  0  2  ]  [  2  3  4  5  0  1  ] ]" );

	// TODO: Verify.
	UnitTest( "AA06: getIndexPermutations() VI, (2,2,2)", &AA06, "[ [  0  1  2  3  4  5  ]  [  0  1  2  4  3  5  ]  [  0  1  2  5  3  4  ]  [  0  1  3  4  2  5  ]  [  0  1  3  5  2  4  ]  [  0  1  4  5  2  3  ]  [  0  2  1  3  4  5  ]  [  0  2  1  4  3  5  ]  [  0  2  1  5  3  4  ]  [  0  2  3  4  1  5  ]  [  0  2  3  5  1  4  ]  [  0  2  4  5  1  3  ]  [  0  3  1  2  4  5  ]  [  0  3  1  4  2  5  ]  [  0  3  1  5  2  4  ]  [  0  3  2  4  1  5  ]  [  0  3  2  5  1  4  ]  [  0  3  4  5  1  2  ]  [  0  4  1  2  3  5  ]  [  0  4  1  3  2  5  ]  [  0  4  1  5  2  3  ]  [  0  4  2  3  1  5  ]  [  0  4  2  5  1  3  ]  [  0  4  3  5  1  2  ]  [  0  5  1  2  3  4  ]  [  0  5  1  3  2  4  ]  [  0  5  1  4  2  3  ]  [  0  5  2  3  1  4  ]  [  0  5  2  4  1  3  ]  [  0  5  3  4  1  2  ]  [  1  2  0  3  4  5  ]  [  1  2  0  4  3  5  ]  [  1  2  0  5  3  4  ]  [  1  2  3  4  0  5  ]  [  1  2  3  5  0  4  ]  [  1  2  4  5  0  3  ]  [  1  3  0  2  4  5  ]  [  1  3  0  4  2  5  ]  [  1  3  0  5  2  4  ]  [  1  3  2  4  0  5  ]  [  1  3  2  5  0  4  ]  [  1  3  4  5  0  2  ]  [  1  4  0  2  3  5  ]  [  1  4  0  3  2  5  ]  [  1  4  0  5  2  3  ]  [  1  4  2  3  0  5  ]  [  1  4  2  5  0  3  ]  [  1  4  3  5  0  2  ]  [  1  5  0  2  3  4  ]  [  1  5  0  3  2  4  ]  [  1  5  0  4  2  3  ]  [  1  5  2  3  0  4  ]  [  1  5  2  4  0  3  ]  [  1  5  3  4  0  2  ]  [  2  3  0  1  4  5  ]  [  2  3  0  4  1  5  ]  [  2  3  0  5  1  4  ]  [  2  3  1  4  0  5  ]  [  2  3  1  5  0  4  ]  [  2  3  4  5  0  1  ]  [  2  4  0  1  3  5  ]  [  2  4  0  3  1  5  ]  [  2  4  0  5  1  3  ]  [  2  4  1  3  0  5  ]  [  2  4  1  5  0  3  ]  [  2  4  3  5  0  1  ]  [  2  5  0  1  3  4  ]  [  2  5  0  3  1  4  ]  [  2  5  0  4  1  3  ]  [  2  5  1  3  0  4  ]  [  2  5  1  4  0  3  ]  [  2  5  3  4  0  1  ]  [  3  4  0  1  2  5  ]  [  3  4  0  2  1  5  ]  [  3  4  0  5  1  2  ]  [  3  4  1  2  0  5  ]  [  3  4  1  5  0  2  ]  [  3  4  2  5  0  1  ]  [  3  5  0  1  2  4  ]  [  3  5  0  2  1  4  ]  [  3  5  0  4  1  2  ]  [  3  5  1  2  0  4  ]  [  3  5  1  4  0  2  ]  [  3  5  2  4  0  1  ]  [  4  5  0  1  2  3  ]  [  4  5  0  2  1  3  ]  [  4  5  0  3  1  2  ]  [  4  5  1  2  0  3  ]  [  4  5  1  3  0  2  ]  [  4  5  2  3  0  1  ] ]" );
	/*
	 * generateSignaturePermutations()
	 */

	UnitTest( "AB01: generateSignaturePermutations() I, (2), operator<< Overload", &AB01, "[ { [ ( 0, 1 ) ] | [] } ]" );

	UnitTest( "AB02: generateSignaturePermutations() II, (4)", &AB02, "[ { [ ( 0, 1 )  ( 1, 2 )  ( 2, 3 ) ] | [] } ]" );

	UnitTest( "AB03: generateSignaturePermutations() III, (2,2)", &AB03, "[ { [ ( 0, 1 )  ( 2, 3 ) ] | [ ( 1, 2 ) ] }  { [ ( 0, 2 )  ( 1, 3 ) ] | [ ( 2, 1 ) ] }  { [ ( 0, 3 )  ( 1, 2 ) ] | [ ( 3, 1 ) ] } ]" );

	UnitTest( "AB04: generateSignaturePermutations() IV, (6)", &AB04, "[ { [ ( 0, 1 )  ( 1, 2 )  ( 2, 3 )  ( 3, 4 )  ( 4, 5 ) ] | [] } ]" );

	UnitTest( "AB05: generateSignaturePermutations() V, (4,2)", &AB05, "[ { [ ( 0, 1 )  ( 1, 2 )  ( 2, 3 )  ( 4, 5 ) ] | [ ( 3, 4 ) ] }  { [ ( 0, 1 )  ( 1, 2 )  ( 2, 4 )  ( 3, 5 ) ] | [ ( 4, 3 ) ] }  { [ ( 0, 1 )  ( 1, 2 )  ( 2, 5 )  ( 3, 4 ) ] | [ ( 5, 3 ) ] }  { [ ( 0, 1 )  ( 1, 3 )  ( 3, 4 )  ( 2, 5 ) ] | [ ( 4, 2 ) ] }  { [ ( 0, 1 )  ( 1, 3 )  ( 3, 5 )  ( 2, 4 ) ] | [ ( 5, 2 ) ] }  { [ ( 0, 1 )  ( 1, 4 )  ( 4, 5 )  ( 2, 3 ) ] | [ ( 5, 2 ) ] }  { [ ( 0, 2 )  ( 2, 3 )  ( 3, 4 )  ( 1, 5 ) ] | [ ( 4, 1 ) ] }  { [ ( 0, 2 )  ( 2, 3 )  ( 3, 5 )  ( 1, 4 ) ] | [ ( 5, 1 ) ] }  { [ ( 0, 2 )  ( 2, 4 )  ( 4, 5 )  ( 1, 3 ) ] | [ ( 5, 1 ) ] }  { [ ( 0, 3 )  ( 3, 4 )  ( 4, 5 )  ( 1, 2 ) ] | [ ( 5, 1 ) ] }  { [ ( 1, 2 )  ( 2, 3 )  ( 3, 4 )  ( 0, 5 ) ] | [ ( 4, 0 ) ] }  { [ ( 1, 2 )  ( 2, 3 )  ( 3, 5 )  ( 0, 4 ) ] | [ ( 5, 0 ) ] }  { [ ( 1, 2 )  ( 2, 4 )  ( 4, 5 )  ( 0, 3 ) ] | [ ( 5, 0 ) ] }  { [ ( 1, 3 )  ( 3, 4 )  ( 4, 5 )  ( 0, 2 ) ] | [ ( 5, 0 ) ] }  { [ ( 2, 3 )  ( 3, 4 )  ( 4, 5 )  ( 0, 1 ) ] | [ ( 5, 0 ) ] } ]" );

	UnitTest( "AB06: generateSignaturePermutations() VI, (2,2,2)", &AB06, "[ { [ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ] | [ ( 1, 2 )  ( 3, 4 ) ] }  { [ ( 0, 1 )  ( 2, 4 )  ( 3, 5 ) ] | [ ( 1, 2 )  ( 4, 3 ) ] }  { [ ( 0, 1 )  ( 2, 5 )  ( 3, 4 ) ] | [ ( 1, 2 )  ( 5, 3 ) ] }  { [ ( 0, 2 )  ( 1, 3 )  ( 4, 5 ) ] | [ ( 2, 1 )  ( 3, 4 ) ] }  { [ ( 0, 2 )  ( 1, 4 )  ( 3, 5 ) ] | [ ( 2, 1 )  ( 4, 3 ) ] }  { [ ( 0, 2 )  ( 1, 5 )  ( 3, 4 ) ] | [ ( 2, 1 )  ( 5, 3 ) ] }  { [ ( 0, 3 )  ( 1, 2 )  ( 4, 5 ) ] | [ ( 3, 1 )  ( 2, 4 ) ] }  { [ ( 0, 3 )  ( 1, 4 )  ( 2, 5 ) ] | [ ( 3, 1 )  ( 4, 2 ) ] }  { [ ( 0, 3 )  ( 1, 5 )  ( 2, 4 ) ] | [ ( 3, 1 )  ( 5, 2 ) ] }  { [ ( 0, 4 )  ( 1, 2 )  ( 3, 5 ) ] | [ ( 4, 1 )  ( 2, 3 ) ] }  { [ ( 0, 4 )  ( 1, 3 )  ( 2, 5 ) ] | [ ( 4, 1 )  ( 3, 2 ) ] }  { [ ( 0, 4 )  ( 1, 5 )  ( 2, 3 ) ] | [ ( 4, 1 )  ( 5, 2 ) ] }  { [ ( 0, 5 )  ( 1, 2 )  ( 3, 4 ) ] | [ ( 5, 1 )  ( 2, 3 ) ] }  { [ ( 0, 5 )  ( 1, 3 )  ( 2, 4 ) ] | [ ( 5, 1 )  ( 3, 2 ) ] }  { [ ( 0, 5 )  ( 1, 4 )  ( 2, 3 ) ] | [ ( 5, 1 )  ( 4, 2 ) ] } ]" );

	/*
	 * TotalSignature
	 */

	UnitTest( "AC01: TotalSignature, areSignaturesDegenerate() I", &AC01, "[ ( 0, 1 ) ]    [ ( 2, 3 ) ]    0" );

	UnitTest( "AC02: TotalSignature, areSignaturesDegenerate() II", &AC02, "[ ( 0, 1 ) ]    [ ( 2, 1 ) ]    0" );

	UnitTest( "AC03: TotalSignature, areSignaturesDegenerate() III", &AC03, "[ ( 2, 3 )  ( 1, 5 ) ]    [ ( 1, 2 )  ( 5, 4 ) ]    0" );

	UnitTest( "AC04: TotalSignature, areSignaturesDegenerate() III, Callee switched", &AC04, "[ ( 2, 3 )  ( 1, 5 ) ]    [ ( 1, 2 )  ( 5, 4 ) ]    0" );

	UnitTest( "AC05: TotalSignature, areSignaturesDegenerate() IV", &AC05, "[ ( 1, 2 )  ( 3, 4 ) ]    [ ( 3, 4 )  ( 1, 2 ) ]    1" );

	/*
	 * calculateAllContractions()
	 */

	UnitTest( "AD01: calculateAllContractions(), n = 2", &AD01, "[ [  2  ] ]" );

	UnitTest( "AD02: calculateAllContractions(), n = 4", &AD02, "[ [  4  ]  [  2  2  ] ]" );

	UnitTest( "AD03: calculateAllContractions(), n = 6", &AD03, "[ [  6  ]  [  2  4  ]  [  2  2  2  ] ]" );

	UnitTest( "AD04: calculateAllContractions(), n = 8", &AD04, "[ [  8  ]  [  2  6  ]  [  2  2  4  ]  [  2  2  2  2  ] ]" );

	UnitTest( "AD05: calculateAllContractions(), n = 10", &AD05, "[ [  10  ]  [  2  8  ]  [  2  2  6  ]  [  2  2  2  4  ]  [  2  2  2  2  2  ] ]" );

	/*
	 * generateCoordinateSpacePathIntegral()
	 */

	UnitTest( "AE01: generateCoordinateSpacePathIntegral(), n = 2", &AE01, " {1 / 2} {Delta( 0, 1 )} " );

	UnitTest( "AE02: generateCoordinateSpacePathIntegral(), n = 4", &AE02, " {3 / 8} {Delta( 0, 1 )} {Delta( 1, 2 )} {Delta( 2, 3 )}  +  {1 / 2} {1 / 2} { {Delta( 0, 1 )} {Delta( 2, 3 )} {1 +  {-1} {Delta( 1, 2 )} }  +  {Delta( 0, 2 )} {Delta( 1, 3 )} {1 +  {-1} {Delta( 2, 1 )} }  +  {Delta( 0, 3 )} {Delta( 1, 2 )} {1 +  {-1} {Delta( 3, 1 )} } } " );

	UnitTest( "AE03: generateCoordinateSpacePathIntegral(), n = 6", &AE03, " {5 / 16} {Delta( 0, 1 )} {Delta( 1, 2 )} {Delta( 2, 3 )} {Delta( 3, 4 )} {Delta( 4, 5 )}  +  {1 / 2} {3 / 8} { {Delta( 0, 1 )} {Delta( 2, 3 )} {Delta( 3, 4 )} {Delta( 4, 5 )} {1 +  {-1} {Delta( 1, 2 )} }  +  {Delta( 0, 2 )} {Delta( 1, 3 )} {Delta( 3, 4 )} {Delta( 4, 5 )} {1 +  {-1} {Delta( 2, 1 )} }  +  {Delta( 0, 3 )} {Delta( 1, 2 )} {Delta( 2, 4 )} {Delta( 4, 5 )} {1 +  {-1} {Delta( 3, 1 )} }  +  {Delta( 0, 4 )} {Delta( 1, 2 )} {Delta( 2, 3 )} {Delta( 3, 5 )} {1 +  {-1} {Delta( 4, 1 )} }  +  {Delta( 0, 5 )} {Delta( 1, 2 )} {Delta( 2, 3 )} {Delta( 3, 4 )} {1 +  {-1} {Delta( 5, 1 )} }  +  {Delta( 1, 2 )} {Delta( 0, 3 )} {Delta( 3, 4 )} {Delta( 4, 5 )} {1 +  {-1} {Delta( 2, 0 )} }  +  {Delta( 1, 3 )} {Delta( 0, 2 )} {Delta( 2, 4 )} {Delta( 4, 5 )} {1 +  {-1} {Delta( 3, 0 )} }  +  {Delta( 1, 4 )} {Delta( 0, 2 )} {Delta( 2, 3 )} {Delta( 3, 5 )} {1 +  {-1} {Delta( 4, 0 )} }  +  {Delta( 1, 5 )} {Delta( 0, 2 )} {Delta( 2, 3 )} {Delta( 3, 4 )} {1 +  {-1} {Delta( 5, 0 )} }  +  {Delta( 2, 3 )} {Delta( 0, 1 )} {Delta( 1, 4 )} {Delta( 4, 5 )} {1 +  {-1} {Delta( 3, 0 )} }  +  {Delta( 2, 4 )} {Delta( 0, 1 )} {Delta( 1, 3 )} {Delta( 3, 5 )} {1 +  {-1} {Delta( 4, 0 )} }  +  {Delta( 2, 5 )} {Delta( 0, 1 )} {Delta( 1, 3 )} {Delta( 3, 4 )} {1 +  {-1} {Delta( 5, 0 )} }  +  {Delta( 3, 4 )} {Delta( 0, 1 )} {Delta( 1, 2 )} {Delta( 2, 5 )} {1 +  {-1} {Delta( 4, 0 )} }  +  {Delta( 3, 5 )} {Delta( 0, 1 )} {Delta( 1, 2 )} {Delta( 2, 4 )} {1 +  {-1} {Delta( 5, 0 )} }  +  {Delta( 4, 5 )} {Delta( 0, 1 )} {Delta( 1, 2 )} {Delta( 2, 3 )} {1 +  {-1} {Delta( 5, 0 )} } }  +  {1 / 2} {1 / 2} {1 / 2} { {Delta( 0, 1 )} {Delta( 2, 3 )} {Delta( 4, 5 )} {1 +  {-1} {Delta( 1, 2 )} } {1 +  {-1} {Delta( 3, 4 )} }  +  {Delta( 0, 1 )} {Delta( 2, 4 )} {Delta( 3, 5 )} {1 +  {-1} {Delta( 1, 2 )} } {1 +  {-1} {Delta( 4, 3 )} }  +  {Delta( 0, 1 )} {Delta( 2, 5 )} {Delta( 3, 4 )} {1 +  {-1} {Delta( 1, 2 )} } {1 +  {-1} {Delta( 5, 3 )} }  +  {Delta( 0, 2 )} {Delta( 1, 3 )} {Delta( 4, 5 )} {1 +  {-1} {Delta( 2, 1 )} } {1 +  {-1} {Delta( 3, 4 )} }  +  {Delta( 0, 2 )} {Delta( 1, 4 )} {Delta( 3, 5 )} {1 +  {-1} {Delta( 2, 1 )} } {1 +  {-1} {Delta( 4, 3 )} }  +  {Delta( 0, 2 )} {Delta( 1, 5 )} {Delta( 3, 4 )} {1 +  {-1} {Delta( 2, 1 )} } {1 +  {-1} {Delta( 5, 3 )} }  +  {Delta( 0, 3 )} {Delta( 1, 2 )} {Delta( 4, 5 )} {1 +  {-1} {Delta( 3, 1 )} } {1 +  {-1} {Delta( 2, 4 )} }  +  {Delta( 0, 3 )} {Delta( 1, 4 )} {Delta( 2, 5 )} {1 +  {-1} {Delta( 3, 1 )} } {1 +  {-1} {Delta( 4, 2 )} }  +  {Delta( 0, 3 )} {Delta( 1, 5 )} {Delta( 2, 4 )} {1 +  {-1} {Delta( 3, 1 )} } {1 +  {-1} {Delta( 5, 2 )} }  +  {Delta( 0, 4 )} {Delta( 1, 2 )} {Delta( 3, 5 )} {1 +  {-1} {Delta( 4, 1 )} } {1 +  {-1} {Delta( 2, 3 )} }  +  {Delta( 0, 4 )} {Delta( 1, 3 )} {Delta( 2, 5 )} {1 +  {-1} {Delta( 4, 1 )} } {1 +  {-1} {Delta( 3, 2 )} }  +  {Delta( 0, 4 )} {Delta( 1, 5 )} {Delta( 2, 3 )} {1 +  {-1} {Delta( 4, 1 )} } {1 +  {-1} {Delta( 5, 2 )} }  +  {Delta( 0, 5 )} {Delta( 1, 2 )} {Delta( 3, 4 )} {1 +  {-1} {Delta( 5, 1 )} } {1 +  {-1} {Delta( 2, 3 )} }  +  {Delta( 0, 5 )} {Delta( 1, 3 )} {Delta( 2, 4 )} {1 +  {-1} {Delta( 5, 1 )} } {1 +  {-1} {Delta( 3, 2 )} }  +  {Delta( 0, 5 )} {Delta( 1, 4 )} {Delta( 2, 3 )} {1 +  {-1} {Delta( 5, 1 )} } {1 +  {-1} {Delta( 4, 2 )} } } " );

	/*
	 * pathIntegrateExpression()
	 */

	UnitTest( "AF01: pathIntegrateExpression() I", &AF01, " {K__( 0, 1 )} {S_(1, 0)} {K__( 2, 3 )} {S_(3, 2)}      {K__( 0, 1 )} {Delta( 1, 0 )} {K__( 2, 3 )} {Delta( 3, 2 )} { {1 / 2} {Delta( 0, 2 )} } " );

	UnitTest( "AF02: pathIntegrateExpression() II", &AF02, " {K__( 0, 1 )} {S_(1, 0)} {K__( 2, 3 )} {S_(3, 2)} {K__( 4, 5 )} {S_(5, 4)} {K__( 6, 7 )} {S_(7, 6)}      {K__( 0, 1 )} {Delta( 1, 0 )} {K__( 2, 3 )} {Delta( 3, 2 )} {K__( 4, 5 )} {Delta( 5, 4 )} {K__( 6, 7 )} {Delta( 7, 6 )} { {3 / 8} {Delta( 0, 2 )} {Delta( 2, 4 )} {Delta( 4, 6 )}  +  {1 / 2} {1 / 2} {Delta( 0, 2 )} {Delta( 4, 6 )} {1}  +  {1 / 2} {1 / 2} {Delta( 0, 2 )} {Delta( 4, 6 )} {-1} {Delta( 2, 4 )}  +  {1 / 2} {1 / 2} {Delta( 0, 4 )} {Delta( 2, 6 )} {1}  +  {1 / 2} {1 / 2} {Delta( 0, 4 )} {Delta( 2, 6 )} {-1} {Delta( 4, 2 )}  +  {1 / 2} {1 / 2} {Delta( 0, 6 )} {Delta( 2, 4 )} {1}  +  {1 / 2} {1 / 2} {Delta( 0, 6 )} {Delta( 2, 4 )} {-1} {Delta( 6, 2 )} } " );

	/*
	 * truncateAOrder()
	 */

	UnitTest( "AG01: truncateAOrder() I", &AG01, " {A}  +  {A} {A}      {A} " );

	/*
	 * truncateOddOrders()
	 */

	UnitTest( "AH01: truncateOddOrders() I", &AH01, " {A}  +  {A} {A}  +  {A} {A} {A}  +  {A} {A} {A} {A}  +  {A} {A} {A} {A} {A}  +  {A} {A} {A} {A} {A} {A}      {A} {A}  +  {A} {A} {A} {A}  +  {A} {A} {A} {A} {A} {A} " );

	/*
	 * getTerminatedContraction(), constructContractionDictionary()
	 */

	UnitTest( "AI01: getTerminatedContraction()" , &AI01, "1: 0    2: 0    3: 0    4: 4    5: 4" );

	UnitTest( "AI02: constructContractionDictionary() I, operator<< Overload", &AI02, "[ ( 0, 1 )  ( 1, 2 )  ( 2, 3 ) ]    [ 0 : 0  1 : 0  2 : 0  3 : 0 ]" );

	UnitTest( "AI03: constructContractionDictionary() II", &AI03, "[ ( 0, 1 )  ( 2, 3 )  ( 0, 2 ) ]    [ 0 : 0  1 : 0  2 : 0  3 : 0 ]" );

	UnitTest( "AI04: constructContractionDictionary() III", &AI04, "[ ( 0, 1 )  ( 2, 3 )  ( 3, 4 ) ]    [ 0 : 0  1 : 0  2 : 2  3 : 2  4 : 2 ]" );

	UnitTest( "AI05: constructContractionDictionary() IV", &AI05, "[ ( 0, 2 )  ( 2, 3 )  ( 0, 1 ) ]    [ 0 : 0  1 : 0  2 : 0  3 : 0 ]" );

	UnitTest( "AI06: constructContractionDictionary() V", &AI06, "[ ( 1, 2 )  ( 0, 1 )  ( 0, 2 ) ]    [ 0 : 0  1 : 0  2 : 0 ]" );

	UnitTest( "AI07: constructContractionDictionary() VI", &AI07, "[ ( 0, 2 )  ( 0, 3 )  ( 3, 2 ) ]    [ 0 : 0  2 : 0  3 : 0 ]" );

	UnitTest( "AI08: constructContractionDictionary() VII", &AI08, "[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 )  ( 6, 7 )  ( 0, 2 )  ( 2, 4 )  ( 4, 6 ) ]    [ 0 : 0  1 : 0  2 : 0  3 : 0  4 : 0  5 : 0  6 : 0  7 : 0 ]" );

	UnitTest( "AI09: constructContractionDictionary() VIII", &AI09, "[ ( 1, 2 )  ( 3, 4 )  ( 0, 5 )  ( 0, 2 )  ( 4, 6 ) ]    [ 0 : 0  1 : 0  2 : 0  3 : 3  4 : 3  5 : 0  6 : 3 ]" );

	UnitTest( "AI10: constructContractionDictionary() IX", &AI10, "[ ( 1, 2 )  ( 3, 4 )  ( 5, 6 )  ( 0, 7 )  ( 0, 4 )  ( 2, 6 ) ]    [ 0 : 0  1 : 1  2 : 1  3 : 0  4 : 0  5 : 1  6 : 1  7 : 0 ]" );

	UnitTest( "AI11: constructContractionDictionary() X", &AI11, "[ ( 1, 2 )  ( 3, 4 )  ( 2, 6 )  ( 5, 6 )  ( 0, 7 )  ( 0, 4 ) ]    [ 0 : 0  1 : 1  2 : 1  3 : 0  4 : 0  5 : 1  6 : 1  7 : 0 ]" );

	UnitTest( "AI12: constructContractionDictionary() XI", &AI12, "[ ( 1, 2 )  ( 3, 4 )  ( 5, 6 )  ( 7, 0 )  ( 2, 4 )  ( 6, 0 )  ( 4, 6 ) ]    [ 0 : 0  1 : 0  2 : 0  3 : 0  4 : 0  5 : 0  6 : 0  7 : 0 ]" );

	/*
	 * fourierTransformExpression()
	 */

	UnitTest( "AJ01: fourierTransformExpression() I", &AJ01, " {K__( 0, 1 )} {K__( 1, 0 )} {Delta( 0, 1 )}      {K__( 0, 1 )} {K__( 1, 0 )} {FourierSum[ ( 0, 0 )  ( 0, 0 ) ]} " );

	/*
	 * areTermsCommon(),
	 */

	UnitTest( "AK01: areTermsCommon() I", &AK01, " {A} {A} {FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]}      {A} {A} {FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]}     1" );

	UnitTest( "AK02: areTermsCommon() II", &AK02, " {A} {A} {A} {FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]}      {A} {A} {FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]}     0" );

	UnitTest( "AK03: areTermsCommon() III", &AK03, " {A} {A} {A} {FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]}      {A} {A} {FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 6, 7 ) ]}     0" );

	/*
	 * combineLikeTerms()
	 */

	UnitTest( "AL01: combineLikeTerms() I", &AL01, " {A} {A} {FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]}  +  {A} {A} {FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 6, 7 ) ]}      {A} {A} {FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 4, 5 ) ]} {1}  +  {A} {A} {FourierSum[ ( 0, 1 )  ( 2, 3 )  ( 6, 7 ) ]} {1} " );

	/*
	 * gcd()
	 */

	UnitTest( "AM01: gcd() I", &AM01, "2" );

	UnitTest( "AM02: gcd() II", &AM02, "3" );

	UnitTest( "AM03: gcd() III", &AM03, "1" );

	UnitTest( "AM04: gcd() IV", &AM04, "4" );

	cout << "----------------------------------------------------------------------" << endl;
	cout << UnitTest::passedTests << " tests PASSED, " << UnitTest::failedTests << " tests FAILED." << endl;
}
