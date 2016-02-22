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
#include "PTSymbolicObjects.h"
#include <sstream>
#include <string>

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

	UnitTest ("i05: GenericTestTerm, getTermID()", &i05, "g" );

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

	cout << "----------------------------------------------------------------------" << endl;
	cout << UnitTest::passedTests << " tests PASSED, " << UnitTest::failedTests << " tests FAILED." << endl;
}
