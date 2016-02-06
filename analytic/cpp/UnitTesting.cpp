/* ***********************************************************************
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * Numerical CUDA Perturbation Theory Fourier Transform Evaluation
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Primary Unit Testing
 *
 * v. 0.1		02 Feb 2016
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at CHapel Hill
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
	GenericTestTerm* B = A.copy();
	ss << *B << " " << A;
	return ss.str();
}

string i04() {
	stringstream ss;
	GenericTestTerm A = GenericTestTerm( 1, 0 );
	GenericTestTerm* B = A.copy();
	Sum C = A.getDerivative();
	ss << *B << " " << C;
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
	TermA* B = A.copy();
	ss << *B;
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
	MatrixM* B = A.copy();
	B->setAsNonInteracting();
	ss << A << " " << *B << " " << A.isTermInteracting() << " " << B->isTermInteracting();
	return ss.str();
}

string C09() {
	stringstream ss;
	MatrixM A = MatrixM( "up" );
	MatrixM* B = A.copy();
	Sum C = B->getDerivative();
	ss << C;
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
	Sum B = A.getDerivative();
	B = B.getDerivative();
	ss << B;
	return ss.str();
}

string J01() {
	stringstream ss;
	Sum A = Sum( new GenericTestTerm(0,0) );
	ss << A;
	return ss.str();
}

string J02() {
	stringstream ss;
	Sum A = Sum();
	A.addTerm( new GenericTestTerm(0,0) );
	A.addTerm( new GenericTestTerm(1,0) );
	A.addTerm( new GenericTestTerm(2,0) );
	ss << A;
	return ss.str();
}

string J03() {
	stringstream ss;
	Sum A = Sum();
	A.addTerm( new GenericTestTerm(0,0) );
	A.addTerm( new GenericTestTerm(1,0) );
	A.addTerm( new GenericTestTerm(2,0) );
	A = A.getDerivative();
	ss << A;
	return ss.str();
}

string J04() {
	stringstream ss;
	Sum A = Sum();
	A.addTerm( new MatrixM( "a" ) );
	A.addTerm( new MatrixM( "b" ) );
	A.addTerm( new MatrixM( "c" ) );
	Sum* B = A.copy();
	B->setAsNonInteracting();
	ss << A << " " << *B;
	return ss.str();
}

string J05() {
	stringstream ss;
	Sum A = Sum();
	A.addTerm( new MatrixM( "a" ) );
	A.addTerm( new MatrixM( "b" ) );
	A.addTerm( new CoefficientFloat( 0.0 ) );
	A.addTerm( new MatrixM( "c" ) );
	ss << A << " ";
	A.simplify();
	ss << A;
	return ss.str();
}

string K01() {
	stringstream ss;
	Product A = Product( new GenericTestTerm(0,0) );
	ss << A;
	return ss.str();
}

string K02() {
	stringstream ss;
	Product A = Product();
	A.addTerm( new GenericTestTerm(1,0) );
	A.addTerm( new GenericTestTerm(2,0) );
	A.addTerm( new GenericTestTerm(3,0) );
	ss << A;
	return ss.str();
}

string K03() {
	stringstream ss;
	Product A = Product();
	A.addTerm( new MatrixM() );
	A.addTerm( new MatrixM() );
	A.addTerm( new MatrixB() );
	A.setAsNonInteracting();
	ss << A;
	return ss.str();
}

string K04() {
	stringstream ss;
	Product A = Product();
	A.addTerm( new MatrixM() );
	A.addTerm( new Sum() );
	A.addTerm( new MatrixB() );
	ss << A.containsSum();
	return ss.str();
}

string K05() {
	stringstream ss;
	Product A = Product();
	A.addTerm( new MatrixM() );
	A.addTerm( new TermA() );
	A.addTerm( new MatrixB() );
	ss << A.containsSum();
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

	/*
	 * SymbolicTerm
	 */

	UnitTest( "A01: SymbolicTerm, Default Constructor", &A01, "<invalid_term>" );

	UnitTest( "A02: SymbolicTerm, getDerivative()", &A02, "0" );

	/*
	 * B: TermA
	 */

	UnitTest( "B01: TermA, Default Constructor", &B01, "A" );

	UnitTest( "B02: TermA, == Overload", &B02, "1" );

	UnitTest( "B03: TermA, copy()", &B03, "A" );

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

	UnitTest( "E06: MatrixB, Second derivative", &E06, "" );

	/*
	 * F: Matrix K
	 */

	/*
	 * G: MatrixS
	 */

	/*
	 * H: DetM
	 */

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

	/*
	 * K: Product
	 */

	UnitTest( "K01: Product, Constructor, SymbolicTerm*", &K01, " {GT_0^0} " );

	UnitTest( "K02: Product, Default Constructor, addTerm()", &K02, " {GT_1^0} {GT_2^0} {GT_3^0} " );

	UnitTest( "K03: Product, setAsNonInteracting()", &K03, " {M0} {M0} {B0} " );

	UnitTest( "K04: Product, containsSum() I", &K04, "1" );

	UnitTest( "K05: Product, containsSum() II", &K05, "0" );

	/*
	 * L: Trace
	 */

	/*
	 * M: Delta
	 */

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

	/*
	 * P: isZeroTrace()
	 */

	cout << "----------------------------------------------------------------------" << endl;
	cout << UnitTest::passedTests << " tests PASSED, " << UnitTest::failedTests << " tests FAILED." << endl;
}
