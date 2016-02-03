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

class GenericTestTerm : public SymbolicTerm {
public:

	GenericTestTerm( int thisId, int thisDerivativeOrder ) : SymbolicTerm() {
		id = thisId;
		derivativeOrder = thisDerivativeOrder;
	}

	const string to_string() const {
		stringstream ss;
		ss << "GT_" << id << "^" << derivativeOrder;
		return ss.str();
	}

	Sum* getDerivative() {
		//SymbolicTerm* D = new GenericTestTerm( id, derivativeOrder + 1 )
		Sum* D = new Sum();
		return D;
	}

	int id;

	int derivativeOrder;
};

std::ostream& operator<<( std::ostream& os, const GenericTestTerm &st ) {
	os << st.to_string();
	return os;
}

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

string A01() {
	stringstream ss;
	TermA A = TermA();
	ss << A;
	return ss.str();
}

int main( int argc, char** argv ) {
	cout << "**********************************************************************" << endl;
	cout << "  Amaunet Primary Unit Testing" << endl;
	cout << "**********************************************************************" << endl << endl;

	/*
	 * Unit test objects unit tests.
	 */

	UnitTest UTi01 = UnitTest( "I01: GenericTestTerm I", &i01, "GT_0^0" );

	/*
	 * A: TermA
	 */

	UnitTest UTA01 = UnitTest( "A01: TermA I", &A01, "A" );

	cout << "----------------------------------------------------------------------" << endl;
	cout << UnitTest::passedTests << " tests PASSED, " << UnitTest::failedTests << " tests FAILED." << endl;
}
