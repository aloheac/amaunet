#include <iostream>
#include <cuda_runtime.h>
#include "PTSymbolicObjects.h"

using namespace std;

int main( int argc, char** argv ) {
	const string WELCOME_HEADER =
			"**********************************************************************\n"
			"    amaunet\n"
			"    --------------------------------------------------------\n"
			"    high-order lattice perturbation theory\n"
			"    for non-relativistic quantum matter\n\n"
			"    Numerical CUDA Perturbation Theory Fourier Transform Evaluation\n"
			"**********************************************************************";

	const string PHYSICAL_SYSTEM =
			"Weak-coupling Expansion for Fermionic Contact Interactions";

	const string VERSION_STRING = "v. 0.1 (alpha)";

	const string BUILD_DATE = "26 Jan 2016";

	const string COMMIT_ID = "7038d07";

	cout << WELCOME_HEADER << endl << endl;
	cout << PHYSICAL_SYSTEM << endl;
	cout << VERSION_STRING << "\t\t" << BUILD_DATE << "\t\t" << COMMIT_ID << endl << endl;

	int NX = 5;
	double BETA = 8;
	int ORDER_IN_A = 4;
	double TAU = 0.05;

	cout << "Loaded parameters:" << endl;
	cout << "\tNX:\t\t" << NX << endl;
	cout << "\tBETA:\t\t" << BETA << endl;
	cout << "\tTAU:\t\t" << TAU << endl;
	cout << "\tORDER_IN_A:\t" << ORDER_IN_A << endl;
	cout << endl;

	int NTAU = BETA / TAU;

	cout << "Calculated parameters:" << endl;
	cout << "\tNTAU:\t\t" << NTAU << endl;
	cout << endl;

	cout << "Checking for valid CUDA devices..." << endl;
	int numDevices;
	cudaGetDeviceCount( &numDevices );
	if ( numDevices > 0 ) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties( &prop, 0 );
		cout << "GPU 0: " << prop.name << endl;
		cout << "       " << "Compute Capability: " << prop.major << '.' << prop.minor << endl;
	} else {
		cout << "***ERROR: (CUDA) No valid CUDA device found." << endl;
	}
	cout << endl;

	cout << "Expression to evaluate:" << endl;
	std::string dat = pt_util::loadExpressionFile( "/home/loheac/cuda-workspace/AmaunetEval/src/DATA.dat" );
	ExpressionInterpreter expinterp = ExpressionInterpreter();
	Sum expr = expinterp.parseExpression( dat );
	cout << expr << endl;

	cout << "Number of terms to evaluate:\t" << expr.getNumberOfTerms() << endl;
	int i = 0;

	for ( vector<Product>::const_iterator iter = expr.products.begin(); iter != expr.products.end(); iter++ ) {
		cout << "\tTerm " << i << "\t\torder: " << (*iter).getOrderInA() << ", unique index count: " << (*iter).getNumberOfUniqueIndices() << endl;
	}

	Factor A = TermA( 7 );
	Factor* pA = &A;
	TermA* ptr = dynamic_cast<TermA*>(pA);
	if (ptr) {
		cout << "ptr good" << endl;
	} else {
		cout << "ptr bad: " << ptr << endl;
	}
	cout << "_______" << endl << ptr << endl;
}




