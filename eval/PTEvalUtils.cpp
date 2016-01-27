/*
 * AmaunetEV: Numerical Evaluation of High-order Lattice Perturbation Theory
 *
 * Utility Classes and Functions Header
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at Chapel Hill
 */

#include "PTEvalUtils.h"
#include <iostream>

using namespace std;

std::string pt_util::str( vector<int> vec ) {
	std::stringstream ss;
	ss << "( ";
	for ( int i = 0; i < vec.size(); i++ ) {
		ss << vec[i] << " ";
	}
	ss << ")";
	return ss.str();
}

std::string pt_util::str( int obj ) {
	std::stringstream ss;
	ss << obj;
	return ss.str();
}

std::string pt_util::str( double obj ) {
	std::stringstream ss;
	ss << obj;
	return ss.str();
}

std::string pt_util::loadExpressionFile( char* filename ) {
	std::ifstream fs;
	fs.open( filename );
	std::string line;
	std::string expr;
	expr = "";

	if ( fs.is_open() ) {
		while ( fs.good() ) {
			fs >> line;
			expr += line;
		}
	} else {
		cout << "***ERROR: Expression input file '" << filename << "' could not be read." << endl;
	}

	return expr;
}
