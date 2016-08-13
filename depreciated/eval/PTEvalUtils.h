/*
 * AmaunetEV: Numerical Evaluation of High-order Lattice Perturbation Theory
 *
 * Utility Classes and Functions Header
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at Chapel Hill
 */

#ifndef PTEVALUTILS_H_
#define PTEVALUTILS_H_

#include <cstdio>
#include <sstream>
#include <fstream>
#include <vector>

namespace pt_util {

class PTSystemParameters {
public:
	int NX;
	int NTAU;
	double BARE_COUPLING;
	double TAU;
	double BETA;
};

std::string loadExpressionFile( char* filename );

std::string str( int );

std::string str( double );

std::string str( std::vector<int> );

} /* namespace pt_util */

#endif /* PTEVALUTILS_H_ */
