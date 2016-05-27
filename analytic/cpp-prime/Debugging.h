/* ***********************************************************************
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * High-order Perturbation Theory Analytics
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Debugging Utility Classes and Functions Header
 *
 * v. 0.1		17 May 2016
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at Chapel Hill
 * ***********************************************************************
 */

#ifndef DEBUGGING_H_
#define DEBUGGING_H_

#include <string>
#include <vector>
#include "PTSymbolicObjects.h"

/*
 * ***********************************************************************
 * CLASS DECLARATIONS
 * ***********************************************************************
 */

class DebugTracer : public SymbolicTerm {

public:

	DebugTracer();

	DebugTracer( int thisCounter );

	const std::string to_string() const;

	SymbolicTermPtr copy();

	static int nextCounter;

	int getCounter();

private:

	int counter;
};

class GenericTestTerm : public SymbolicTerm {
public:

	GenericTestTerm( int thisId );

	const std::string to_string() const;

	SymbolicTermPtr copy();

	int id;
};

/*
 * ***********************************************************************
 * DEBUGGING FUNCTIONS
 * ***********************************************************************
 */

void injectDebuggingTracers( Sum &expr );

Sum handpickTerms( Sum &expr, std::vector<int> termIDs );

#endif /* DEBUGGING_H_ */