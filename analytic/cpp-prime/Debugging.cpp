/*
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * High-order Perturbation Theory Analytics
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Debugging Utility Classes and Functions Header Implementation
 *
 * v. 0.1		17 May 2016
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at Chapel Hill
 */

#include <sstream>
#include "Debugging.h"

using namespace std;

/*
 * ***********************************************************************
 * TYPEDEF DECLARATIONS
 * ***********************************************************************
 */

typedef std::shared_ptr<DebugTracer> DebugTracerPtr;

/*
 * ***********************************************************************
 * CLASS IMPLEMENTATIONS
 * ***********************************************************************
 */

/*
 * DebugTracer
 */

int DebugTracer::nextCounter = 0;  // Static field of class DebugTracer.

DebugTracer::DebugTracer() {
	counter = nextCounter;
	nextCounter++;
	termID = TermTypes::DEBUG_TRACE;
}

DebugTracer::DebugTracer( int thisCounter ) {
	counter = thisCounter;
	termID = TermTypes::DEBUG_TRACE;
}

const string DebugTracer::to_string() const {
	stringstream ss;
	ss << "T" << counter;
	return ss.str();
}

SymbolicTermPtr DebugTracer::copy() {
	return SymbolicTermPtr( new DebugTracer( counter ) );
}

int DebugTracer::getCounter() {
	return counter;
}

/*
 * ***********************************************************************
 * DEBUGGING FUNCTIONS
 * ***********************************************************************
 */

Sum handpickTerms( Sum &expr, vector<int> termIDs ) {
	Sum handpickedSum;

	for ( vector<int>::iterator id = termIDs.begin(); id != termIDs.end(); ++id ) {
		for ( vector<SymbolicTermPtr>::iterator term = expr.getIteratorBegin(); term != expr.getIteratorEnd(); ++term ) {
			ProductPtr castProduct = static_pointer_cast<Product>( *term );
			for ( vector<SymbolicTermPtr>::iterator factor = castProduct->getIteratorBegin(); factor != castProduct->getIteratorEnd(); ++factor ) {
				if ( (*factor)->getTermID() == TermTypes::DEBUG_TRACE ) {
					DebugTracerPtr castDebugTracer = static_pointer_cast<DebugTracer>( *factor );
					if ( castDebugTracer->getCounter() == *id ) {
						handpickedSum.addTerm( (*term)->copy() );
					}
				}
			}
		}
	}

	return handpickedSum;
}