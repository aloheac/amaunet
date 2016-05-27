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
#include <iostream>
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
 * GenericTestTerm
 */

GenericTestTerm::GenericTestTerm( int thisId ) : SymbolicTerm() {
	id = thisId;
	termID = TermTypes::GENERIC_TEST_TERM ;
}

const string GenericTestTerm::to_string() const {
	stringstream ss;
	ss << "GT_" << id;
	return ss.str();
}

SymbolicTermPtr GenericTestTerm::copy() {
	GenericTestTermPtr cpy( new GenericTestTerm( 0 ) );
	cpy->id = id;
	cpy->termID = TermTypes::GENERIC_TEST_TERM;

	return SymbolicTermPtr( static_pointer_cast<SymbolicTerm>( cpy ) );
}

std::ostream& operator<<( std::ostream& os, const GenericTestTerm &st ) {
	os << st.to_string();
	return os;
}

/*
 * ***********************************************************************
 * DEBUGGING FUNCTIONS
 * ***********************************************************************
 */

void injectDebuggingTracers( Sum &expr ) {
	for ( vector<SymbolicTermPtr>::iterator iter = expr.getIteratorBegin(); iter != expr.getIteratorEnd(); ++iter ) {
		if ( (*iter)->getTermID() != TermTypes::PRODUCT ) {
			if ( (*iter)->to_string() != "0" and (*iter)->to_string() != "1"  and (*iter)->to_string() != "1 / 0"  and (*iter)->to_string() != "1 / 1" ) {
				cout << "***WARNING: (WA1) A term other then a product, zero, or one was encountered when injecting debugging tracers. The solution may still be correct, but should be inspected." << endl;
			}

			(*iter) = Product( *iter ).copy();
		}

		ProductPtr castProduct = static_pointer_cast<Product>( *iter );
		castProduct->addTerm( SymbolicTermPtr( new DebugTracer() ) );
	}
}

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