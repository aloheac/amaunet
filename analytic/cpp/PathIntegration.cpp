/* ***********************************************************************
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * High-order Perturbation Theory Analytics
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Symbolic Path Integration and Fourier Transform Functions for
 * Contact Interactions Implementation
 *
 * v. 0.1		22 Feb 2016
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at Chapel Hill
 * ***********************************************************************
 */

#include <sstream>
#include "PathIntegration.h"

using namespace std;

/*
 * IndexContraction
 */

IndexContraction::IndexContraction() {
    i = 0;
    j = 0;
}

IndexContraction::IndexContraction(int a, int b) {
    i = a;
    j = b;
}

/*
 * DeltaContractionSet
 */

DeltaContractionSet::DeltaContractionSet() { }  // Default constructor is sufficient.

void DeltaContractionSet::addContraction( IndexContraction newContraction ) {
    contractions.push_back( newContraction );
}

unsigned int DeltaContractionSet::getNumContractions() {
    return (int)contractions.size();
}

string DeltaContractionSet::to_string() const {
    stringstream ss;
    ss << "[";
    for ( vector<IndexContraction>::const_iterator iter = contractions.begin(); iter !=contractions.end(); ++iter ) {
        ss << " ( " << (*iter).i << ", " << (*iter).j << " ) ";
    }
    ss << "]";

    return ss.str();
}

/*
 * DeltaSignature
 */

DeltaSignature::DeltaSignature() { }  // Default constructor is sufficient.

void DeltaSignature::addContractionSet( DeltaContractionSet contractionSet ) {
    signatureContractions.push_back( contractionSet );
}

string DeltaSignature::to_string() const {
    stringstream ss;
    ss << "[";
    for ( vector<DeltaContractionSet>::const_iterator iter = signatureContractions.begin(); iter != signatureContractions.end(); ++iter ) {
        ss << " " << (*iter).to_string() << " ";
    }
    ss << "]";

    return ss.str();
}

map<int, CoefficientFraction> Amaunet::SINE_PATH_INTEGRALS;  // Initialized by initializeStaticReferences().

void initializeStaticReferences() {
    using namespace Amaunet;

    // Initialize reference values for the result of the definite integral of various powers of sine:
    //
    //      \int_{ -\pi, \pi} dx \sin^n( x )
    //
    // where here the key of the dictionary is the integer n, and the value is the numerical result
    // of the integral expressed as a fraction (of type CoefficientFraction). Notice that the integral
    // for odd values of n vanish.
    SINE_PATH_INTEGRALS = map<int, CoefficientFraction>();
    SINE_PATH_INTEGRALS[ 0 ] = CoefficientFraction( 0, 1 );
    SINE_PATH_INTEGRALS[ 1 ] = CoefficientFraction( 0, 1 );
    SINE_PATH_INTEGRALS[ 2 ] = CoefficientFraction( 1, 2 );
    SINE_PATH_INTEGRALS[ 3 ] = CoefficientFraction( 0, 1 );
    SINE_PATH_INTEGRALS[ 4 ] = CoefficientFraction( 3, 8 );
    SINE_PATH_INTEGRALS[ 5 ] = CoefficientFraction( 0, 1 );
    SINE_PATH_INTEGRALS[ 6 ] = CoefficientFraction( 5, 16 );
    SINE_PATH_INTEGRALS[ 7 ] = CoefficientFraction( 0, 1 );
    SINE_PATH_INTEGRALS[ 8 ] = CoefficientFraction( 35, 128 );
    SINE_PATH_INTEGRALS[ 9 ] = CoefficientFraction( 0, 1 );
    SINE_PATH_INTEGRALS[ 10 ] = CoefficientFraction( 63, 256 );
}

TotalSignature getDeltaSignature( vector<int> contraction ) {
    DeltaContractionSet deltas;
    DeltaContractionSet deltaBars;
    int nextDeltaIndex = 0;

    for ( int i = 0; i < contraction.size(); i++ ) {
        for ( int j = nextDeltaIndex; j < nextDeltaIndex + contraction[i] - 1; j++ ) {
            deltas.addContraction( IndexContraction( j, j + 1 ) );
        }

        if ( i != contraction.size() - 1 ) {
            deltaBars.addContraction( IndexContraction( nextDeltaIndex + contraction[i] - 1, nextDeltaIndex + contraction[i] ) );
        }

        nextDeltaIndex += contraction[i];
    }

    TotalSignature signature;
    signature.deltas = deltas;
    signature.deltaBars = deltaBars;

    return signature;
}

//def combination( lst, k ):
//if k == 1:
//combos = []
//    for e in lst:
//    combos.append( [e] )
//    return combos
//    elif len( lst ) == k:
//    return [lst]
//    else: # len( lst ) > k
//    subcombinations = combination( lst[1:], k-1 )
//    combos = []
//        for sc in subcombinations:
//        combos.append( [lst[0]] + sc )
//
//        subcombinations = combination( lst[1:], k )
//        for sc in subcombinations:
//        combos.append( sc )
//
//        return combos

vector< vector<int> > combinations( vector<int> list, int k ) {
    // TODO: Throw proper exceptions here.
    if ( k == 1 ) {  // Simply return the list with each element placed in a vector<int>. n choose 1 is always n.
        vector< vector<int> > combos;
        for( vector<int>::iterator element = list.begin(); element != list.end(); ++element ) {
            vector<int> nextCombination;
            nextCombination.push_back( *element );
            combos.push_back( nextCombination );
        }

        return combos;
    } else if ( list.size() == k ) {  // n choose n is always 1. Simply return the list inside a vector<int>/
        vector< vector<int> > combos;
        combos.push_back( list );

        return combos;
    } else {  // list.size() > k. We have a non-trivial case; recursively call.
        int firstElement = list[ 0 ];
        list.erase( list.begin() );  // Remove first element of list. All references to list below depend only on the
                                     // tail of the list.
        vector< vector<int> > combos;

        // Consider the combinations with the first element and choose k-1 combinations with the tail of the list.
        vector< vector<int> > subcombinations = combinations( list, k - 1 );
        for ( vector< vector<int> >::iterator sc = subcombinations.begin(); sc != subcombinations.end(); ++sc ) {
            vector<int> nextCombination;
            nextCombination.push_back( firstElement );
            nextCombination.insert( nextCombination.end(), (*sc).begin(), (*sc).end() );
            combos.push_back( nextCombination );
        }

        // Consider choose k combinations with the tail of the list; the head element is not included.
        subcombinations = combinations( list, k );
        for ( vector< vector<int> >::iterator sc = subcombinations.begin(); sc != subcombinations.end(); ++sc ) {
            combos.push_back( *sc );
        }

        return combos;
    }
}
/*
 * ***********************************************************************
 * INPUT REDIRECTION OPERATOR OVERLOADS
 * ***********************************************************************
 */

ostream& operator<<( std::ostream& os, const DeltaContractionSet &obj ) {
    os << obj.to_string();
    return os;
}

ostream& operator<<( std::ostream& os, const DeltaSignature &obj ) {
    os << obj.to_string();
    return os;
}

ostream& operator<<( std::ostream& os, const vector< vector<int> > &obj ) {
    os << "[";
    for( vector< vector<int> >::const_iterator vec = obj.begin(); vec != obj.end(); ++vec ) {
        os << " [ ";
        for ( vector<int>::const_iterator element = (*vec).begin(); element != (*vec).end(); ++element ) {
            os << " " << *element << " ";
        }
        os << " ] ";
    }
    os << "]";
}