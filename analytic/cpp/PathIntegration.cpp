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
#include <assert.h>
#include <algorithm>
#include <map>
#include <boost/bimap.hpp>
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

bool IndexContraction::operator<( const IndexContraction& rhs ) const {
    if ( i == rhs.i ) {
        return j < rhs.j;
    } else {
        return i < rhs.i;
    }
}

bool IndexContraction::operator>( const IndexContraction& rhs ) const {
    if ( i == rhs.i ) {
        return j > rhs.j;
    } else {
        return i > rhs.i;
    }
}



bool IndexContraction::operator==( const IndexContraction& rhs ) const {
    return i == rhs.i and j == rhs.j or i == rhs.j and j == rhs.i;
}

bool IndexContraction::containsIndex( int index ) {
    return i == index or j == index;
}

/*
 * DeltaContractionSet
 */

DeltaContractionSet::DeltaContractionSet() { }  // Default constructor is sufficient.

void DeltaContractionSet::addContraction( IndexContraction newContraction ) {
    contractions.push_back( newContraction );
}

void DeltaContractionSet::addContractionSet( DeltaContractionSet set ) {
    contractions.insert( contractions.end(), set.contractions.begin(), set.contractions.end() );
}

unsigned int DeltaContractionSet::getNumContractions() {
    return (int)contractions.size();
}

bool DeltaContractionSet::containsIndex( int index ) {
    for ( vector<IndexContraction>::iterator contraction = contractions.begin(); contraction != contractions.end(); ++contraction ) {
        if ( (*contraction).containsIndex( index ) ) return true;
    }

    return false;
}

void DeltaContractionSet::orderContractionIndices() {
    for ( vector<IndexContraction>::iterator indexPair = contractions.begin(); indexPair != contractions.end(); ++indexPair ) {
        if ( indexPair->i > indexPair->j ) {
            // Swap i and j.
            int tmp = indexPair->i;
            indexPair->i = indexPair->j;
            indexPair->j = tmp;
        }
    }
}

void DeltaContractionSet::sortContractions() {
    sort( contractions.begin(), contractions.end() );
}

std::vector<IndexContraction>::iterator DeltaContractionSet::getIteratorBegin() {
    return contractions.begin();
}

std::vector<IndexContraction>::iterator DeltaContractionSet::getIteratorEnd() {
    return contractions.end();
}

std::vector<IndexContraction>::const_iterator DeltaContractionSet::getIteratorBegin() const {
    return contractions.begin();
}

std::vector<IndexContraction>::const_iterator DeltaContractionSet::getIteratorEnd() const {
    return contractions.end();
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

bool DeltaContractionSet::operator==( const DeltaContractionSet &rhs ) const {
    // If the vector sizes of the left- and right-hand sides are not equal, the two object will not be equivalent.
    if ( contractions.size() != rhs.contractions.size() ) return false;

    // It is known that lhs.size() == rhs.size().
    for ( int i = 0; i < contractions.size(); i++ ) {
        // Overload for operator== of types IndexContraction is implemented.
        if ( not ( contractions[ i ] == rhs.contractions[ i ] ) ) return false;
    }

    return true;
}

/*
 * TotalSignature
 */

TotalSignature::TotalSignature() {
    deltas = DeltaContractionSet();     // If the default constructor is called, initialize deltas and deltaBars with
    deltaBars = DeltaContractionSet();  // empty DeltaContractionSet objects.
}

bool TotalSignature::areSignaturesDegenerate( const TotalSignature &other) {
    set<IndexContraction> thisDeltaContractions, otherDeltaContractions;

    for ( vector<IndexContraction>::const_iterator indexPair = deltas.getIteratorBegin(); indexPair != deltas.getIteratorEnd(); ++indexPair ) {
        thisDeltaContractions.insert( IndexContraction( *indexPair ) );
    }

    for ( vector<IndexContraction>::const_iterator indexPair = other.deltas.getIteratorBegin(); indexPair != other.deltas.getIteratorEnd(); ++indexPair ) {
        otherDeltaContractions.insert( IndexContraction( *indexPair ) );
    }

    return thisDeltaContractions == otherDeltaContractions;
}

bool TotalSignature::isValidSignature() {
    for ( vector<IndexContraction>::iterator indexPair = deltas.getIteratorBegin(); indexPair != deltas.getIteratorEnd(); ++indexPair ) {
        if ( indexPair->i == indexPair-> j) return false;
    }

    for ( vector<IndexContraction>::iterator indexPair = deltaBars.getIteratorBegin(); indexPair != deltaBars.getIteratorEnd(); ++indexPair ) {
        if ( indexPair->i == indexPair-> j) return false;
    }

    return true;
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

vector< vector<int> > combinations( vector<int> list, int k ) {
    assert(k <= list.size());

    // TODO: Throw proper exceptions here.
    if ( k == 0 ) {
        return vector< vector<int> >();  // Return empty vector, choose nothing.

    } else if ( k == 1 ) {  // Simply return the list with each element placed in a vector<int>. n choose 1 is always n.
        vector< vector<int> > combos;
        for( vector<int>::iterator element = list.begin(); element != list.end(); ++element ) {
            vector<int> nextCombination;
            nextCombination.push_back( *element );
            combos.push_back( nextCombination );
        }

        return combos;
    } else if ( list.size() == k ) {  // n choose n is always 1. Simply return the list inside a vector<int>.
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

vector< vector<int> > getIndexPermutations( vector<int> contraction, vector<int> list ) {
    // Declare data structure to be returned.
    vector< vector<int> > indexPermutations;

    // Check for the trivial base case. If met, return result.
    if ( contraction.size() == 1 ) {
        indexPermutations.push_back( list );
        return indexPermutations;
    }

    // Get the total signature (which includes \delta and \bar{\delta} index contractions) corresponding to the passed
    // contraction.
    TotalSignature signature = getDeltaSignature( contraction );

    vector<int> difference;
    vector<int>::iterator iter;
    int groupSize = contraction[0];  // Work with first element of contraction, recursively call on the rest.
    vector< vector<int> > subpermutations;

    // Get all combinations of elements of list where groupSize elements are chosen.
    vector< vector<int> > listCombinations = combinations( list, groupSize );

    // Iterate over listCombinations, the set of ways to choose groupSize indices to be contracted together.
    for ( vector< vector<int> >::iterator combination = listCombinations.begin(); combination != listCombinations.end(); ++combination ) {
        // Calculate difference between sets list and *combination to find remaining indices to use in recursive
        // call. Difference is computed using set_difference which is a member of the STL in header <algorithm>.
        difference = vector<int>( 2 * list.size() );
        iter = set_difference( list.begin(), list.end(), combination->begin(), combination->end(), difference.begin() );
        difference.resize( iter - difference.begin() );

        // Get index permutations for the remaining contractions.
        vector<int> subcontraction( ++contraction.begin(), contraction.end() );
        subpermutations = getIndexPermutations( subcontraction, difference );

        // Concatenate *combination with each subpermutation, and push the result on to indexPermutations.
        for ( vector< vector<int> >::iterator subpermutation = subpermutations.begin(); subpermutation != subpermutations.end(); ++subpermutation ) {
            vector<int> nextPermutation( *combination );
            nextPermutation.insert( nextPermutation.end(), subpermutation->begin(), subpermutation->end() );
            indexPermutations.push_back( nextPermutation );
        }
    }

    return indexPermutations;
}

vector< vector<int> > getIndexPermutations( vector<int> contraction ) {
    // Calculate the order of the contraction, n.
    unsigned int n = 0;
    for ( vector<int>::iterator groupSize = contraction.begin(); groupSize != contraction.end(); ++groupSize ) {
        n += *groupSize;
    }

    vector<int> list;
    for ( int i = 0; i < n; i++ ) {
        list.push_back( i );
    }

    return getIndexPermutations( contraction, list );
}

vector<TotalSignature> generateSignaturePermutations( vector< vector<int> > indexPermutations, TotalSignature signature ) {
    vector<TotalSignature> signatureSet;

    // Iterate over all permutations generated by getIndexPermutations.
    TotalSignature nextSignature;
    DeltaContractionSet deltaSet;
    DeltaContractionSet deltaBarSet;
    bool signatureIsNondegenerate;

    for ( vector< vector<int> >::iterator permutation = indexPermutations.begin(); permutation != indexPermutations.end(); ++permutation ) {
        deltaSet = DeltaContractionSet();
        for ( vector<IndexContraction>::iterator indexPair = signature.deltas.getIteratorBegin(); indexPair != signature.deltas.getIteratorEnd(); ++indexPair ) {
            deltaSet.addContraction( IndexContraction( permutation->at( indexPair->i ), permutation->at( indexPair->j ) ) );
        }

        deltaBarSet = DeltaContractionSet();
        for ( vector<IndexContraction>::iterator indexPair = signature.deltaBars.getIteratorBegin(); indexPair != signature.deltaBars.getIteratorEnd(); ++indexPair ) {
            deltaBarSet.addContraction( IndexContraction( permutation->at( indexPair->i ), permutation->at( indexPair->j ) ) );
        }

        nextSignature.deltas = deltaSet;
        nextSignature.deltaBars = deltaBarSet;

        signatureIsNondegenerate = true;
        for ( vector<TotalSignature>::iterator pushedSignature = signatureSet.begin(); pushedSignature != signatureSet.end(); ++pushedSignature ) {
            if ( nextSignature.areSignaturesDegenerate( *pushedSignature ) ) {
                signatureIsNondegenerate = false;
                break;
            }
        }

        if ( signatureIsNondegenerate ) signatureSet.push_back( nextSignature );
    }

    return signatureSet;
}

std::vector< std::vector<int> > calculateAllContractions( int n ) {
    vector< vector<int> > contractions;
    vector<int> nextContraction;
    vector< vector<int> > subcontractions;

    if ( n == 2 ) {
        nextContraction.push_back( 2 );
        contractions.push_back( nextContraction );
        return contractions;
    } else {
        nextContraction.push_back( n );
        contractions.push_back( nextContraction );

        subcontractions = calculateAllContractions( n - 2 );

        for ( vector< vector<int> >::iterator contraction = subcontractions.begin(); contraction != subcontractions.end(); ++contraction ) {
            nextContraction = vector<int>();
            nextContraction.push_back( 2 );
            nextContraction.insert( nextContraction.end(), contraction->begin(), contraction->end() );
            contractions.push_back( nextContraction );
        }
    }

    return contractions;
}

Sum generateCoordinateSpacePathIntegral( int n ) {
    Sum pathIntegral;

    vector< vector<int> > contractions = calculateAllContractions( n );
    TotalSignature signature;
    vector< vector<int> > indexPermutations;
    vector<TotalSignature> signaturePermutations;

    for ( vector< vector<int> >::iterator contraction = contractions.begin(); contraction != contractions.end(); ++contraction ) {
        Product nextPathIntegralTerm;
        for ( vector<int>::iterator contractedGroup = contraction->begin(); contractedGroup != contraction->end(); ++contractedGroup ) {
            nextPathIntegralTerm.addTerm( SymbolicTermPtr( Amaunet::SINE_PATH_INTEGRALS[ *contractedGroup ].copy() ) );
        }

        // Note: If you would like the signatures to presented such that the largest group appears first, reverse the
        // contraction vector here.

        Sum vertexIntegrals;
        signature = getDeltaSignature( *contraction );
        indexPermutations = getIndexPermutations( *contraction );
        signaturePermutations = generateSignaturePermutations( indexPermutations, signature );

        for ( vector<TotalSignature>::iterator permutation = signaturePermutations.begin(); permutation != signaturePermutations.end(); ++permutation ) {
            Product nextDeltaProduct;
            for ( vector<IndexContraction>::iterator indexPair = permutation->deltas.getIteratorBegin(); indexPair != permutation->deltas.getIteratorEnd(); ++indexPair ) {
                nextDeltaProduct.addTerm( SymbolicTermPtr( new Delta(  indexPair->i, indexPair->j ) ) );
            }

            for ( vector<IndexContraction>::iterator indexPair = permutation->deltaBars.getIteratorBegin(); indexPair != permutation->deltaBars.getIteratorEnd(); ++indexPair ) {
                Sum deltaBarSum;
                Product negativeDelta;

                negativeDelta.addTerm( SymbolicTermPtr( new CoefficientFloat( -1.0 ) ) );
                negativeDelta.addTerm( SymbolicTermPtr( new Delta( indexPair->i, indexPair->j ) ) );
                deltaBarSum.addTerm( SymbolicTermPtr( new CoefficientFloat( 1.0 ) ) );
                deltaBarSum.addTerm( negativeDelta.copy() );
                nextDeltaProduct.addTerm( deltaBarSum.copy() );
            }

            vertexIntegrals.addTerm( nextDeltaProduct.copy() );
        }

        nextPathIntegralTerm.addTerm( vertexIntegrals.copy() );
        pathIntegral.addTerm( nextPathIntegralTerm.copy() );
    }

    pathIntegral.reduceTree();
    return pathIntegral;
}

Sum pathIntegrateExpression( SymbolicTermPtr expr ) {
    Sum integratedExpression;

    if ( expr->getTermID() != TermTypes::SUM ) {
        return Sum();  // TODO: Raise exception.
    }

    SumPtr castExpr = static_pointer_cast<Sum>( expr );
    for ( vector<SymbolicTermPtr>::iterator term = castExpr->getIteratorBegin(); term != castExpr->getIteratorEnd(); ++term ) {
        if ( (*term)->getTermID() != TermTypes::PRODUCT ) {
            return Sum();  // TODO: Raise exception.
        }

        ProductPtr castTerm = static_pointer_cast<Product>( *term );
        Product integratedProduct;
        int orderInSigma = 0;
        vector<int> secondMatrixSIndices;

        for ( vector<SymbolicTermPtr>::iterator factor = castTerm->getIteratorBegin(); factor != castTerm->getIteratorEnd(); ++factor ) {
            if ( (*factor)->getTermID() == TermTypes::MATRIX_S ) {
                orderInSigma++;
                integratedProduct.addTerm( SymbolicTermPtr( new Delta( (*factor)->getIndices()[0], (*factor)->getIndices()[1] ) ) );
                secondMatrixSIndices.push_back( (*factor)->getIndices()[1] );
            } else {
                integratedProduct.addTerm( (*factor)->copy() );
            }
        }

        map<int, int> expressionToSignatureIndexMapping;
        for ( int i = 0; i < secondMatrixSIndices.size(); i++ ) {
            expressionToSignatureIndexMapping[ i ] = secondMatrixSIndices.at( i );
        }

        if ( orderInSigma > 1 ) {
            if ( orderInSigma % 2 == 0 ) {
                Sum spatialPathIntegral = generateCoordinateSpacePathIntegral( orderInSigma );
                spatialPathIntegral = spatialPathIntegral.getExpandedExpr();
                spatialPathIntegral.reduceTree();

                for ( vector<SymbolicTermPtr>::iterator pathIntegralTerm = spatialPathIntegral.getIteratorBegin(); pathIntegralTerm != spatialPathIntegral.getIteratorEnd(); ++pathIntegralTerm ) {
                    if ( (*pathIntegralTerm)->getTermID() != TermTypes::PRODUCT ) {
                        return Sum();  // TODO: Raise exception.
                    }

                    ProductPtr castPathIntegralTerm = static_pointer_cast<Product>( *pathIntegralTerm );

                    for ( vector<SymbolicTermPtr>::iterator pathIntegralFactor = castPathIntegralTerm->getIteratorBegin(); pathIntegralFactor != castPathIntegralTerm->getIteratorEnd(); ++pathIntegralFactor ) {
                        if ( (*pathIntegralFactor)->getTermID() == TermTypes::DELTA ) {
                            int* indices;  // Size of assigned array is 2.
                            indices = (*pathIntegralFactor)->getIndices();
                            indices[0] = expressionToSignatureIndexMapping[ indices[0] ];
                            indices[1] = expressionToSignatureIndexMapping[ indices[1] ];
                        }
                    }
                }

                integratedProduct.addTerm( spatialPathIntegral.copy() );
            } else {
                integratedProduct.addTerm( SymbolicTermPtr( new CoefficientFloat( 0.0 ) ) );
            }
        }

        integratedExpression.addTerm( integratedProduct.copy() );
    }

    return integratedExpression;
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

ostream& operator<<( std::ostream& os, const vector<DeltaContractionSet> &obj ) {
    os << "[";
    for( vector<DeltaContractionSet>::const_iterator vec = obj.begin(); vec != obj.end(); ++vec ) {
        os << " [ ";
        for ( vector<IndexContraction>::const_iterator element = (*vec).getIteratorBegin(); element != (*vec).getIteratorEnd(); ++element ) {
            os << "( " << element->i << ", " << element->j << " )";
        }
        os << " ] ";
    }
    os << "]";
}

ostream& operator<<( ostream& os, const vector<TotalSignature> &obj ) {
    os << "[";
    for ( vector<TotalSignature>::const_iterator vec = obj.begin(); vec != obj.end(); ++vec ) {
        os << " { " << vec->deltas << " | " << vec->deltaBars << " } ";
    }
    os << "]";
}