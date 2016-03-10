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
#include "PathIntegration.h"
#include <iostream>
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
    return i == rhs.i and j == rhs.j;
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

vector<DeltaContractionSet> generatePairedPermutations( vector<int> combination ) {
    assert( combination.size() % 2 == 0 );  // n must be even for contact interactions.

    vector<DeltaContractionSet> pairedPermutations;

    if ( combination.size() == 2 ) {  // Trivial base case.
        DeltaContractionSet permutation;

        // Impose convention that the larger valued index appears second in the contraction, e.g. ( 1, 2 ).
        if ( combination[0] < combination[1] ) {
            permutation.addContraction( IndexContraction( combination[0], combination[1] ) );
        } else {
            permutation.addContraction( IndexContraction( combination[1], combination[0] ) );
        }

        pairedPermutations.push_back( permutation );

    } else {
        vector< vector<int> > possiblePairs = combinations( combination, 2 );

        for ( vector< vector<int> >::iterator pair = possiblePairs.begin(); pair != possiblePairs.end(); ++pair ) {
            DeltaContractionSet permutation;

            // Impose convention that the larger valued index appears second in the contraction, e.g. ( 1, 2 ).
            if ( (*pair)[0] < (*pair)[1] ) {
                permutation.addContraction( IndexContraction( (*pair)[0], (*pair)[1] ) ); // Append selected pair combination
                                                                                     // to the next permutation.
            } else {
                permutation.addContraction( IndexContraction( (*pair)[1], (*pair)[0] ) );
            }

            vector<int> subcombination( combination );  // Copy combination from which we will remove the pair that was
                                                        // already added to the next permutation; make a recursive call
                                                        // that will generate the paired permutations for this subset
                                                        // (the 'tail' of the data structure).

            // Remove elements in *pair from subcombination. Choose to do each removal in separate loops, where we
            // break from the loop after the element is found and erased.
            for ( vector<int>::iterator subcombinationValue = subcombination.begin(); subcombinationValue != subcombination.end(); ++subcombinationValue ) {
                if ( (*subcombinationValue) == (*pair)[0] ) {
                    subcombination.erase( subcombinationValue );
                    break;
                }
            }

            for ( vector<int>::iterator subcombinationValue = subcombination.begin(); subcombinationValue != subcombination.end(); ++subcombinationValue ) {
                if ( (*subcombinationValue) == (*pair)[1] ) {
                    subcombination.erase( subcombinationValue );
                    break;
                }
            }

            vector<DeltaContractionSet> tailPairedPermutations = generatePairedPermutations( subcombination );

            // For each possible set of paired permutations for the tail, generate a new permutation, concatenate pairs,
            // and add the permutation to pairedPermutations, which is the data (a set of permutations) to be returned.
            for ( vector<DeltaContractionSet>::iterator tailPermutation = tailPairedPermutations.begin(); tailPermutation != tailPairedPermutations.end(); ++tailPermutation ) {
                DeltaContractionSet nextPermutation( permutation );  // Copy vector permutation.

                // Concatenate head and tail vectors.
                nextPermutation.addContractionSet( *tailPermutation );

                // Sort the IndexContraction elements of nextPermutation according to the operator< overload.
                sort( nextPermutation.getIteratorBegin(), nextPermutation.getIteratorEnd() );

                // Add nextPermutation to pairedPermutations only if pairedPermutations does not already contain an
                // identically equivalent object. Note that std::find is declared in <algorithm>. std::find is defined
                // such that here the iterator pairedPermutations.end() is returned if nextPermutation is not contained
                // in pairedPermutations.
                if ( find( pairedPermutations.begin(), pairedPermutations.end(), nextPermutation ) == pairedPermutations.end() ) {
                    pairedPermutations.push_back( nextPermutation );
                }
            }
        }
    }

    return pairedPermutations;
}

vector<TotalSignature> getIndexPermutations( TotalSignature signature, int n ) {
    // TODO: Add appropriate exception throws.

    vector<TotalSignature> indexPermutations;

    // Check for a trivial case - if the length of the signature.deltaBars vector is zero, the only possible
    // permutation is the passed signature. Return it. If the length of signature.deltaBars is non-zero, this
    // case will be included by the nature of finding all paired combinations.
    if ( signature.deltaBars.getNumContractions() == 0 ) {
        indexPermutations.push_back( signature );
        return indexPermutations;
    }


    // Generate a vector of integers ranging from 0 to n. Permutations of this list will be generated to provide all
    // possible combinations of indices for a given delta signature.
    vector<int> list;
    for ( int i = 0; i < n; i++ ) {
        list.push_back( i );
    }

    // Obtain all combinations of available indices, where we choose the appropriate number of \bar{\delta} indices
    // out of n indices. Recall that DeltaContractionSet::getNumContractions() gives the number of delta functions
    // in that set, therefore doubling that result gives the number of indices.
    vector< vector<int> > deltaBarIndexCombinations = combinations( list, 2 * signature.deltaBars.getNumContractions() );

    for ( vector< vector<int> >::iterator combination = deltaBarIndexCombinations.begin(); combination != deltaBarIndexCombinations.end(); ++combination ) {
        // Generate all possible sets of \bar{\delta} index pairs from the current combination.
        vector<DeltaContractionSet> pairedDeltaBarPermutations = generatePairedPermutations( *combination );

        // Iterate over all possible contraction sets of \bar{\delta} indices.
        for ( vector<DeltaContractionSet>::iterator contractionSet = pairedDeltaBarPermutations.begin(); contractionSet != pairedDeltaBarPermutations.end(); ++contractionSet ) {
            TotalSignature nextPermutation;
            nextPermutation.deltaBars = *contractionSet;  // Delta bar contractions already computed by
                                                          // generatePairedPermutations.
            nextPermutation.deltas = DeltaContractionSet();  // Delta contractions to be added below.

            map<int, int> signatureIndexMapping;

            // Generate the mapping from the signature index to the \bar{\delta} indices for this permutation. Loop
            // over IndexContraction objects of the \bar{\delta} signature.
            vector<IndexContraction>::iterator permutationIndexPair = contractionSet->getIteratorBegin();

            for( vector<IndexContraction>::iterator indexPair = signature.deltaBars.getIteratorBegin(); indexPair != signature.deltaBars.getIteratorEnd(); ++indexPair ) {
                assert( signature.deltaBars.getNumContractions() == contractionSet->getNumContractions() );

                if ( signatureIndexMapping.count( indexPair-> i ) == 0 ) {  // Index i is not in the map; add mapping.
                    signatureIndexMapping[ indexPair->i ] = permutationIndexPair->i;
                }

                if ( signatureIndexMapping.count( indexPair-> j ) == 0 ) {  // Index j is not in the map; add mapping.
                    signatureIndexMapping[ indexPair->j ] = permutationIndexPair->j;
                }

                ++permutationIndexPair;
            }

            // Generate the mapping from the signature index to the \delta indices for this permutation. Loop over
            // IndexContraction objects of the \delta signature.
            int currentFreeIndex = 0;
            for ( vector<IndexContraction>::iterator indexPair = signature.deltas.getIteratorBegin(); indexPair != signature.deltas.getIteratorEnd(); ++indexPair ) {
                // Advance the current free index to the next available index that is not a \bar{\delta} index.
                while ( contractionSet->containsIndex( currentFreeIndex ) ) {
                    currentFreeIndex++;
                }

                if ( signatureIndexMapping.count( indexPair-> i ) == 0 ) {  // Index i is not in the map; add mapping.
                    signatureIndexMapping[ indexPair->i ] = currentFreeIndex;
                    currentFreeIndex++;
                }

                // Advance the current free index to the next available index that is not a \bar{\delta} index.
                while ( contractionSet->containsIndex( currentFreeIndex ) ) {
                    currentFreeIndex++;
                }

                if ( signatureIndexMapping.count( indexPair-> j ) == 0 ) {  // Index j is not in the map; add mapping.
                    signatureIndexMapping[ indexPair->j ] = currentFreeIndex;
                    currentFreeIndex++;
                }

                assert( currentFreeIndex < n + 1 );
            }

            // Generate the DeltaContractionSet object for the \delta indices for this permutation. Loop over
            // IndexContraction objects of the \delta signature.
            for ( vector<IndexContraction>::iterator indexPair = signature.deltas.getIteratorBegin(); indexPair != signature.deltas.getIteratorEnd(); ++indexPair ) {
                nextPermutation.deltas.addContraction( IndexContraction( signatureIndexMapping[ indexPair->i ], signatureIndexMapping[ indexPair->j ] ) );
            }


            // Since indices of a delta function commute, adding every generated permutation will double the number of
            // terms in each contribution to the spatial path integral. To eliminate double counting, use the
            // convention that the permutation is added to the returned data structure only if for the first
            // IndexContraction, i < j.
            assert( nextPermutation.deltas.getNumContractions() > 0 );

            if ( nextPermutation.deltas.getIteratorBegin()->i < nextPermutation.deltas.getIteratorBegin()->j ) {
                indexPermutations.push_back( nextPermutation );
            }
        }
    }

    return indexPermutations;
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