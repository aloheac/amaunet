/* ***********************************************************************
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * High-order Perturbation Theory Analytics
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Symbolic Path Integration and Fourier Transform Functions for
 * Contact Interactions Header
 *
 * v. 0.1		22 Feb 2016
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at Chapel Hill
 * ***********************************************************************
 */

#ifndef AMAUNETC_PATHINTEGRATION_H
#define AMAUNETC_PATHINTEGRATION_H

#include <map>
#include <vector>
#include "PTSymbolicObjects.h"

/*
 * ***********************************************************************
 * CLASS AND STRUCT DEFINITIONS
 * ***********************************************************************
 */

struct IndexContraction{

    IndexContraction();

    IndexContraction( int a, int b );

    int i;

    int j;

    bool operator<( const IndexContraction& rhs ) const;

    bool operator>( const IndexContraction& rhs ) const;

    bool operator==( const IndexContraction& rhs ) const;

    bool containsIndex( int index );

};

class DeltaContractionSet {

public:

    DeltaContractionSet();

    void addContraction( IndexContraction newContraction );

    void addContractionSet( DeltaContractionSet set );

    unsigned int getNumContractions();

    bool containsIndex( int index );

    std::vector<IndexContraction>::iterator getIteratorBegin();

    std::vector<IndexContraction>::iterator getIteratorEnd();

    std::vector<IndexContraction>::const_iterator getIteratorBegin() const;

    std::vector<IndexContraction>::const_iterator getIteratorEnd() const;

    std::string to_string() const;

    bool operator==( const DeltaContractionSet &rhs ) const;

private:

    std::vector<IndexContraction> contractions;

};

struct TotalSignature {

    TotalSignature();

    DeltaContractionSet deltas;

    DeltaContractionSet deltaBars;

    bool areSignaturesDegenerate( const TotalSignature &signature );

    bool isValidSignature();

};

class DeltaSignature {

public:

    DeltaSignature();

    void addContractionSet( DeltaContractionSet contractionSet );

    std::string to_string() const;

private:

    std::vector<DeltaContractionSet> signatureContractions;

};

/*
 * ***********************************************************************
 * STATIC VARIABLES
 * ***********************************************************************
 */

namespace Amaunet {
    extern std::map<int, CoefficientFraction> SINE_PATH_INTEGRALS;
}

/*
 * ***********************************************************************
 * FUNCTION DECLARATIONS
 * ***********************************************************************
 */

void initializeStaticReferences();

TotalSignature getDeltaSignature( std::vector<int> contraction );

std::vector< std::vector<int> > getIndexPermutations( std::vector<int> contraction, std::vector<int> list );

std::vector< std::vector<int> > getIndexPermutations( std::vector<int> contraction );

std::vector<TotalSignature> generateSignaturePermutations( std::vector< std::vector<int> > indexPermutations, TotalSignature signature );

std::vector< std::vector<int> > calculateAllContractions( int n );

Sum generateCoordinateSpacePathIntegral( int n );

Sum pathIntegrateExpression( SymbolicTermPtr expr );

std::vector< std::vector<int> > combinations( std::vector<int> list, int k );

/*
 * ***********************************************************************
 * INPUT REDIRECTION OPERATOR OVERLOADS
 * ***********************************************************************
 */

std::ostream& operator<<( std::ostream& os, const DeltaContractionSet &obj );

std::ostream& operator<<( std::ostream& os, const DeltaSignature &obj );

std::ostream& operator<<( std::ostream& os, const std::vector< std::vector<int> > &obj );

std::ostream& operator<<( std::ostream& os, const std::vector<DeltaContractionSet> &obj );

std::ostream& operator<<( std::ostream& os, const std::vector<TotalSignature> &obj );

#endif //AMAUNETC_PATHINTEGRATION_H
