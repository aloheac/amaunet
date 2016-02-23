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

#ifndef AMAUNETC_PATHINTEGRATION_Hs
#define AMAUNETC_PATHINTEGRATION_H

#include <map>
#include <vector>
#include "PTSymbolicObjects.h"

namespace Amaunet {
    extern std::map<int, CoefficientFraction> SINE_PATH_INTEGRALS;
}

void initializeStaticReferences();

std::vector<int*> getIndexPermutations( std::vector<int*> signature, int n );

std::vector<int*> calculateAllContractions( int n );

Sum generateCoordinateSpacePathIntegral( int n );

Sum pathIntegrateExpression( SymbolicTermPtr expr );

#endif //AMAUNETC_PATHINTEGRATION_H
