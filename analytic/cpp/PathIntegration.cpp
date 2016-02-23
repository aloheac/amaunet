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

#include "PathIntegration.h"

using namespace std;

map<int, CoefficientFraction> Amaunet::SINE_PATH_INTEGRALS;  // Initialized by initializeStaticReferences().

void initializeStaticReferences() {
    using namespace Amaunet;

    // Initialize reference values for the result of the definite integral of various powers of sine:
    //
    //      \int_{ -\pi, \pi} dx \sin^n( x )
    //
    // where here the key of the dictionary is the integer n, and the value is the numerical result
    // of the integral expressed as a fraction (of type CoefficientFraction).
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
