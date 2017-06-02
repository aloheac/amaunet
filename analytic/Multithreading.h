/* ***********************************************************************
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * High-order Perturbation Theory Analytics
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Multithreaded and Performance Function Wrappers for Symbolic Term
 * Manipulation Headers
 *
 * v. 0.1		08 Apr 2016
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at Chapel Hill
 * ***********************************************************************
 */

#ifndef AMAUNETC_MULTITHREADING_H
#define AMAUNETC_MULTITHREADING_H

#include "PTSymbolicObjects.h"

Sum getDualExpansionByParts( SumPtr exprA, SumPtr exprB );

SymbolicTermPtr fullyEvaluatePartialExpression( SumPtr expr, int EXPANSION_ORDER_IN_A, int POOL_SIZE );

SumPtr fullyEvaluateExpressionByParts( SumPtr expr, int EXPANSION_ORDER_IN_A, int POOL_SIZE );

SumPtr expandAndEvaluateExpressionByParts( SumPtr exprA, SumPtr exprB, int EXPANSION_ORDER_IN_A, int POOL_SIZE );

SumPtr multithreaded_expandAndEvaluateExpressionByParts( SumPtr exprA, SumPtr exprB, int EXPANSION_ORDER_IN_A, int POOL_SIZE, int NUM_THREADS );

SumPtr multithreaded_getDualExpansionByParts( SumPtr exprA, SumPtr exprB, int NUM_THREADS );

#endif //AMAUNETC_MULTITHREADING_H
