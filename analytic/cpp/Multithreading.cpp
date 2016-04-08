/* ***********************************************************************
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * High-order Perturbation Theory Analytics
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Multithreaded and Performance Function Wrappers for Symbolic Term
 * Manipulation Source Implementation
 *
 * v. 0.1		08 Apr 2016
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at Chapel Hill
 * ***********************************************************************
 */

#include <iostream>
#include "PTSymbolicObjects.h"

using namespace std;

Sum getDualExpansionByParts( SumPtr exprA, SumPtr exprB ) {
    exprA->reduceTree();
    exprB->reduceTree();

    Sum expandedExpression;
    SymbolicTermPtr exprBCopy = exprB->copy();

    int termID = 1;
    for ( vector<SymbolicTermPtr>::iterator term = exprA->getIteratorBegin(); term != exprA->getIteratorEnd(); ++term ) {
        cout << ">> Performing dual expression expansion for term " << termID << " of " << exprB->getNumberOfTerms() << "..." << endl;

        Product nextExpansion;
        nextExpansion.addTerm( (*term)->copy() );
        nextExpansion.addTerm( exprBCopy );

        SumPtr expanded = static_pointer_cast<Sum>( nextExpansion.getExpandedExpr().copy() );
        expanded->reduceTree();

        expandedExpression.addTerm( static_pointer_cast<SymbolicTerm>( expanded ) );
        nextExpansion.clear();

        termID++;
    }

    cout << ">> Dual expansion complete. Reducing expression tree..." << endl;
    expandedExpression.reduceTree();

    return expandedExpression;
}
