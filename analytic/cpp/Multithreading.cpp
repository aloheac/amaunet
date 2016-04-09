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
#include "PathIntegration.h"

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

SymbolicTermPtr fullyEvaluatePartialExpression( SumPtr expr, int EXPANSION_ORDER_IN_A, int POOL_SIZE ) {
    cout << ">> Reducing expression tree..." << endl;
    expr->reduceTree();

    cout << ">> Truncating high orders in A of expansion..." << endl;
    expr = static_pointer_cast<Sum>( truncateAOrder( expr, EXPANSION_ORDER_IN_A ).copy() );

    cout << ">> Truncating odd orders in A of expansion..." << endl;
    expr = static_pointer_cast<Sum>( truncateOddOrders( expr ).copy() );

    cout << ">> Setting all objects in expansion as non-interacting..." << endl;
    expr->setAsNonInteracting();

    cout << ">> Rewriting expansion in terms of KS formalism..." << endl;
    rewriteSumInKSFormalism( expr );

    cout << ">> Indexing terms in expansion..." << endl;
    indexExpression( expr );

    cout << ">> Reducing expression tree..." << endl;
    expr->reduceTree();

    cout << ">> Computing path integral of expression..." << endl;
    Sum pathIntegral = pathIntegrateExpression( expr );

    cout << ">> Reducing expression tree..." << endl;
    pathIntegral.reduceTree();

    cout << ">> Performing trivial mathematical simplification..." << endl;
    pathIntegral.simplify();

    cout << ">> Expanding integrated expression..." << endl;
    pathIntegral = pathIntegral.getExpandedExpr();

    cout << ">> Reducing expression tree..." << endl;
    pathIntegral.reduceTree();

    cout << ">> Computing analytic Fourier transform of expression..." << endl;
    pathIntegral = fourierTransformExpression( pathIntegral.copy() );

    cout << ">> Reducing dummy indices of Fourier transformation..." << endl;
    pathIntegral.reduceFourierSumIndices();

    cout << ">> Combining like terms..." << endl;
    pathIntegral = combineLikeTerms( pathIntegral, POOL_SIZE );

    cout << ">> Performing trivial mathematical simplification..." << endl;
    pathIntegral.simplify();

    return pathIntegral.copy();
}

SumPtr fullyEvaluateExpressionByParts( SumPtr expr, int EXPANSION_ORDER_IN_A, int POOL_SIZE ) {
    if ( expr->getNumberOfTerms() <= POOL_SIZE ) {
        return static_pointer_cast<Sum>( fullyEvaluatePartialExpression( expr, EXPANSION_ORDER_IN_A, POOL_SIZE ) );
    } else {
        Sum evaluatedExpression;
        int i = 0;
        while ( i < expr->getNumberOfTerms() ) {
            cout << "Processing expansion for term range " << i << " to " << i + POOL_SIZE << " of " << expr->getNumberOfTerms() << " terms..." << endl;

            vector<SymbolicTermPtr>::iterator endingIterator;
            if ( i + POOL_SIZE > expr->getNumberOfTerms() ) {
                endingIterator = expr->getIteratorEnd();
            } else {
                endingIterator = expr->getIteratorBegin();
                endingIterator += i + POOL_SIZE;
            }

            vector<SymbolicTermPtr> nextGroupOfTerms( expr->getIteratorBegin() += i, endingIterator );
            SumPtr nextExpressionToEvaluate( new Sum( nextGroupOfTerms ) );

            evaluatedExpression.addTerm( fullyEvaluatePartialExpression( nextExpressionToEvaluate, EXPANSION_ORDER_IN_A, POOL_SIZE ) );

            i += POOL_SIZE;
        }

        evaluatedExpression.reduceTree();
        evaluatedExpression = combineLikeTerms( evaluatedExpression, POOL_SIZE );
        return static_pointer_cast<Sum>( evaluatedExpression.copy() );
    }
}

SumPtr expandAndEvaluateExpressionByParts( SumPtr exprA, SumPtr exprB, int EXPANSION_ORDER_IN_A, int POOL_SIZE ) {
    exprA->reduceTree();
    exprB->reduceTree();

    Sum expandedExpression;
    SymbolicTermPtr exprBCopy = exprB->copy();

    int termID = 1;
    for ( vector<SymbolicTermPtr>::iterator term = exprA->getIteratorBegin(); term != exprA->getIteratorEnd(); ++term ) {
        cout << ">> Performing expression expansion and evaluation for term " << termID << " of " << exprB->getNumberOfTerms() << "..." << endl;

        Product nextExpansion;
        nextExpansion.addTerm( (*term)->copy() );
        nextExpansion.addTerm( exprBCopy );

        SumPtr expanded = static_pointer_cast<Sum>( nextExpansion.getExpandedExpr().copy() );
        expanded->reduceTree();

        expandedExpression.addTerm( fullyEvaluateExpressionByParts( expanded, EXPANSION_ORDER_IN_A, POOL_SIZE ) );
        nextExpansion.clear();

        termID++;
    }

    cout << ">> Dual expansion complete. Reducing expression tree and combining like terms..." << endl;
    expandedExpression.reduceTree();
    expandedExpression = combineLikeTerms( expandedExpression, POOL_SIZE );
    return static_pointer_cast<Sum>( expandedExpression.copy() );
}