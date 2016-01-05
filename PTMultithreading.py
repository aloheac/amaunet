"""
Amaunet - High-order Perturbation Theory Expansion for Interacting Quantum Matter
Weak-coupling Expansion for Fermionic Systems with Contact Interactions

Multithreading Interface for High Performance Perturbation Theory Manipulation

Andrew C. Loheac, Joaquin E. Drut
Department of Physics and Astronomy
University of North Carolina at Chapel Hill
"""

import multiprocessing as mp
import PTSymbolicObjects as pt

global N_PROCS

N_PROCS = 8

"""
Basic implementation of an atomic counter for multithreaded function call.
A counter is incremented upon each call to tick().
"""
class atomicCounter:
    """
    Constructor for atomicCounter. Takes no arguments.
    """
    def __init__( self ):
        self.value = 0
        self.lock = mp.Lock()
        
    """
    Increments the counter by one unit.
    """
    def tick( self ):
        with self.lock:
            self.value += 1
    
    """
    Gets the current value of the counter. Requires main thread to possess the
    multithreading lock.
    @return: Value of this atomic counter.
    """
    def getValue( self ):
        with self.lock:
            return self.value
        
def _paralleleval_getExpandedExpr( args ):
    i = args[0]
    nTerms = args[1]
    product = args[2]
    sumOfProducts = args[3]
    print ">> [MT] Evaluating expansion for " + str( i ) + " of " + str( nTerms ) + "."

    expandedSum = pt.Sum([])
    if isinstance( product, pt.Product ):
        for term in sumOfProducts.terms:
            if not isinstance( term, pt.Product ) or pt.getDualProductAOrder( product, term ) <= N_PROCS:
                expandedSum.addTerm( pt.Product( [product, term] ).getExpandedExpr() )
    else:
        for term in sumOfProducts.terms:
            expandedSum.addTerm( pt.Product( [ product, term ] ).getExpandedExpr() )
    
    return expandedSum

def execParallelDualExpansion( expr1, expr2 ):
    print ">> [MT] Evaluating " + str( len(expr1.terms) * len(expr2.terms) ) + " terms in parallel."
    procs = mp.Pool( processes=N_PROCS )
    counter = atomicCounter()
    
    parallelExprs = []
    for i in range( 0, len( expr1.terms ) ):
        parallelExprs.append( [ i, len( expr1.terms ), expr1.terms[i], expr2 ] )
    
    expandedTerms = procs.map( _paralleleval_getExpandedExpr, parallelExprs )
    
    return pt.Sum( expandedTerms )

def _paralleleval_distributeAllTraces( args ):
    i = args[0]
    nTerms = args[1]
    term = args[2]
    print ">> [MT] Evaluating trace distribution and simplification for " + str( i ) + " of " + str( nTerms ) + "."
    
    return pt.distributeAllTraces( term )

def execParallelTraceDistribution( expr ):
    if not isinstance( expr, pt.Sum ):
        raise pt.PTSymbolicException( "A Sum instance must be passed to parallel distributeAllTraces()." )
    
    print ">> [MT] Evaluating trace distribution and simplification for " + str( len( expr.terms ) ) + "."
    procs = mp.Pool( processes=N_PROCS )
    
    parallelExprs = []
    for i in range( 0, len( expr.terms ) ):
        parallelExprs.append( [ i, len( expr.terms ), expr.terms[i] ] )
        
    distributedTerms = procs.map( _paralleleval_distributeAllTraces, parallelExprs )
    result = pt.Sum( distributedTerms )
    result.reduceTree()
    return result