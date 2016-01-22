"""
Amaunet - High-order Lattice Perturbation Theory Expansion 
          for Interacting Quantum Matter
Weak-coupling Expansion for Fermionic Systems with Contact Interactions

Classes and Methods for Symbolic Perturbation Theory Manipulation

Andrew C. Loheac, Joaquin E. Drut
Department of Physics and Astronomy
University of North Carolina at Chapel Hill
"""

from copy import deepcopy
import gc
from itertools import combinations

# Static variables that define the behavior of symbolic calculations, system parameters,
# and output.
global NX, NTAU, GC_COLLECT_INTERVAL, num_memory_intensive_calls

NX = 2                      # Spatial volume of the system.
NTAU = 2                    # Temporal volume of the system.

SPLIT_SUMS_BY_LINE = True   # If true, each term of a sum will be output on a separate late.
GC_COLLECT_INTERVAL = 5000   # Number of calls to memory-intensive functions before garbage
                             # collection is explicitly called.

num_memory_intensive_calls = 0  # Local global copy for each process. Tracks the number of calls to
                                # memory intensive functions.
                                
# ****************************************************************************
#   EXCEPTION CLASSES
# ****************************************************************************

"""
General exception class that is to be raised when a user error would produce
an invalid symbolic calculation.
"""
class PTSymbolicException(Exception):
    """
    Constructor for PTSymbolicException.
    @param: value Error message to be passed out to stderr.
    """
    def __init__(self, value):
        self.value = value
        
# ****************************************************************************
#   SCALAR AND MATRIX/TENSOR TERM CLASSES
# ****************************************************************************
"""
Representation of the matrix M in the defined perturbation theory language.
Elements are in a basis such that they are indexed by two coupled (in space
and time) indices, where M has a square dimension of N_tau N_x^d.

The data structure for this matrix is optimized assuming that M is a sparse
matrix. Only non-zero elements are stored in memory.
"""
class MatrixM:
    def __init__( self, isInteracting=True, flavorLabel="" ):
        # Order of the derivative with respect to A for this instance of the matrix.
        # Must be an int.
        self.derivativeOrder = 0
        
        # Boolean flag indicating whether this matrix is interacting (e.g. whether
        # there is any dependence on A, or if A is set to zero).
        self.isInteracting = isInteracting
        
        # String indicating a unique identifier for each particle flavor.
        self.flavorLabel = flavorLabel
        
        
    """
    Pretty prints string representation of this matrix.
    """   
    def __str__( self ):
        if self.flavorLabel == "":
            if self.isInteracting:
                if self.derivativeOrder == 0:
                    return "M"
                elif self.derivativeOrder == 1:
                    return "dM / dA"
                else:
                    #return "d^" + str( self.derivativeOrder ) + " M / dA^" + str( self.derivativeOrder )
                    return "0.0"
            else:
                if self.derivativeOrder == 0:
                    return "M0"
                else:
                    return "dM0 / dA"
        else:
            if self.isInteracting:
                if self.derivativeOrder == 0:
                    return "M_" + str( self.flavorLabel )
                elif self.derivativeOrder == 1:
                    return "dM_" + str( self.flavorLabel ) + " / dA"
                else:
                    #return "d^" + str( self.derivativeOrder ) + " M_" + str( self.flavorLabel ) + " / dA^" + str( self.derivativeOrder )
                    return "0.0"
            else:
                if self.derivativeOrder == 0:
                    return "M0_" + str( self.flavorLabel )
                else:
                    return "dM0_" + str( self.flavorLabel ) + " / dA" 
    
    """
    Operator overload for equality comparison ('==').
    @param: other MatrixM object on the right hand side of operator.
    @return: True if the right hand side is equivalent to the left hand side.
    """
    def __eq__( self, other ):
        if isinstance( other, MatrixM ):
            return self.derivativeOrder == other.derivativeOrder \
                and self.flavorLabel == other.flavorLabel
        else:
            return False
    
    """
    Calculates the derivative of all matrix elements, and returns a new
    MatrixM object which is the derivative of this instance.
    @return: The derivative of this matrix.
    """
    def derivative( self ):
        D = MatrixM( self.isInteracting, self.flavorLabel )
        D.derivativeOrder = self.derivativeOrder + 1
        
        return D
    
    """
    Matrix M cannot be simplified -- do nothing.
    """
    def simplify( self ):
        pass
           
    """
    Sets this matrix as a non-interacting object (evaluate A to zero).
    """ 
    def setAsNoninteracting( self ):
        self.isInteracting = False
        
"""
Representation of the matrix B, or the inverse of the matrix M, in the defined
perturbation theory language.
"""
class MatrixB:
    def __init__( self, isInteracting=True, flavorLabel="" ):
        self.derivativeOrder = 0
        self.isInteracting = isInteracting
        self.flavorLabel = flavorLabel
        
    def __str__( self ):
        if self.flavorLabel == "":
            if self.isInteracting:
                if self.derivativeOrder == 0:
                    return "B"
                elif self.derivativeOrder == 1:
                    return "dB / dA"
                else:
                    return "d^" + str( self.derivativeOrder ) + " B / dA^" + str( self.derivativeOrder )
            else:
                if self.derivativeOrder == 0:
                    return "B0"
                else:
                    return "0.0"
        else:
            if self.isInteracting:
                if self.derivativeOrder == 0:
                    return "B_" + str( self.flavorLabel )
                elif self.derivativeOrder == 1:
                    return "dB_" + str( self.flavorLabel ) + " / dA"
                else:
                    return "d^" + str( self.derivativeOrder ) + " B_" + str( self.flavorLabel ) + " / dA^" + str( self.derivativeOrder )
            else:
                if self.derivativeOrder == 0:
                    return "B0_" + str( self.flavorLabel )
                else:
                    return "0.0"
    
    def __eq__( self, other ):
        if isinstance( other, MatrixB ):
            return self.derivativeOrder == other.derivativeOrder \
                and self.flavorLabel == other.flavorLabel
                
        else:
            return False
        
    def derivative( self ):
        if self.isInteracting:
            dM = MatrixM( self.isInteracting, self.flavorLabel )
            dM.derivativeOrder = self.derivativeOrder + 1
        
            return Product( [ CoefficientFloat( -1.0 ), MatrixB( self.isInteracting, self.flavorLabel ), dM, MatrixB( self.isInteracting, self.flavorLabel ) ] )
        
        else:
            return CoefficientFloat( 0.0 )
          
    def simplify( self ):
        pass
    
    def setAsNoninteracting( self ):
        self.isInteracting = False

class MatrixK:
    def __init__( self, flavorLabel = "" ):
        self.indices = None
        self.isFourierTransformed = False
        self.flavorLabel = flavorLabel
        
    def __str__( self ):
        if self.isFourierTransformed:
            return "D_" + str( self.flavorLabel ) + "_" + str( self.indices )
        else:
            if self.indices == None:
                return "K_" + str( self.flavorLabel )
            else:
                return "K_" + str( self.flavorLabel ) + "_"  + str( self.indices )
        
    def __eq__( self, other ):
        if isinstance( other, MatrixK ):
            return self.indices == other.indices and self.isFourierTransformed == other.isFourierTransformed and self.flavorLabel == other.flavorLabel
        else:
            return False
            
    def simplify( self ):
        pass
    
    def fourierTransform( self ):
        self.isFourierTransformed = True
        
class MatrixS:
    def __init__( self ):
        self.indices = None
        
    def __str__( self ):
        if self.indices == None:
            return "S"
        else:
            return "S_" + str( self.indices )
        
    def simplify(self):
        pass
"""
Representation of the determinant of M.
"""
class DetM:
    """
    Constructor for DetM.
    @param: M MatrixM which is the argument of the determinant.
    @param: flavorLabel Unique flavor label for the particle species of this determinant.
    @param: isInteracting Boolean flag indicating whether this determinant is dependent on A.
    @param: isInverted Boolean flag indicating whether this expression is inverted.
    """
    def __init__(self, M, isInteracting=True, isInverted=False ):
        self.M = M
        self.isInteracting = isInteracting
        self.isInverted = isInverted

    def __str__(self):
        if self.M.flavorLabel == "":
            if self.isInteracting:
                if self.isInverted:
                    return "1 / Det[ M ]"
                else:
                    return "Det[ M ]"
            else:
                if self.isInverted:
                    return "1 / Det[ M0 ]"
                else:
                    return "Det[ M0 ]"
        else:
            if self.isInteracting:
                if self.isInverted:
                    return "1 / Det[ M_" + str( self.M.flavorLabel ) + " ]"
                else:
                    return "Det[ M_" + str( self.M.flavorLabel ) + " ]"
            else:
                if self.isInverted:
                    return "1 / Det[ M0_" + str( self.M.flavorLabel ) + " ]"
                else:
                    return "Det[ M0_" + str( self.M.flavorLabel ) + " ]"
    
    def __eq__( self, other ):
        if isinstance( other, DetM):
            return self.M == other.M and self.isInverted == other.isInverted
        else:
            return False
    
    """
    Calculates the analytic derivative of this expression and returns the
    expression.
    @return: The derivative of this expression.
    """        
    def derivative( self ):
        if self.isInteracting:
            if self.isInverted:
                # Leave the below result in case it is useful for future implementations.
                # For now, this path should not be permitted; raise an exception.
                # return Product( [ CoefficientFloat( -1.0 ), Trace( Product( [ DetM( self.flavorLabel, True, True ), self.M.derivative() ] ) ), DetM( self.flavorLabel, True, True ) ] )
                raise PTSymbolicException( "The derivative of an inverted, interacting determinant should not be taken in this context." )
            else:
                return Product( [ DetM( self.M, True, False ), Trace( Product( [ MatrixB( True, self.M.flavorLabel ), self.M.derivative() ] ) ) ] )
        else:
            if self.isInverted:
                return "1 / Det[ M0 ]"
            else:
                return "Det[ M0 ]"
    
    """
    Mathematically simplify the expression. DetM cannot be further simplified,
    so do nothing (method is implemented for polymorphism).
    """        
    def simplify( self ):
        pass
    
    """
    Sets the determinant as a non-interacting object (independent of the coupling).
    """
    def setAsNoninteracting( self ):
        self.isInteracting = False
        
    
"""
Representation of a floating-point coefficient.
"""        
class CoefficientFloat:
    """
    Constructor for CoefficientFloat.
    @param: value Floating-point value of the coefficient.
    """
    def __init__( self, value ):
        self.value = value
        
    def __str__( self ):
        return str( float( self.value ) )
    
    """
    Returns the derivative of a constant, which is always zero.
    """
    def derivative( self ):
        return CoefficientFloat( 0.0 ) 
   
    """
    A floating-point number cannot be further simplified, so this
    method has no action.
    """
    def simplify( self ):
        pass
    
    def eval( self ):
        return float( self.value )
    
"""
Representation of a coefficient fraction, or ratio of two integers.
"""
class CoefficientFraction:
    """
    Constructor for CoefficientFraction.
    @param: num Value of the numerator.
    @param: den Value of the denominator.
    """
    def __init__( self, num, den ):
        self.num = num
        self.den = den
        
    def __str__( self ):
        return str( self.num ) + " / " + str( self.den )
    
    """
    Returns the derivative of a constant, which is always zero.
    """
    def derivative(self):
        return CoefficientFloat( 0.0 )
    
    """
    Simplifies the fraction.
    TODO: Make this method reduce the fraction to lowest terms.
    """
    def simplify( self ):
        pass
    
    """
    Evaluates the fraction and returns its floating-point value.
    """
    def eval( self ):
        return float( self.num ) / float( self.den )
 
"""
 Representation of the term A. This class is simply used to keep track of the
 order of a term in the perturbation series expansion.
 """
class TermA:
    """
    A is just A -- it cannot be simplified, so do nothing.
    """
    def simplify(self):
        pass
    
    def __str__(self):
        return "A"
    
    def __eq__( self, other ):
        if isinstance( other, TermA ):
            return True
        else:
            return False
                
# ****************************************************************************
#   OPERATOR AND EXPRESSION CLASSES
# ****************************************************************************
"""
Representation of a sum of an arbitrary number of terms.
"""
class Sum:
    """
    Constructor for Sum.
    @param: terms List of terms in the sum.
    """
    def __init__( self, terms ):
        self.terms = terms
        
    def __str__( self ):
        s = ""
        for i in range( 0, len( self.terms ) ):
            s += str( self.terms[ i ] )
            if not i == len( self.terms ) - 1:
                if SPLIT_SUMS_BY_LINE:
                    s += "\n\n"
                s += " + "
                
        return s
    
    def __eq__( self, other ):
        for term in self.terms:
            if not term in other.terms:
                return False
            
        return True
        
    """
    Calculate the derivative of a sum of terms.
    @returns: A sum that is the derivative of the calling sum.
    """
    def derivative( self ):
        self.reduceTree()
        D = deepcopy( self )
        for i in range( 0, len( self.terms ) ):
            D.terms[i] = self.terms[i].derivative()

        D.simplify()
        return D
    
    """
    Mathematically simplifies this expression. If any terms evaluate to zero,
    that term is removed from the expression.
    """ 
    def simplify( self ):     
        simplifiedTerms = []
        for term in self.terms:
            term.simplify()
            term = unpackTrivialExpr( term )
            if not ( str( term ).strip() == '0.0' or str( term ).strip() == '-0.0' or  isZeroTrace( term ) ):
                simplifiedTerms.append( term )
                
        if simplifiedTerms == []:
            simplifiedTerms = [CoefficientFloat( 0.0 )]
            
        self.terms = simplifiedTerms
        
        if len( self.terms ) == 0 :
            self.terms = [CoefficientFloat( 0.0 )]
    
    """
    Reduces the data structure representation of this instance to the most
    simple form. Note that some methods require this representation.
    """        
    def reduceTree( self ):    
        reducedTerms = []
        for s in self.terms:
            if isinstance( s, Sum ):
                s.reduceTree()
                for i in range( 0, len( s.terms ) ):
                    reducedTerms.append( unpackTrivialExpr( s.terms[ i ] ) )
            elif isinstance( s, Product ) or isinstance( s, Trace ):
                s.reduceTree()
                reducedTerms.append( unpackTrivialExpr( s ) )
            else:
                reducedTerms.append( unpackTrivialExpr( s ) )
                
        self.terms = unpackTrivialExpr( reducedTerms )
    
    """
    Mathematically expands the product of sums contained in this expression.
    @return: The expanded expression (a Sum).
    """        
    def getExpandedExpr( self ):
        expandedSum = Sum([])
        #self.reduceTree()
        for term in self.terms:
            if isinstance( term, Product ) or isinstance( term, Sum ):
                expandedSum.addTerm( term.getExpandedExpr() )
            else:
                expandedSum.addTerm( term )
                
        return expandedSum
          
    """
    Adds a term to the end of this Sum expression.
    @param: term Term to add to the end of this expression.
    """    
    def addTerm( self, term ):
        self.terms.append( term )
     
    """    
    Sets this expression as non-interacting (independent of the coupling). In
    the current context A is set to zero.
    """   
    def setAsNoninteracting( self ):
        for term in self.terms:
            if isinstance( term, Sum ) or isinstance( term, Product ) or isinstance( term, Trace ) or isinstance( term, DetM ) or isinstance( term, MatrixM ) or isinstance( term, MatrixB ): 
                term.setAsNoninteracting()
"""
Representation of a product of an arbitrary number of terms. All objects are
treated as noncommutative.
"""   
class Product:
    """
    Constructor for Product.
    @param: terms List of terms in the product.
    """
    def __init__( self, terms ): 
        self.terms = terms
        self.startingIndex = None
        self.endingIndex = None
    
    """
    @returns: Returns a string representation of the expression.
    """   
    def __str__( self ):
        s = " "
        for i in range( 0, len( self.terms ) ):
            s += "{" + str( self.terms[ i ] ) + "}"
            if not i == len( self.terms ) - 1:
                s += "  "
                
        return s
        
    """
    Calculate the derivative of a product of terms according to the product rule.
    @returns: A sum that is the derivative of the calling product.
    """
    def derivative( self ):
        self.reduceTree()
        D = Sum([])
        for i in range( 0, len( self.terms ) ):
            newTerm = Product( deepcopy( self.terms ) )
            newTerm.terms[ i ] = self.terms[ i ].derivative()
            D.terms.append( newTerm )
            
        D.simplify()
        return D
    
    """
    Mathematically simplifies this expression. If any term is zero, the entire
    expression is evaluated to zero. Any constants that evaluate to 1.0 are
    removed.
    """
    def simplify( self ):
        simplifiedTerms = []
        for term in self.terms:
            term.simplify()
            term = unpackTrivialExpr( term )
            if str( term ) == '0.0' or str( term ) == '-0.0' or isZeroTrace( term ):
                self.terms = [CoefficientFloat( 0.0 )]
                return
            elif not str( term ) == '1.0':
                simplifiedTerms.append( term )

        self.terms = simplifiedTerms
    
    """
    Reduces the data structure representation of this instance to the most
    simple form. Note that some methods require this representation.
    """    
    def reduceTree( self ):
        # TODO: Optimize necessary calls to unpackTrivialExpr().           
        reducedTerms = []
        for term in self.terms:
            if isinstance( term, Product ):
                term.reduceTree()
                for i in range( 0, len( term.terms ) ):
                    reducedTerms.append( unpackTrivialExpr( term.terms[ i ] ) )
            elif isinstance( term, Sum ) or isinstance( term, Trace ):
                term.reduceTree()
                reducedTerms.append( unpackTrivialExpr( term ) )    
            else:
                reducedTerms.append( unpackTrivialExpr( term ) )    
                
        self.terms = unpackTrivialExpr( reducedTerms )    
        
    """
    Recursively expand a product of sums to a explicit sum of products; i.e. the
    result expression is in the form AB... + CD... + ..., where applicable.
    @returns: The expanded product (a Sum) of this Product.
    """
    def getExpandedExpr( self ):
        checkGarbageCollection()      
        if len( self.terms ) == 0 or len( self.terms ) == 1:
            # No expansion to be done; return the original object.
            return self
        
        elif len( self.terms ) == 2:
            if isinstance( unpackTrivialExpr( self.terms[0] ), Sum ):  # Case 1: First term is a Sum, or both terms are a Sum.
                expandedSum = Sum([])
                factorA = self.terms[1]
                if isinstance( factorA, Product ) and factorA.containsSum():
                    factorA = factorA.getExpandedExpr()
                for i in range( 0, len( self.terms[0].terms ) ):
                    factorB = self.terms[0].terms[i]
                    if isinstance( factorB, Product ) and factorB.containsSum():
                        factorB = factorB.getExpandedExpr()
                    expandedSum.terms.append( Product( [ deepcopy( factorB ), deepcopy( factorA ) ] ).getExpandedExpr() )
                return expandedSum
            
            elif isinstance( unpackTrivialExpr( self.terms[1] ), Sum ):  # Case 2: Second term is a Sum.
                expandedSum = Sum([])
                factorA = self.terms[0]
                if isinstance( factorA, Product ) and factorA.containsSum():
                    factorA = factorA.getExpandedExpr()
                for i in range( 0, len( self.terms[1].terms ) ):
                    factorB = self.terms[1].terms[i]
                    if isinstance( factorB, Product ) and factorB.containsSum():
                        factorB = factorB.getExpandedExpr()
                    expandedSum.terms.append( Product( [ deepcopy( factorA ), deepcopy( factorB ) ] ).getExpandedExpr() )
                return expandedSum 
            else:  # Case 3: Neither term is a sum, so no expansion is necessary.
                return self
            
        else:  # Recursively expand the product.
            return Product( [ Product( self.terms[0:2] ).getExpandedExpr(), Product( self.terms[2:] ).getExpandedExpr() ] ).getExpandedExpr()
    
    """
    Helper function for determining if an instance of a Product contains a
    Sum as a term and therefore if expansion would be necessary. Calling this
    function prior to calling the getExpandedExpr() method on a product reduces
    unnecessary calls to the computationally expensive method.
    @param: product Product to be checked for a Sum.
    @return: True if expansion would be necessary for this Product, False otherwise.
    """
    def containsSum( self ):
        for term in self.terms:
            if isinstance( term, Sum ):
                return True
            
        return False
    
    """
    Adds a term to the Product expression.
    @param: term The term to add to the end of the expression.
    """
    def addTerm( self, term ):
        self.terms.append( term )
        
    """    
    Sets this expression as non-interacting (independent of the coupling). In
    the current context A is set to zero.
    """
    def setAsNoninteracting( self ):
        for term in self.terms:
            if isinstance( term, Sum ) or isinstance( term, Product ) or isinstance( term, Trace ) or isinstance( term, DetM ) or isinstance( term, MatrixM ) or isinstance( term, MatrixB ): 
                term.setAsNoninteracting()
        
"""
Representation of the trace of an expression.
"""
class Trace:
    """
    Constructor for Trace.
    @param: expr Expression that the trace is taken over.
    """
    def __init__( self, expr ):
        self.expr = expr
    
    def __str__( self ):
        return "Trace[ " + str( self.expr ) + " ]"
    
    """
    Calculates and returns the symbolic derivative with respect to the
    expansion parameter. Since the trace is a linear operator, this
    method simply takes the derivative of the argument.
    @return: The derivative of this expression.
    """
    def derivative( self ):
        return Trace( self.expr.derivative() )
    
    """
    Simplifies the expression within the trace.
    @return: The simplified expression.
    """
    def simplify( self ):
        self.expr.simplify()

    """
    Reduces the data structure representation of this instance to the most
    simple form. Note that some methods require this representation.
    """
    def reduceTree( self ):
        if isinstance( self.expr, Product ) or isinstance( self.expr, Sum ):
            self.expr.reduceTree()
    
    """    
    Sets this expression as non-interacting (independent of the coupling). In
    the current context A is set to zero.
    """        
    def setAsNoninteracting( self ):
        for term in self.expr.terms:
            if isinstance( term, Sum ) or isinstance( term, Product ) or isinstance( term, Trace ) or isinstance( term, DetM ) or isinstance( term, MatrixM ) or isinstance( term, MatrixB ): 
                term.setAsNoninteracting()
                
    def rewriteTraceInKSFormalism( self ):
        if not isinstance( self.expr, Product ):
            raise PTSymbolicException( "Argument of a Trace must be a Product before rewriting in KS formalism. Distribute trace operators first." )
        
        newExpr = Product([])
        for term in self.expr.terms:
            if isinstance( term, MatrixB ):
                newExpr.addTerm( MatrixK( term.flavorLabel ) )
            elif isinstance( term, MatrixM ):
                newExpr.addTerm( MatrixS() )
                
        self.expr = newExpr            
                
"""
Representation of a slightly modified definition of the Kronecker delta. A 
Delta object with two assigned indices behaves the same as \delta_{i,j}; an
object with more then two indices is "\delta 'bar'", which is written in
terms of a combination of standard Kronecker deltas. This class is intended
to be used in conjuction with an IndexedTrace.
"""            
class Delta:
    """
    Constructor for Delta.
    @param: indices Contracted indices for this Kronecker delta.
    """
    def __init__( self, indices, isSpatial=True, isTemporal=False, isBar = False ):
        self.indices = indices
        self.isBar = isBar
        self.isSpatial = isSpatial
        self.isTemporal = isTemporal
        
        if not len( indices ) == 2:
            raise PTSymbolicException( "Two indices must be specified in the constructor of Delta (" + str( len( indices ) ) + " was specified)." )
    
        if indices[1] > indices[0]:
            self.indices = ( indices[0], indices[1] )
        else:
            self.indices = ( indices[1], indices[0] )
            
    def __str__( self ):
        s = "Delta"
        if self.isBar:
            s += "Bar"
        s += "<"
        if self.isSpatial:
            s += "X"
        if self.isTemporal:
            s+= "T"
        
        s +=">( "
        for i in range( 0, len( self.indices ) - 1):
            s += str( self.indices[i] ) + ", "
        s += str( self.indices[-1] )
        s += " )"
        return s

    def __eq__( self, other ):
        if isinstance( other, Delta ):
            return self.indices == other.indices and self.isBar == other.isBar
        else:
            return False
        
    def simplify( self ):
        pass

class FourierSum:
    def __init__( self, indices, KOrder ):
        self.indices = indices
        self.KOrder = KOrder
        
    def __str__( self ):
        return "FourierSum<XT>" + str( self.indices )
    
    def __eq__( self, other ):
        if isinstance( other, FourierSum ):
            return str( set( self.indices ) ) == str( set( other.indices ) )
        else:
            return False
        
    def simplify( self ):
        pass
    
# ****************************************************************************
#   GENERIC HELPER FUNCTIONS
# ****************************************************************************

"""
Reduce a trivial expression to the child expression in the case where either
a Sum or a Product contains a single term.
@returns: The single child term of the expression, or None.
"""
def unpackTrivialExpr( expr ):
    if isinstance( expr, Sum ) or isinstance( expr, Product ):
        if len( expr.terms ) == 1:
            return expr.terms[0]
        else:
            return expr
    else:
        return expr
 
"""
Determines whether the expression within a trace trivially evaluates to zero.
@return: True if the expression trivially evaluates to zero, False otherwise.
"""
def isZeroTrace( expr ):
    if isinstance( expr, Trace ):
        if isinstance( expr.expr, Sum ) or isinstance( expr.expr, Product ):
            if len( expr.expr.terms ) == 0:
                return True
        
        if str( expr.expr ).strip() == "0.0" or str( expr.expr ).strip() == "0" or str( expr.expr ).strip() == "{0.0}" or str( expr.expr ).strip() == "{0}":
            return True
    
    return False

"""
Utility function to check if garbage collection should be explicitly called.
A global counter is used to track how many calls to predetermined memory
intensive functions are made by invoking this method. If garbage collection
is not yet to be run, this counter is simply incremented. Takes no arguments,
returns no value.
"""
def checkGarbageCollection():
    global num_memory_intensive_calls
    num_memory_intensive_calls += 1
    if num_memory_intensive_calls > GC_COLLECT_INTERVAL:
        print ">> [GC] Running garbage collection..."
        gc.collect()
        num_memory_intensive_calls = 0

"""
Helper function to get the total order of A in an expression where the product
of two factors are taken together. This method is called to determine if an
expression should be evaluated in full or not, considering that high-order
terms will be truncated at the end of the calculation anyway.
@param: prodA The first product to be evaluated (instance of Product)
@param: prodB The second product to be evaluated (instance of Product)
@return: An int indicating the order in A of the total expression (prodA * prodB).
"""
def getDualProductAOrder( prodA, prodB ):
    prodAOrder = 0
    if not isinstance( prodA, Product ):
        raise PTSymbolicException( "All expressions passed to getDualProductAOrder() must be a Product; prodA is not." )
    
    for term in prodA.terms:
        if isinstance( term, TermA ):
            prodAOrder += 1
            
    prodBOrder = 0
    if not isinstance( prodB, Product ):
        raise PTSymbolicException( "All expressions passed to getDualProductAOrder() must be a Product; prodB is not." )
    
    for term in prodB.terms:
        if isinstance( term, TermA ):
            prodBOrder += 1
            
    return prodAOrder + prodBOrder
  
def getAOrder( prod ):
    prodOrder = 0
    if not isinstance( prod, Product ):
        raise PTSymbolicException( "All expressions passed to getAOrder() must be a Product." )
    
    for term in prod.terms:
        if isinstance( term, TermA ):
            prodOrder += 1
            
    return prodOrder
      
# ****************************************************************************
#   EXPRESSION MANIPULATION FUNCTIONS
# ****************************************************************************

"""
Mathematically simplifies the expression of a trace. The trace operator is
distributed over a sum (e.g. Tr[ a + b + ... ] --> Tr[a] + Tr[b] + ...), 
and constant coefficients are also factored out of the trace. (e.g.
Tr[ {-1} a ] --> {-1} Tr[ a ]). Depending on the input, either a Sum
or Product object is returned, whichever is appropriate. In a typical
situation this method is not intended to be used outside of distributeCompleteExpr().
@return: A Sum or Product instance whose expression contains fully
         distributed traces.
"""
def distributeTrace( tr ):
    checkGarbageCollection()
    
    if not isinstance( tr, Trace ):
        raise PTSymbolicException( "An object other then a Trace was passed to distributeTrace()." )
    
    # Simple case: a single product inside a trace has been passed in.
    if isinstance( tr.expr, Product ):
        # Consider relatively more complicated cases, where the argument to the
        # trace is an expression like a ( b c + d e ); we must expand the product
        # and reduce the representation.
        expandedArgument = tr.expr  # Set scope; expression will be expanded next
                                    # if deemed necessary.
        if tr.expr.containsSum():
            expandedArgument = tr.expr.getExpandedExpr()
            expandedArgument.reduceTree()
        
        # If the expansion resulted in the product turning into a sum, evaluate
        # it appropriately through a recursive call.
        if isinstance( expandedArgument, Sum ):
            return distributeTrace( Trace( expandedArgument ) )
        else:  # Otherwise, if the expression is simple and a single product, proceed.
            simplifiedProduct = Product([])
            matrixTerms = []
            for term in expandedArgument.terms:
                if isinstance( term, CoefficientFloat ) or isinstance( term, CoefficientFraction ):
                    simplifiedProduct.terms.append( term )  # There used to be a deepcopy call here.
                elif isinstance( term, MatrixB ) or isinstance( term, MatrixM ):
                    matrixTerms.append( term )  # There used to be a deepcopy call here.
    
            simplifiedProduct.terms.append( Trace( unpackTrivialExpr( Product( matrixTerms ) ) ) )
            return simplifiedProduct
    
    # Typical case: a sum of terms inside a trace has been passed in to be evaluated.
    elif isinstance( tr.expr, Sum ):
        distributedSum = Sum([])
        for term in tr.expr.terms:
            distributedSum.addTerm( distributeTrace( Trace( term ) ) )  # There used to be a deepcopy call here.
            
        return distributedSum
    
    # Something trivial has been passed in...just spit it back!
    else:
        return tr

"""
Public method intended to be used as part of the integration and reduction
process. Performs the task of simplifying expressions containing traces;
the trace operator is distributed over sums, and constant coefficients are
pulled out of traces.
@param: expr 
""" 
def distributeAllTraces( expr ):
    # First, before doing anything, check if the expression can
    # be distributed as a whole, but do so only if necessary since
    # calling an expansion is computationally intensive:
    if isinstance( expr, Product ):
        if expr.containsSum():
            expr = expr.getExpandedExpr()
    
    # Reduce the expression tree.
    expr.reduceTree()
        
    distributedExpr = Sum([])
    if isinstance( expr, Sum ):
        for term in expr.terms:
            distributedExpr.addTerm( distributeAllTraces( term ) )
            
    elif isinstance( expr, Product ):
        distributedProduct = Product([])
        for factor in expr.terms:
            if isinstance( factor, Trace ):
                distributedTrace = distributeTrace( factor )
                distributedProduct.addTerm( distributedTrace )
            else:
                distributedProduct.addTerm( factor )
                
        if distributedProduct.containsSum():        
            distributedProduct = distributedProduct.getExpandedExpr()
            
        distributedExpr.addTerm( distributedProduct )
                
    elif isinstance( expr, Trace ):
        distributedExpr.addTerm( distributeTrace( expr ) )
           
    distributedExpr.reduceTree()                         
    return unpackTrivialExpr( distributedExpr )

"""
Truncate the order of A of an expansion.
@param: sum The expression (a Sum) to truncate.
@param: highestOrder The highest order in A which should appear in the truncated expression (an int).
@return: The truncated sum.
"""
def truncateAOrder( sum, highestOrder ):
    truncatedSum = Sum([])
    for term in sum.terms:
        AOrder = 0
        if isinstance( term, Sum  ) or isinstance( term, Product ):
            for factor in term.terms:
                if isinstance( factor, TermA ):
                    AOrder +=1
        
        if not AOrder > highestOrder:
            truncatedSum.terms.append( term )
            
    return truncatedSum

"""
Truncate terms in the expansion where A is of odd order. For contact
interactions, it has been shown that odd-ordered terms vanish, therefore
we eliminate such terms prior to computing the path integral to reduce
computational demands.
@param: sum Expression to truncate.
@return: Truncated expression.
"""
def truncateOddOrders( sum ):
    truncatedSum = Sum([])
    for term in sum.terms:
        AOrder = 0
        if isinstance( term, Sum  ) or isinstance( term, Product ):
            for factor in term.terms:
                if isinstance( factor, TermA ):
                    AOrder +=1
                    
        if AOrder % 2 == 0:
            truncatedSum.addTerm( term )
            
    return truncatedSum

def rewriteExprInKSFormalism( expr ):
    if not isinstance( expr, Sum ):
        raise PTSymbolicException( "Expression passed to rewriteExprInKSFormalism() must be an instance of a Sum." )

    for term in expr.terms:
        if isinstance( term, Product ):
            for factor in term.terms:
                if isinstance( factor, Trace ):
                    factor.rewriteTraceInKSFormalism()
        elif isinstance( term, Trace ):
            term.rewriteTraceInKSFormalism()
            
    return expr
             
def indexExpr( expr ):
    indexedExpr = Sum([])
    if isinstance( expr, Sum ):
        for term in expr.terms:
            if isinstance( term, Product ):
                if getAOrder( term ) == 2:
                    print "shtap!"
                nextIndex = 0
                indexedProduct = Product([])
                for factor in term.terms:
                    if isinstance( factor, Trace ):
                        for i in range( 0, len( factor.expr.terms ) - 1 ):
                            factor.expr.terms[i].indices = ( nextIndex + i, nextIndex + i + 1 )
                            indexedProduct.addTerm( factor.expr.terms[i] )
                        factor.expr.terms[-1].indices = ( nextIndex + i + 1, nextIndex )
                        indexedProduct.addTerm( factor.expr.terms[-1] )
                        nextIndex += i + 2
                    else:
                        indexedProduct.addTerm( factor )
                indexedExpr.addTerm( indexedProduct )
            else:
                indexedExpr.addTerm( term )
    else:
        raise PTSymbolicException( "indexExpr() expects a Sum as the passed expression." )
    
    return indexedExpr

def _getTerminatedContraction( index, deltas ):
    print (index, deltas[index])
    if index in deltas:
        if deltas[ index ] in deltas and not index == deltas[ index ]:
            return _getTerminatedContraction( deltas[ index ], deltas )
        else:
            return deltas[ index ]
    else:
        return index

def _constructContractionDict( contractions ):
    indexMapping = dict()
    for indexPair in contractions:
        if not indexPair[1] in indexMapping:
            if not indexPair[0] in indexMapping:
                indexMapping[ indexPair[1] ] = indexPair[0]
                indexMapping[ indexPair[0] ] = indexPair[0]
            else:
                if indexPair[0] < indexMapping[ indexPair[0] ]:
                    indexMapping[ indexPair[0] ] = indexPair[0]
                    indexMapping[ indexPair[1] ] = indexPair[0]
                else:
                    indexMapping[ indexPair[1] ] = indexMapping[ indexPair[0] ]
        else:
            if indexPair[0] < indexMapping[ indexPair[1] ]:
                indexMapping[ indexPair[0] ] = indexPair[0]
                indexMapping[ indexMapping[ indexPair[1] ] ] = indexPair[0]
                indexMapping[ indexPair[1] ] = indexPair[0]
                
                #if not indexPair[0] in indexMapping:
                    #indexMapping[ indexPair[0] ] = indexPair[0]
            elif not indexPair[0] in indexMapping:
                indexMapping[ indexPair[0] ] = indexMapping[ indexPair[1] ]
                
        if indexPair[1] in indexMapping.values():
                for index in indexMapping:
                    if indexMapping[ index ] > indexPair[0]:
                        indexMapping[ index ] = indexPair[0]
                        
    return indexMapping


def fourierTransformExpr( expr ):
    if not isinstance( expr, Sum ):
        raise PTSymbolicException( "Expression passed to fourierTransformExpr() must be an instance of a Sum." )

    transformedSum = Sum([])    
    for term in expr.terms:
        if not isinstance( term, Product ):
            raise PTSymbolicException( "Sum passed to fourierTransformExpr() must be fully distributed and reduced." )
        
        transformedProduct = Product([])
        KOrder = 0
        deltas = dict()
        fourierIndices = []
        indexPairsToBeContracted = []
        for factor in term.terms:
            if isinstance( factor, MatrixK ):
                factor.fourierTransform()
                #factor.indices = factor.indices[0]
                transformedProduct.addTerm( factor )
                fourierIndices.append( factor.indices )
                KOrder += 1
            elif isinstance( factor, Delta ) and factor.isBar == False:
                indexPairsToBeContracted.append( factor.indices )
                #if factor.indices[0] in deltas:
                #    deltas[ factor.indices[1] ] = deltas[ factor.indices[0] ]
                #else:
                #    deltas[ factor.indices[1] ] = factor.indices[0]
            else:
                transformedProduct.addTerm( factor )
            
        deltas = _constructContractionDict( indexPairsToBeContracted )
        if KOrder > 0:    
            #fourierIndices = range( 0, KOrder * 2, 2 )
            for contractedIndex in deltas:
                for i in range( 0, len( fourierIndices ) ):
                    if fourierIndices[i][0] == contractedIndex and not fourierIndices[i][1] == contractedIndex:
                        fourierIndices[i] = (deltas[ contractedIndex ], fourierIndices[i][1])
                    elif fourierIndices[i][1] == contractedIndex and not fourierIndices[i][0] == contractedIndex:
                        fourierIndices[i] = (fourierIndices[i][0], deltas[ contractedIndex ])
                    elif fourierIndices[i][1] == contractedIndex and fourierIndices[i][0] == contractedIndex:
                        fourierIndices[i] = (deltas[ contractedIndex ], deltas[ contractedIndex ])
                    # else: do nothing! The indices do not need to be contracted.
                        
                for term in transformedProduct.terms:
                    if isinstance( term, MatrixK ):
                        if term.indices[0] == contractedIndex and not term.indices[1] == contractedIndex:
                            term.indices = (deltas[ contractedIndex ], term.indices[1])
                        elif term.indices[1] == contractedIndex and not term.indices[0] == contractedIndex:
                            term.indices = (term.indices[0], deltas[ contractedIndex ])
                        elif term.indices[0] == contractedIndex and term.indices[1] == contractedIndex:
                            term.indices = (deltas[ contractedIndex ], deltas[ contractedIndex ])
                        # else: do nothing! The indices do not need to be contracted.
        
            transformedProduct.addTerm( FourierSum( fourierIndices, KOrder ) )
            
        transformedSum.addTerm( transformedProduct )
    
    return transformedSum    

def areTermsCommon( termA, termB ):
    #if isinstance( termA, CoefficientFloat ) or isinstance( termB, CoefficientFloat ):
    #    return False 
    if len( termA.terms ) == 1 or len( termB.terms ) == 1:
        return False
    
    determinantsA = dict()
    determinantsB = dict()
    AorderA = 0
    AorderB = 0
    DtermsA = dict()
    DtermsB = dict()
    fourierSumA = None
    fourierSumB = None
    
    for term in termA.terms:
        if isinstance( term, DetM ):
            if not term.M.flavorLabel in determinantsA:
                determinantsA[ term.M.flavorLabel ] = 1
            else:
                determinantsA[ term.M.flavorLabel ] += 1
        elif isinstance( term, TermA ):
            AorderA += 1
        elif isinstance( term, MatrixK ):
            if not term.flavorLabel in DtermsA:
                DtermsA[ term.flavorLabel ] = 1
            else:
                DtermsA[ term.flavorLabel ] += 1
        elif isinstance( term, FourierSum ):
            fourierSumA = term
            
    for term in termB.terms:
        if isinstance( term, DetM ):
            if not term.M.flavorLabel in determinantsB:
                determinantsB[ term.M.flavorLabel ] = 1
            else:
                determinantsB[ term.M.flavorLabel ] += 1
        elif isinstance( term, TermA ):
            AorderB += 1
        elif isinstance( term, MatrixK ):
            if not term.flavorLabel in DtermsB:
                DtermsB[ term.flavorLabel ] = 1
            else:
                DtermsB[ term.flavorLabel ] += 1
        elif isinstance( term, FourierSum ):
            fourierSumB = term

    if not AorderA == AorderB:
        return False
    
    if not fourierSumA == fourierSumB:
        return False
    
    if not DtermsA == DtermsB:
        return False
    
    if not determinantsA == determinantsB:
        return False
    
    return True                     
   
def combineLikeTerms( expr ):
    if not isinstance( expr, Sum ):
        raise PTSymbolicException( "Expression passed to combineLikeTerms() must be an instance of a Sum." )
    
    for i in range( 0, len( expr.terms ) ):
        runningLikeTermsCoefficient = 0
        
        termATotalCoefficient = 1
        combinedFactorA = Product([])
        for factor in expr.terms[i].terms:
            if isinstance( factor, CoefficientFloat ) or isinstance( factor, CoefficientFraction ):
                termATotalCoefficient *= factor.eval()
            else:
                combinedFactorA.addTerm( factor )
                
        for j in range( 0, len( expr.terms ) ):
            termBTotalCoefficient = 1
            if not i == j and areTermsCommon( expr.terms[i], expr.terms[j] ):       
                termBTotalCoefficient = 1
                for factor in expr.terms[j].terms:
                    if isinstance( factor, CoefficientFloat ) or isinstance( factor, CoefficientFraction ):
                        termBTotalCoefficient *= factor.eval()
                
                runningLikeTermsCoefficient += termBTotalCoefficient        
                expr.terms[j] = Product( [CoefficientFloat( 0.0 )] )
                
        combinedFactorA.addTerm( CoefficientFloat( round( termATotalCoefficient + runningLikeTermsCoefficient, 8 ) ) )
        expr.terms[i] = combinedFactorA
        
             
    return expr 

# ****************************************************************************
#   NUMERICAL LOOKUP TABLES
# ****************************************************************************

# Lookup tables required for calculations. All tables are implemented as dictionaries.

# Table for the result of the one-dimensional integral \int_{-\pi}^{\pi} \sin^n(x) dx,
# where n is the dictionary key, and the evaluated definite integral is the value. This
# set of integrals is used when evaluating path integrals for contact interactions.
SINE_PATH_INTEGRALS = { 1: CoefficientFloat( 0 ), 2: CoefficientFraction( 1, 2 ), 3: CoefficientFloat( 0 ), 4: CoefficientFraction( 3, 8 ), 5: CoefficientFloat( 0 ), 6: CoefficientFraction( 5, 16 ), 7: CoefficientFloat( 0 ), 8: CoefficientFraction( 35, 128 ), 9: CoefficientFloat( 0 ), 10: CoefficientFraction( 63, 256 ) }

# Table for the 1/n! coefficients that appear in a standard Taylor's series.
# Dictionary key refers to n, and the value refers to the evaluation of 1/n!.
# Having this table before execution avoids the need to compute n! for large
# values of n; however this also implies that that the program cannot compute
# expansions at arbitrary orders; although at that point you have much bigger
# problems to worry about.
TAYLORS_SERIES_COEFFICIENTS = { 0: CoefficientFloat( 1.0 ), 1: CoefficientFloat( 1.0 ), 2: CoefficientFraction( 1, 2 ), 3: CoefficientFraction( 1, 6 ), 4: CoefficientFraction( 1, 24 ), 5: CoefficientFraction( 1, 120 ), 6: CoefficientFraction( 1, 720 ), 7: CoefficientFraction( 1, 5040 ), 8: CoefficientFraction( 1, 40320 ), 9: CoefficientFraction( 1, 362880 ), 10: CoefficientFraction( 1, 3628800 ) }
   
# ****************************************************************************
#   PATH INTEGRATION AND ANALYTIC FOURIER TRANSFORM FUNCTIONS
# ****************************************************************************

def getDeltaSignature( contraction ):
    deltas = []
    deltaBars = []
    nextDeltaIndex = 0
    for i in range( 0, len( contraction ) ):
        for j in range( nextDeltaIndex, nextDeltaIndex + contraction[i] - 1 ):
            deltas.append( (j, j + 1) )
        nextDeltaIndex += 2
        
        if not i == len( contraction ) - 1:
            deltaBars.append( (j + 1, j + 2) )
        
    return [ deltas, deltaBars ]

def getIndexPermutations( signature, n ):
    if len( signature ) == 0:
        raise PTSymbolicException( "Signature length must be greater then zero." )
    
    indexPermutations = []
    for c in combinations( range( 0, n ), len( signature[1] ) * 2  ):
        newPermutation = [ -1 for i in range( 0, n ) ]
        deltaBarPositions = []
        i = 0
        for deltaBarIndices in signature[1]:
            newPermutation[ deltaBarIndices[0] ] = c[i]
            newPermutation[ deltaBarIndices[1] ] = c[i + 1]
            deltaBarPositions.append( deltaBarIndices[0] )
            deltaBarPositions.append( deltaBarIndices[1] )
            i += 2
            
        j = 0
        for i in range( 0, n ):
            if not i in deltaBarPositions:
                while j in newPermutation:
                    j += 1
                newPermutation[i] = j
        
        if newPermutation[0] < newPermutation[1]:        
            indexPermutations.append( newPermutation )
        
    return indexPermutations
    
def calculateAllContractions( n ):
    if not n % 2 == 0 or n <= 0:
        raise PTSymbolicException( "Parameter 'n' (index vector dimension) must be an even, positive, non-zero integer." )
    
    if n == 2:
        return [ (2,) ]
    else:
        listOfContractions = [ (n,) ]
        subContractions = calculateAllContractions( n - 2 )
        
        for contraction in subContractions:
            listOfContractions.append( (2,) + contraction )
            
        return listOfContractions
            
def generateCoordinateSpacePathIntegral( n ):
    listOfContractions = calculateAllContractions( n )
    
    integralResult = Sum([])
    for contraction in listOfContractions:
        integralResultTerm = Product([])
        
        # Append numerical coefficients that result from the one-dimensional
        # integral of powers of sine.
        for numberOfMatchingSigmas in contraction:
            integralResultTerm.addTerm( SINE_PATH_INTEGRALS[ numberOfMatchingSigmas ] )
        
        # Generate the sum of delta functions taking into account all applicable
        # permutations of indices.        
        deltaSum = Sum([])
        signature = getDeltaSignature( contraction )
        for p in getIndexPermutations( signature, n ):
            deltaProduct = Product([])
            
            # Add delta functions to the next product.
            for deltaIndices in signature[0]:
                deltaProduct.addTerm( Delta( ( p[ deltaIndices[0] ] * 2, p[ deltaIndices[1] ] * 2 ), True, True ) )
                
            # Add delta bar functions to the next product.
            for deltaBarIndices in signature[1]:
                deltaProduct.addTerm( Delta( ( p[ deltaBarIndices[0] ] * 2, p[ deltaBarIndices[1] ] * 2 ), True, True, True ) )
                
            deltaSum.addTerm( deltaProduct )
        
        integralResultTerm.addTerm( deltaSum )
        integralResult.addTerm( integralResultTerm )
    
    return integralResult

def pathIntegrateExpression( expr ):
    integratedSum = Sum([])
    if not isinstance( expr, Sum ):
        raise PTSymbolicException( "Expression passed to pathIntegrateExpression() must be an instance of a Sum." )
    
    for term in expr.terms:
        if not isinstance( term, Product ):
            raise PTSymbolicException( "Each term in the Sum passed to pathIntegrateExpression() must be an instance of a Product." )
        
        integratedProduct = Product([])
        sigmaOrder = 0
        for factor in term.terms:
            if isinstance( factor, MatrixS ):
                sigmaOrder +=1
                if factor.indices[1] > factor.indices[0]:
                    integratedProduct.addTerm( Delta( (factor.indices[0], factor.indices[1]) ) )
                else:
                    integratedProduct.addTerm( Delta( (factor.indices[1], factor.indices[0]) ) )
            else:
                integratedProduct.addTerm( factor )
        
        if sigmaOrder > 0:
            if sigmaOrder % 2 == 0:        
                integratedProduct.addTerm( generateCoordinateSpacePathIntegral( sigmaOrder ) )
            else:
                integratedProduct.addTerm( CoefficientFloat( 0.0 ) )
            
        integratedSum.addTerm( integratedProduct )
        
    return integratedSum