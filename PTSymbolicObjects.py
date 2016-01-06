"""
Amaunet - High-order Perturbation Theory Expansion for Interacting Quantum Matter
Weak-coupling Expansion for Fermionic Systems with Contact Interactions

Classes and Methods for Symbolic Perturbation Theory Manipulation

Andrew C. Loheac, Joaquin E. Drut
Department of Physics and Astronomy
University of North Carolina at Chapel Hill
"""

from copy import deepcopy
import gc

# Static variables that define the behavior of symbolic calculations, system parameters,
# and output.
global NX, NTAU, MINIMAL_U_DISPLAY, MINIMAL_M_DISPLAY, MINIMAL_GAMMA_DISPLAY, GC_COLLECT_INTERVAL, num_memory_intensive_calls

NX = 2                      # Spatial volume of the system.
NTAU = 2                    # Temporal volume of the system.

MINIMAL_U_DISPLAY = False      # If true, the matrix U will be output using a shorthand notation.
MINIMAL_M_DISPLAY = True       # If true, the matrix M will be output using a shorthand notation.
MINIMAL_GAMMA_DISPLAY = False  # If true, Gamma expressions will be output using a shorthand notation.
SPLIT_SUMS_BY_LINE = False     # If true, each term of a sum will be output on a separate late.
GC_COLLECT_INTERVAL = 1000     # Number of calls to memory-intensive functions before garbage
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
    def __init__( self, Nx, Ntau, spatialDimension, isInteracting=True, flavorLabel="", constructMatrix=False ):
        # Dictionary that holds non-zero matrix elements. The key is given as a
        # two element tuple with temporal indices ( i, j ). The value may be any
        # valid expression.
        self.elements = dict()
        
        # Spatial length independent of the volume. Note that the spatial
        # dimension is specified separately; Nx**spatialDimension gives
        # the total spatial volume. Must be an int.
        self.Nx = Nx
        
        # Temporal volume. Must be an int.
        self.Ntau = Ntau
        
        # Spatial dimension; Nx**spatialDimension gives the total spatial volume.
        # Must be an int.
        self.spatialDimension = spatialDimension

        # Order of the derivative with respect to A for this instance of the matrix.
        # Must be an int.
        self.derivativeOrder = 0
        
        # Boolean flag indicating whether this matrix is interacting (e.g. whether
        # there is any dependence on A, or if A is set to zero).
        self.isInteracting = isInteracting
        
        # String indicating a unique identifier for each particle flavor.
        self.flavorLabel = flavorLabel
        
        # Boolean indicating if the matrix representation was explicitly
        # constructed. For some steps of the calculation this is not
        # necessary, so resource utilization is reduced somewhat.
        self.matrixIsConstructed = constructMatrix
        
        if constructMatrix:
            # Place identity elements along the matrix diagonal.
            for i in range( 0, self.Ntau):
                self.setElement(i, i, CoefficientFloat( 1.0 ) )
                
            # Place U matrices along lower triangular off diagonal, indexed
            # explicitly by a temporal index. The spatial index is left
            # unsummed, as denoted by the index "S$".
            for i in range( 0, self.Ntau - 1 ):
                self.setElement( i + 1, i, Product( [ CoefficientFloat( -1.0 ), MatrixU( "$",  i + 1, Nx, spatialDimension ) ] ) )
                
            # Place final U matrix at upper right hand corner element.
            self.setElement( 0, self.Ntau - 1, MatrixU("$", self.Ntau, Nx, spatialDimension ) )
        
    """
    Pretty prints string representation of this matrix.
    """   
    def __str__( self ):
        if not MINIMAL_M_DISPLAY and self.matrixIsConstructed:
            columnWidth = 60        
            elementValues = []
            s = ""
            for i in range( 0, self.Ntau ):
                for j in range( 0, self.Ntau ):
                    s += "{:>" + str( columnWidth ) + "}"
                    if ( i, j ) in self.elements:
                        elementValues.append( str( self.elements[ ( i, j ) ] ) )   
                    else:
                        elementValues.append( '0.0' )
                    s += '\t'
                s += '\n'
            
            return s.format( *elementValues )
        else:
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
                        return "M0_" + str( self.flavorLabel )
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
    Gets the matrix elements at temporal indices ( i, j ). Returns zero if the
    indices are not found in the internal dictionary.
    @param: i First temporal index (row) of the matrix element.
    @param: j Second temporal index (column) of the matrix element.
    @return: The expression of the matrix element at indices ( i, j ).
    """            
    def getElement( self, i, j ):
        if not self.matrixIsConstructed:
            raise PTSymbolicException( "getElement() was called on an instance of MatrixM whose matrix was not constructed." )
        
        if ( i, j ) in self.elements:
            return self.elements[ ( i, j ) ]
        else:
            return CoefficientFloat( 0.0 )
    
    """
    Sets the matrix element at temporal indices ( i, j ).
    @param i First temporal index (row) of the matrix element to set.
    @param j Second temporal index (column) of the matrix element to set.
    @param entry Expression of the matrix element to set.
    """
    def setElement(self, i, j, entry ):
        if not self.matrixIsConstructed:
            raise PTSymbolicException( "setElement() was called on an instance of MatrixM whose matrix was not constructed." )
        
        self.elements[ ( i, j ) ] = entry
    
    """
    Calculates the derivative of all matrix elements, and returns a new
    MatrixM object which is the derivative of this instance.
    @return: The derivative of this matrix.
    """
    def derivative( self ):
        D = MatrixM( self.Nx, self.Ntau, self.spatialDimension, self.isInteracting, self.flavorLabel, self.matrixIsConstructed )
        
        if self.matrixIsConstructed:
            for key in self.elements:
                D.setElement( key[0], key[1], self.elements[ key ].derivative() )
            D.simplify()
            
        D.derivativeOrder = self.derivativeOrder + 1
        
        return D
    
    """
    Mathematically simplifies the expressions for all matrix elements.
    """
    def simplify( self ):
        if self.matrixIsConstructed:
            for key in self.elements:
                self.elements[ key ].simplify()
            
    def setAsNoninteracting( self ):
        self.isInteracting = False
        
"""
Representation of the matrix B, or the inverse of the matrix M, in the defined
perturbation theory language.
"""
class MatrixB:
    def __init__( self, Nx, Ntau, spatialDimension, isInteracting=True, flavorLabel="" ):
        self.elements = dict()
        self.Nx = Nx
        self.Ntau = Ntau
        self.spatialDimension = spatialDimension
        self.isInteracting = isInteracting
        self.derivativeOrder = 0
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
        
    def derivative( self ):
        dM = MatrixM( self.Nx, self.Ntau, self.spatialDimension, self.isInteracting, self.flavorLabel )
        dM.derivativeOrder = self.derivativeOrder + 1
        
        return Product( [ CoefficientFloat( -1.0 ), MatrixB( self.Nx, self.Ntau, self.spatialDimension, self.isInteracting, self.flavorLabel ), dM, MatrixB( self.Nx, self.Ntau, self.spatialDimension, self.isInteracting, self.flavorLabel ) ] )
          
    def simplify( self ):
        pass
    
    def setAsNoninteracting( self ):
        self.isInteracting = False
    
"""
Representation of the matrix U in the defined perturbation theory language.
"""
class MatrixU:
    """
    Constructor for MatrixU.
    @param: spatialIndex Spatial index of this matrix.
    @param: temporalIndex Temporal index of this matrix.
    @param: Nx Spatial lattice size of this matrix.
    @param: spatialDimension Spatial dimension of this matrix, such that the
                             total spatial volume is Nx**d.
    """
    def __init__( self, spatialIndex, temporalIndex, Nx, spatialDimension ):
        self.derivativeOrder = 0
        self.spatialIndex = spatialIndex
        self.temporalIndex = temporalIndex
        self.Nx = Nx
        self.spatialDimension = spatialDimension
     
    def __str__( self ):
        # The reported expression depends on the derivative order. Implemented
        # for contact interactions.
        if not MINIMAL_U_DISPLAY:
            if self.derivativeOrder == 0:
                return "Exp(-T/2) (1+A Sin[Sigma_{ T" + str( self.temporalIndex ) + " }( S" + str( self.spatialIndex ) + " )]) Exp(-T/2)"
            elif self.derivativeOrder == 1:
                return "Exp(-T/2) Sin[Sigma_{ T" + str( self.temporalIndex ) + " }( S" + str( self.spatialIndex ) + " )] Exp(-T/2)"
            elif self.derivativeOrder >= 2:
                return "0.0"
            else:
                raise PTSymbolicException( "Invalid derivative order specifier in MatrixU with indices ( " + str( self.spatialIndex ) + str( self.temporalIndex ) + " )." )
        else:
            if self.derivativeOrder < 2:
                return "U^(" + str( self.derivativeOrder ) + ")_{ T" + str( self.temporalIndex ) + " }( S" + str( self.spatialIndex ) + " )"
            else:
                return "0.0"

    """
    Evaluates the derivative of U, which in practice increments the derivative
    order counter. The reported expression depends on the specific interaction.
    @return: self with the derivative counter incremented.
    """
    def derivative( self ):
        D = MatrixU( self.spatialIndex, self.temporalIndex, self.Nx, self.spatialDimension )
        D.derivativeOrder = self.derivativeOrder + 1
        return D
    
    """
    Mathematically simplify this expression. MatrixU cannot be further
    simplified, so do nothing.
    """
    def simplify( self ):
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
                return Product( [ DetM( self.M, True, False ), Trace( Product( [ MatrixB( self.M.Nx, self.M.Ntau, self.M.spatialDimension, True, self.M.flavorLabel ), self.M.derivative() ] ) ) ] )
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
            if not ( str( term ).strip() == '0.0' or  isZeroTrace( term ) ):
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
            if str( term ) == '0.0' or isZeroTrace( term ):
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
                
"""
Representation of an expression where all terms within a trace are assigned a
set of indices. Note that the expression that is passed to the constructor
must be a single product; the complete expression being evaluated must
be fully expanded prior to using this object. At present this assertion
is not enforced.
"""               
class IndexedTrace:
    """
    Constructor for IndexedTrace.
    @param: expr Expression that is the argument to the trace. Must be a single product.
    @param: startingIndex Starting integer for indexing elements of the trace, which is
                          by definition cyclic. Default value is 0.
    """
    def __init__( self, expr, startingIndex=0 ):
        self.expr = expr
        self.indices = []
        self.startingIndex = startingIndex
        self.endingIndex = None  # Set scope.
        
        # Index all terms within traces. Note that the expression passed to the
        # constructor is expected to be already distributed, in the sense that
        # all constant coefficients are pulled out of the trace, and the trace
        # operator has been linearly distributed over all terms in a sum.
        self._getExprIndices()
    
    """
    Private method which indexes all elements within a trace. In principle this method may
    be called when updating the argument of the trace, although this behavior is not
    preferred.
    """    
    def _getExprIndices( self ):
        # Check for trivial cases and set the appropriate result.
        if len( self.expr.terms ) == 0:
            self.indices = []
            self.endingIndex = None
        elif len( self.expr.terms ) == 1:
            self.indices = [ (self.startingIndex, self.startingIndex) ]
            self.endingIndex = self.startingIndex
        else:  # For non-trivial products, calculate the connected indices.
            indices = []
            for i in range( 0, len( self.expr.terms ) - 1):
                indices.append( (self.startingIndex + i, self.startingIndex + i + 1) )
            indices.append( (self.startingIndex + i + 1, self.startingIndex) )
            
            self.indices = indices
            self.endingIndex = self.startingIndex + i + 1
         
    def __str__( self ):
        s = "Trace[ "
        for i in range( 0, len( self.expr.terms ) ):
            s += "{ " + str( self.expr.terms[i] ) + "_" + str( self.indices[i] ) + " }"
        s += " ]"
        return s

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
    def __init__( self, indices ):
        self.indices = indices
        self.isBar = None
        
        if len( indices ) <= 1:
            raise PTSymbolicException( "Two or more indices must be specified in the constructor of Delta (" + str( len( indices ) ) + " was specified)." )
        elif len( indices ) == 2:
            self.isBar = False
        else:
            self.isBar = True
    
    def __str__( self ):
        s = "Delta( "
        for index in self.indices:
            s += str( index ) + ", "
        s += " )"
        return s
            
    def getAllPairings( self ):
        pairings = []
        for i in range( 0, len( self.indices ) ):
            for j in range( i + 1, len( self.indices ) ):
                pairings.append( (self.indices[i], self.indices[j]) )
        return pairings
 
class IndexedProduct( Product ):
    def __init__( self, terms, indices ):
        Product.__init__( self, terms )  # Call constructor of parent class Product.
        self.indices = indices
        
    """
    @returns: Returns a string representation of the expression.
    """   
    def __str__( self ):
        s = " "
        for i in range( 0, len( self.terms ) ):
            s += "{" + str( self.terms[ i ] )
            
            if not self.indices[ i ] == None:
                s += "_" + str( self.indices[ i ] )
            
            s += "}"
            
            if not i == len( self.terms ) - 1:
                s += "  "
                
        return s
    
    def addTerm( self, term, index ):
        self.terms.append( term )
        self.indices.append( index )

class MatrixT:
    def __init__( self, spatialIndices ):
        self.spatialIndices = spatialIndices
        
    def __str__( self ):
        return "T_" + str( self.spatialIndices )
    
"""
Representation of the mathematical Gamma object that appears in perturbation
theory under the current context. Indices are managed separately by this object;
the expression in this container should be a Product.
"""
class Gamma:
    """
    Constructor for Gamma. Note that this object does not construct any indices
    for the contained expressions; the list of indices corresponding to each
    factor should be generated elsewhere (as required by the corresponding math).
    @param: expr A non-indexed Product that is an argument to this object.
    @param: indices A list of two-element tuples that correspond to the indices
                    of each factor in the contained Product.
    """
    def __init__( self, expr, indices ):
        self.expr = expr  # A non-indexed Product.
        self.indices = indices
        self.uniqueIndices = []
        self.integratedTensor = None  # Tensor must be constructed later with a
                                      # call to constructIntegratedTensor().
        self.gammaOrder = str( len( self.expr.terms) )
              
    def __str__( self ):
        s = "Gamma<" + str( len( self.expr.terms ) ) + ">_" + self._getUniqueIndexStr()
        
        if not MINIMAL_GAMMA_DISPLAY:
            s += "[ "
            for i in range( 0, len( self.expr.terms ) ):
                s += "{ " + str( self.expr.terms[ i ] ) + "_" + str( self.indices[ i ] ) + " }"
            s += " ]"
            
        return s  
    
    """
    Private helper method that generates a string of unique indices
    referenced in this object.
    """    
    def _getUniqueIndexStr( self ):
        self.updateUniqueIndices()
        
        s = "( "
        for index in self.uniqueIndices:
            s += str( index ) + " "
            
        return s + ")"
    
    def updateUniqueIndices( self ):
        for index in self.indices:
            if not index[0] in self.uniqueIndices:
                self.uniqueIndices.append( index[0] )
                
            if not index[1] in self.uniqueIndices:
                self.uniqueIndices.append( index[1] )
                
    def getIntegratedTensorElement(self, elementIndices ):
        # Check if element will be trivially non-zero from the known matrix
        # structure of M. Mathematically equivalent to enforcing Kronocker
        # deltas \sum_i \delta_{i_\tau, i_{\tau + 1}}. Checking for structure
        # of temporal indices at this point only. Note that the index pair on
        # each \frac{\partial M}{\partial A} term is a set of nested tuples,
        # such that the particular indices are ordered by (spatial, temporal),
        # e.g. ( (x1, t1), (x2, t2) ), where of course in practice xi and ti
        # are integers that determine a tensor element.
        for indexPair in elementIndices:
            if not indexPair[1][1] - indexPair[0][1] == 1 :  # Remember, we are grabbing
                return CoefficientFloat( 0.0 )                 # the temporal indices here.
                    
        # Perform path integral for contact interactions. Note that a pair of
        # temporal indices must coincide for the term to be non-vanishing;
        # since the interaction matrix is diagonal in coordinate space, the
        # spatial components of both indices are implied as equal. Count the
        # number of spacetime (x2 t2) indices that coincide:
        spacetimeIndexCount = dict()
        for indexPair in elementIndices:
            if indexPair[1] in spacetimeIndexCount:
                spacetimeIndexCount[ indexPair[1] ] += 1
            else:
                spacetimeIndexCount[ indexPair[1] ] = 1

        # For each unique index, determine if the path integral of the
        # corresponding power of sine vanishes or is non-trivial, and
        # generate the appropriate expression.
        integralResult = Product([])         
        for spacetimeIndex in spacetimeIndexCount:
            spacetimeCount = spacetimeIndexCount[ spacetimeIndex ]
            if spacetimeCount % 2 == 1:  # If temporalCount is odd.
                return CoefficientFloat( 0.0 )  # Integral vanishes.
            else:
                integralResult.addTerm(  SINE_PATH_INTEGRALS[ spacetimeCount ] )
                
        # If we haven't returned yet, the integral is non-vanishing, and we
        # need to add the kinetic energy matrices which were implicitly
        # factored out in the previous steps. The indicies of these matrices
        # are given by the corresponding spatial components (e.g. (x1 x2)).
        for indexPair in elementIndices:
            integralResult.addTerm( MatrixT( (indexPair[0][0], indexPair[1][0]) ) )
            
        return integralResult
          
    """
    Add a term to the product contained in this object.
    @param: term Term to add to the end of the product.
    @param: index Two-element tuple that references the indices of this term.
    """
    def addTerm( self, term, index ):
        self.expr.addTerm( term )
        self.indices.append( index )
    
# ****************************************************************************
#   PRIVATE BASIC HELPER AND EXPRESSION MANIPULATION FUNCTIONS
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
        
"""
Multiplies two M matrices as a binary product.
@param: X First MatrixM term.
@param: Y Second MatrixM term.
@return: The product of X and Y.
"""   
def multiplyMMatrices( X, Y ):
    # Enforce that the dimensions of matrices X and Y must be the same. If not,
    # raise an exception.
    if not X.Nx == Y.Nx:
        raise PTSymbolicException( "Spatial volume Nx of M matrices do not match in matrix multiplication." )
    elif not X.Ntau == Y.Ntau:
        raise PTSymbolicException( "Temporal volume Ntau of M matrices do not match in matrix multiplication." )
    elif not X.spatialDimension == Y.spatialDimension:
        raise PTSymbolicException( "Spatial dimension of M matrices do not match in matrix multiplication." )
    
    # Carry on with matrix multiplication.
    Z = MatrixM( X.Nx, X.Ntau, X.spatialDimension )
    
    for i in range( 0, X.Ntau ):
        for j in range( 0, X.Ntau ):
            newMatrixElement = Sum([])
            for k in range( 0, X.Ntau ):
                newMatrixElement.addTerm( Product( [X.getElement( i, k ), Y.getElement( k, j )] ) )
            
            newMatrixElement.simplify()
            Z.setElement( i, j, unpackTrivialExpr( newMatrixElement ) )
            
    return Z

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
Converts all traces in an expression to indexed traces.
@param: expr An expression that is a Sum.
@return: A Sum() whose Trace objects are now IndexedTrace objects.
"""
def indexExpr( expr ):
    indexedExpr = Sum([])
    if isinstance( expr, Sum ):
        for term in expr.terms:
            if isinstance( term, Product ):
                nextIndexedTerm = Product([])
                nextIndexedTrace = None
                nextStartingIndex = 0
                for factor in term.terms:
                    if isinstance( factor, Trace ):
                        nextIndexedTrace = IndexedTrace( factor.expr, nextStartingIndex )
                        nextStartingIndex = nextIndexedTrace.endingIndex + 1
                        nextIndexedTerm.addTerm( nextIndexedTrace )
                    else:
                        nextIndexedTerm.addTerm( factor )
                indexedExpr.addTerm( nextIndexedTerm )
            elif isinstance( term, Trace ):
                indexedExpr.addTerm( IndexedTrace( term.expr ) )
            else:
                indexedExpr.addTerm( term )
    else:
        raise PTSymbolicException( "indexExpr() expects a Sum as the passed expression." )
    
    return indexedExpr

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

"""
Private helper function to determine whether all elements in a list are equal.
Typical use is within truncateSingleFlavorTerms().
@param: lst List to check contents of.
@return: True if all elements in lst are equal, False otherwise.
"""
def _elementsAreAllEqual( lst ):
    if len( lst ) == 0 or len( lst ) == 1:
        return True
    
    for i in range( 0, len( lst ) - 1 ):
        if not lst[ i ] == lst[ i + 1 ]:
            return False
        
    return True

"""
Terms which have matrices under a single particle flavor are shown to vanish
under the path integral, therefore we truncate them from the expansion.
@param: sum Expansion to truncate terms from.
@return: The truncated expansion.
"""        
def truncateSingleFlavorTerms( sum ):
    truncatedSum = Sum([])
    
    for term in sum.terms:
        # Obtain all flavor labels for each matrix that appears in this term.
        flavorLabels = []
        for factor in term.terms:
            if isinstance( factor, Trace ):
                for subterm in factor.expr.terms:
                    if isinstance( subterm, MatrixM ) or isinstance( subterm, MatrixB ):
                        flavorLabels.append( subterm.flavorLabel )
            
        # Check list of flavor labels; if all labels are equivalent, omit the
        # term from the sum. Keep terms that do not have matrices (the
        # non-interacting contribution).
        if not _elementsAreAllEqual( flavorLabels ) or len( flavorLabels ) == 0:
            truncatedSum.addTerm( term )
            
    return truncatedSum

"""
Public method to evaluate the product of all coefficients in an expression.
Passed expression can either be a Sum or a Product; the method is called
recursively as appropriate. If a Sum is passed, it is expected to be fully
distributed first.
@param: expr Expression whose coefficients are to be evaluated.
@return: Expression where each product contains at most one evaluated coefficient. 
"""
def evaluateCoefficients( expr ):
    if isinstance( expr, Product ):
        simplifiedProduct = Product([])
        coefficient = 1.0
        for term in expr.terms:
            if isinstance( term, CoefficientFloat ):
                coefficient *= term.value
            elif isinstance( term, CoefficientFraction ):
                coefficient *= term.eval()
            else:
                simplifiedProduct.addTerm( term )
                
        # If the coefficient is 1.0, don't bother adding it the expression.
        if not coefficient == 1.0:
            simplifiedProduct.addTerm( CoefficientFloat( coefficient ) )
           
        return simplifiedProduct
    
    elif isinstance( expr, Sum ):
        simplifiedSum = Sum([])
        for term in expr.terms:
            simplifiedSum.addTerm( evaluateCoefficients( term ) )
            
        return simplifiedSum

"""
Factors inverted M matrices (MatrixB) from IndexedTraces and produces an
IndexedProduct with all associated indices carried along. Scalar objects
do not carry any indices (they are assigned None in the index data
structure).
"""
def _factorInvertedMatrices( product ):
    factoredProduct = IndexedProduct([], [])
    gammaProduct = Gamma(Product([]), [])
    for factor in product.terms:
        if isinstance( factor, IndexedTrace ):
            # Check for an error case - raise exception if needed.
            if not isinstance( factor.expr, Product ):
                raise PTSymbolicException( "All traces passed to _factorInvertedMatrices must be fully distributed; arguments of traces must be a product." )
            
            i = 0    
            for element in factor.expr.terms:
                if isinstance( element, MatrixB ):
                    factoredProduct.addTerm( element, factor.indices[ i ] )
                elif isinstance( element, MatrixM ):
                    gammaProduct.addTerm( element, factor.indices[ i ] )
                i += 1
                
        else:
            factoredProduct.addTerm( factor, None )
    
    if not len( gammaProduct.expr.terms ) == 0:
        factoredProduct.addTerm( gammaProduct, None )       
         
    return factoredProduct

"""
Public method to factor B matrices from indexed traces and place all M matrices
in a Gamma object. Indices are carried throughout. The passed expression must
be a fully distributed sum, and all traces must already be indexed.
@param: expr Sum to be modified.
@return: The modified expression (a Sum).
"""    
def generateGammaExpression( expr ):
    if isinstance( expr, Sum ):
        newSum = Sum([])
        for term in expr.terms:
            if not isinstance( term, Product ):
                raise PTSymbolicException( "Expression passed to generateGammaExpression must be fully distributed; expr must be a Sum of Products." )
            newSum.addTerm( _factorInvertedMatrices( term ))
            
    return newSum


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
#   INTEGRATION ROUTINE HELPER FUNCTIONS
# ****************************************************************************

def calculateAllContractions( n ):
    if not n % 2:
        raise PTSymbolicException( "Parameter 'n' (index vector dimension) must be an even integer." )
    
    


                