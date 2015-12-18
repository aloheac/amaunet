"""
Amaunet - High-order Perturbation Theory Expansion for Interacting Quantum Matter
Fermionic Systems with Contact Interactions

Andrew C. Loheac, Joaquin E. Drut
Department of Physics and Astronomy
University of North Carolina at Chapel Hill
"""

from copy import deepcopy

# Static variables that define the behavior of symbolic calculations, system parameters,
# and output.
global NX, NTAU, MINIMAL_U_DISPLAY, MINIMAL_M_DISPLAY

NX = 2                      # Spatial volume of the system.
NTAU = 2                    # Temporal volume of the system.

MINIMAL_U_DISPLAY = False    # If true, the matrix U will be output using a shorthand notation.
MINIMAL_M_DISPLAY = True   # If true, the matrix M will be output using a shorthand notation.

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
    def __init__( self, Nx, Ntau, spatialDimension, isInteracting=True ):
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
        if not MINIMAL_M_DISPLAY:
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
            if self.isInteracting:
                if self.derivativeOrder == 0:
                    return "M"
                elif self.derivativeOrder == 1:
                    return "dM / dA"
                else:
                    return "d^" + str( self.derivativeOrder ) + " M / dA^" + str( self.derivativeOrder )
                    #return "0.0"
            else:
                if self.derivativeOrder == 0:
                    return "M0"
                else:
                    return "0.0"
    
    """
    Gets the matrix elements at temporal indices ( i, j ). Returns zero if the
    indices are not found in the internal dictionary.
    @param: i First temporal index (row) of the matrix element.
    @param: j Second temporal index (column) of the matrix element.
    @return: The expression of the matrix element at indices ( i, j ).
    """            
    def getElement( self, i, j ):
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
        self.elements[ ( i, j ) ] = entry
    
    """
    Calculates the derivative of all matrix elements, and returns a new
    MatrixM object which is the derivative of this instance.
    @return: The derivative of this matrix.
    """
    def derivative( self ):
        D = MatrixM( self.Nx, self.Ntau, self.spatialDimension )
        for key in self.elements:
            D.setElement( key[0], key[1], self.elements[ key ].derivative() )
        D.simplify()
        D.derivativeOrder = self.derivativeOrder + 1
        
        return D
    
    """
    Mathematically simplifies the expressions for all matrix elements.
    """
    def simplify( self ):
        for key in self.elements:
            self.elements[ key ].simplify()
            
"""
Representation of the matrix B, or the inverse of the matrix M, in the defined
perturbation theory language.
"""
class MatrixB:
    def __init__( self, Nx, Ntau, spatialDimension, isInteracting=True ):
        self.elements = dict()
        self.Nx = Nx
        self.Ntau = Ntau
        self.spatialDimension = spatialDimension
        self.isInteracting = isInteracting
        self.derivativeOrder = 0
        
    def __str__( self ):
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
        
    def derivative( self ):
        dM = MatrixM( self.Nx, self.Ntau, self.spatialDimension, self.isInteracting )
        dM.derivativeOrder = self.derivativeOrder + 1
        
        return Product( [ CoefficientFloat( -1.0 ), MatrixB( self.Nx, self.Ntau, self.spatialDimension ), dM, MatrixB( self.Nx, self.Ntau, self.spatialDimension ) ] )
          
    def simplify( self ):
        pass
    
"""
Representation of the matrix U in the defined perturbation theory language.
"""
class MatrixU:
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
    def __init__(self, M, flavorLabel="", isInteracting=True, isInverted=False ):
        self.M = M
        self.flavorLabel = ""
        self.isInteracting = isInteracting
        self.isInverted = isInverted

    def __str__(self):
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
            
    def derivative( self ):
        if self.isInteracting:
            if self.isInverted:
                # Leave the below result in case it is useful for future implementations.
                # For now, this path should not be permitted; raise an exception.
                # return Product( [ CoefficientFloat( -1.0 ), Trace( Product( [ DetM( self.flavorLabel, True, True ), self.M.derivative() ] ) ), DetM( self.flavorLabel, True, True ) ] )
                raise PTSymbolicException( "The derivative of an inverted, interacting determinant should not be taken in this context." )
            else:
                return Product( [ DetM( self.M, self.flavorLabel, True, False ), Trace( Product( [ MatrixB( self.M.Nx, self.M.Ntau, self.M.spatialDimension ), self.M.derivative() ] ) ) ] )
        else:
            if self.isInverted:
                return "1 / Det[ M0 ]"
            else:
                return "Det[ M0 ]"
            
    def simplify( self ):
        pass
        
    
"""
Representation of a floating-point coefficient.
"""        
class CoefficientFloat:
    """
    Constructor for CoefficientFloat.
    @param: value Floating-point value of the coefficent.
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
        return "{ " + str( self.num ) + " / " + str( self.den ) + " }"
    
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
        self.terms = []
        
    def __str__( self ):
        s = ""
        for i in range( 0, len( self.terms ) ):
            s += str( self.terms[ i ] )
            if not i == len( self.terms ) - 1:
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
            
    def reduceTree( self ):    
        reducedTerms = []
        for s in self.terms:
            if isinstance( s, Sum ):
                s.reduceTree()
                for i in range( 0, len( s.terms ) ):
                    reducedTerms.append( s.terms[ i ] )
            elif isinstance( s, Product ):
                s.reduceTree()
                reducedTerms.append( s )
            else:
                reducedTerms.append( s )
                
        self.terms = reducedTerms
            
    def addTerm( self, term ):
        self.terms.append( term )
        
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
        
    def reduceTree( self ):           
        reducedTerms = []
        for term in self.terms:
            if isinstance( term, Product ):
                term.reduceTree()
                for i in range( 0, len( term.terms ) ):
                    reducedTerms.append( term.terms[ i ] )
            elif isinstance( term, Sum ):
                term.reduceTree()
                reducedTerms.append( term )    
            else:
                reducedTerms.append( term )    
                
        self.terms = reducedTerms    
        
    """
    Recursively expand a product of sums to a explicit sum of products; i.e. the
    result expression is in the form AB... + CD... + ..., where applicable.
    @returns: The expanded product (a Sum) of this Product.
    """
    def getExpandedExpr( self ):      
        if len( self.terms ) == 0 or len( self.terms ) == 1:
            # No expansion to be done; return the original object.
            return self
        
        elif len( self.terms ) == 2:
            if isinstance( unpackTrivialExpr( self.terms[0] ), Sum ):  # Case 1: First term is a Sum, or both terms are a Sum.
                expandedSum = Sum([])
                for i in range( 0, len( self.terms[0].terms ) ):
                    expandedSum.terms.append( Product( [ self.terms[0].terms[i], self.terms[1] ] ).getExpandedExpr() )
                return expandedSum
            
            elif isinstance( unpackTrivialExpr( self.terms[1] ), Sum ):  # Case 2: Second term is a Sum.
                expandedSum = Sum([])
                for i in range( 0, len( self.terms[1].terms ) ):
                    expandedSum.terms.append( Product( [ self.terms[0], self.terms[1].terms[i] ] ).getExpandedExpr() )
                return expandedSum 
            else:  # Case 3: Neither term is a sum, so no expansion is necessary.
                return self
            
        else:  # Recursively expand the product.
            return Product( [ Product( self.terms[0:2] ).getExpandedExpr(), Product( self.terms[2:] ).getExpandedExpr() ] ).getExpandedExpr()
 
    def addTerm( self, term ):
        self.terms.append( term )
        
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
    
    def derivative( self ):
        return Trace( self.expr.derivative() )
    
    """
    Simplifies the expression within the trace.
    @return: The simplified expression.
    """
    def simplify( self ):
        self.expr.simplify()
               
# ****************************************************************************
#   PRIVATE HELPER FUNCTIONS
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
            newMatrixElement = Sum()
            for k in range( 0, X.Ntau ):
                newMatrixElement.addTerm( Product( [X.getElement( i, k ), Y.getElement( k, j )] ) )
            
            newMatrixElement.simplify()
            Z.setElement( i, j, unpackTrivialExpr( newMatrixElement ) )
            
    return Z