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
SPLIT_SUMS_BY_LINE = False # If true, each term of a sum will be output on a separate late.

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
    def __init__( self, Nx, Ntau, spatialDimension, isInteracting=True, flavorLabel="" ):
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
        D = MatrixM( self.Nx, self.Ntau, self.spatialDimension, self.isInteracting, self.flavorLabel )
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
            
    def simplify( self ):
        pass
    
    def setAsNoninteracting( self ):
        self.isInteracting = False
        
    
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
        self.reduceTree()
        for term in self.terms:
            if isinstance( term, Product ):
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
        if len( self.terms ) == 0 or len( self.terms ) == 1:
            # No expansion to be done; return the original object.
            return self
        
        elif len( self.terms ) == 2:
            if isinstance( unpackTrivialExpr( self.terms[0] ), Sum ):  # Case 1: First term is a Sum, or both terms are a Sum.
                expandedSum = Sum([])
                for i in range( 0, len( self.terms[0].terms ) ):
                    expandedSum.terms.append( Product( [ deepcopy( self.terms[0].terms[i] ), deepcopy( self.terms[1] ) ] ).getExpandedExpr() )
                return expandedSum
            
            elif isinstance( unpackTrivialExpr( self.terms[1] ), Sum ):  # Case 2: Second term is a Sum.
                expandedSum = Sum([])
                for i in range( 0, len( self.terms[1].terms ) ):
                    expandedSum.terms.append( Product( [ deepcopy( self.terms[0] ), deepcopy( self.terms[1].terms[i] ) ] ).getExpandedExpr() )
                return expandedSum 
            else:  # Case 3: Neither term is a sum, so no expansion is necessary.
                return self
            
        else:  # Recursively expand the product.
            return Product( [ Product( self.terms[0:2] ).getExpandedExpr(), Product( self.terms[2:] ).getExpandedExpr() ] ).getExpandedExpr()
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

    def reduceTree( self ):
        if isinstance( self.expr, Product ) or isinstance( self.expr, Sum ):
            self.expr.reduceTree()
            
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
        self.startingIndex = 0
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
            self.endingIndex = i + 1
         
    def __str__( self ):
        s = "IndTrace[ "
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
    if not isinstance( tr, Trace ):
        raise PTSymbolicException( "An object other then a Trace was passed to distributeTrace()." )
    
    # Simple case: a single product inside a trace has been passed in.
    if isinstance( tr.expr, Product ):
        # Consider relatively more complicated cases, where the argument to the
        # trace is an expression like a ( b c + d e ); we must expand the product
        # and reduce the representation.
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
                    simplifiedProduct.terms.append( deepcopy( term ) )
                elif isinstance( term, MatrixB ) or isinstance( term, MatrixM ):
                    matrixTerms.append( deepcopy( term ) )
    
            simplifiedProduct.terms.append( Trace( unpackTrivialExpr( Product( matrixTerms ) ) ) )
            return simplifiedProduct
    
    # Typical case: a sum of terms inside a trace has been passed in to be evaluated.
    elif isinstance( tr.expr, Sum ):
        distributedSum = Sum([])
        for term in tr.expr.terms:
            distributedSum.addTerm( distributeTrace( Trace( deepcopy( term ) ) ) )
            
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
    # TODO: Need to rewrite/improve this algorithm.
    # First, before doing anything, check if the expression can
    # be distributed as a whole:
    if isinstance( expr, Product ):
        expr = expr.getExpandedExpr()
        expr.reduceTree()
        
    distributedExpr = Sum([])
    if isinstance( expr, Sum ):
        for term in expr.terms:
            if isinstance( term, Product ):
                distributedProduct = Product([])
                for factor in term.terms:
                    if isinstance( factor, Trace ):
                        distributedProduct.addTerm( distributeTrace( factor ) )
                    else:
                        distributedProduct.addTerm( factor )
                distributedProduct = distributedProduct.getExpandedExpr()
                distributedExpr.addTerm( distributedProduct )
            elif isinstance( term, Trace ):
                distributedExpr.addTerm( distributeTrace( term ) )
                distributedExpr.reduceTree()
            else:
                distributedExpr.addTerm( term )
                
    elif isinstance( expr, Product ):
        distributedExpr = Product([])
        for factor in expr.terms:
            if isinstance( factor, Trace ):
                distributedExpr.addTerm( distributeTrace( factor ) )
            else:
                distributedExpr.addTerm( factor )
        distributedExpr = distributedExpr.getExpandedExpr()
        
    else:
        return expr
    
    distributedExpr.reduceTree()
    return distributedExpr 

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
          
# ****************************************************************************
#   INTEGRATION ROUTINE HELPER FUNCTIONS
# ****************************************************************************

def calculateAllContractions( n ):
    if not n % 2:
        raise PTSymbolicException( "Parameter 'n' (index vector dimension) must be an even integer." )
    
    


                