"""
Amaunet - High-order Lattice Perturbation Theory Expansion 
          for Interacting Quantum Matter
Weak-coupling Expansion for Fermionic Systems with Contact Interactions

Analytic Expression Output Interface

Andrew C. Loheac, Joaquin E. Drut
Department of Physics and Astronomy
University of North Carolina at Chapel Hill
"""

import PTSymbolicObjects as pt

def outputExprToInterpreterFormat( expr ):
    if not isinstance( expr, pt.Sum ):
        raise pt.PTSymbolicException( "Expression passed to outputExprToInterpreterFormat() must be an instance of a fully-distributed Sum." )
    
    s = ""
    for term in expr.terms:
        if not isinstance( term, pt.Product ):
            raise pt.PTSymbolicException( "Each term in the Sum passed to outputExprToInterpreterFormat() must be a single Product." )
        
        line = ""
        AOrder = 0
        for factor in term.terms:
            if isinstance( factor, pt.TermA ):
                AOrder += 1
            elif isinstance( factor, pt.MatrixK ):
                line += "D," + str( factor.flavorLabel ) + "," + str( factor.indices[0] ) + "," + str( factor.indices[1] ) + "/"
            elif isinstance( factor, pt.FourierSum ):
                line += "F,"
                for index in factor.indices:
                    line += str( index[0] ) + "," + str( index[1] ) + ","
                line = line[:-1] + "/"  # Remove the trailing comma on the last index.
            elif isinstance( factor, pt.CoefficientFloat ):
                line += "C," + str( factor.eval() ) + "/"
            elif isinstance( factor, pt.Delta ) and factor.isBar == True:
                line += "B," + str( factor.indices[0] ) + "," + str( factor.indices[1] ) + "/"
            
        s += "A," + str( AOrder ) + "/" + line + ';'
    
    return s
            
def outputExprToLatexFormat( expr ):
    if not isinstance( expr, pt.Sum ):
        raise pt.PTSymbolicException( "Expression passed to outputExprToLatexFormat() must be an instance of a fully-distributed Sum." )
    
    s = ""
    for term in expr.terms:
        if not isinstance( term, pt.Product ):
            raise pt.PTSymbolicException( "Each term in the Sum passed to outputExprToLatexFormat() must be a single Product." )
        
        line = ""
        AOrder = 0
        for factor in term.terms:
            if isinstance( factor, pt.TermA ):
                AOrder += 1
            elif isinstance( factor, pt.MatrixK ):
                line += " D" + str( factor.flavorLabel ) + "," + str( factor.indices[0] ) + "," + str( factor.indices[1] ) + "} "
            elif isinstance( factor, pt.FourierSum ):
                line += "{F,"
                for index in factor.indices:
                    line += str( index[0] ) + "," + str( index[1] ) + ","
                line = line[:-1] + "} "  # Remove the trailing comma on the last index.
            elif isinstance( factor, pt.CoefficientFloat ):
                line += "{C," + str( factor.eval() ) + "} "
            
        s += "{A," + str( AOrder ) + "} " + line + '\n'
    
    return s