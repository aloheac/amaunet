import PTSymbolicObjects as pt
import PTMultithreading as ptmt
import sys
import PTOutput as ptout 

Mup = pt.MatrixM( True, "" )
Mdn = pt.MatrixM( True, "" )

detMup = pt.DetM( Mup )
detMdn = pt.DetM( Mdn )

print "Calculating first derivative..."
D1up = detMup.derivative()
D1dn = detMdn.derivative()
D1up.reduceTree()
D1dn.reduceTree()
D1up = D1up.getExpandedExpr()
D1dn = D1dn.getExpandedExpr()
D1up.reduceTree()
D1dn.reduceTree()
D1up = pt.distributeAllTraces( D1up )
D1dn = pt.distributeAllTraces( D1dn )
D1up.reduceTree()
D1dn.reduceTree()
D1up = D1up.getExpandedExpr()
D1dn = D1dn.getExpandedExpr()
D1up.reduceTree()
D1dn.reduceTree()

print "Calculating second derivative..."
D2up = D1up.derivative()
D2dn = D1dn.derivative()
D2up.reduceTree()
D2dn.reduceTree()
D2up = D2up.getExpandedExpr()
D2dn = D2dn.getExpandedExpr()
D2up.reduceTree()
D2dn.reduceTree()
D2up = pt.distributeAllTraces( D2up )
D2dn = pt.distributeAllTraces( D2dn )
D2up.reduceTree()
D2dn.reduceTree()
D2up.reduceTree()
D2dn.reduceTree()
D2up = D2up.getExpandedExpr()
D2dn = D2dn.getExpandedExpr()

print "Calculating third derivative..."
D3up = D2up.derivative()
D3dn = D2dn.derivative()
D3up.reduceTree()
D3dn.reduceTree()
D3up = D3up.getExpandedExpr()
D3dn = D3dn.getExpandedExpr()
D3up.reduceTree()
D3dn.reduceTree()
D3up = pt.distributeAllTraces( D3up )
D3dn = pt.distributeAllTraces( D3dn )
D3up.reduceTree()
D3dn.reduceTree()
D3up = D3up.getExpandedExpr()
D3dn = D3dn.getExpandedExpr()
D3up.reduceTree()
D3dn.reduceTree()

print "Calculating fourth derivative..."
D4up = D3up.derivative()
D4dn = D3dn.derivative()
D4up.reduceTree()
D4dn.reduceTree()
D4up = D4up.getExpandedExpr()
D4dn = D4dn.getExpandedExpr()
D4up.reduceTree()
D4dn.reduceTree()
D4up = pt.distributeAllTraces( D4up )
D4dn = pt.distributeAllTraces( D4dn )
D4up.reduceTree()
D4dn.reduceTree()
D4up = D4up.getExpandedExpr()
D4dn = D4dn.getExpandedExpr()
D4up.reduceTree()
D4dn.reduceTree()

# print "Calculating fifth derivative..."
# D5up = D4up.derivative()
# D5dn = D4dn.derivative()
# D5up.reduceTree()
# D5dn.reduceTree()
# D5up = D5up.getExpandedExpr()
# D5dn = D5dn.getExpandedExpr()
# D5up.reduceTree()
# D5dn.reduceTree()
# D5up = pt.distributeAllTraces( D5up )
# D5dn = pt.distributeAllTraces( D5dn )
# D5up.reduceTree()
# D5dn.reduceTree()
     
# print "Calculating sixth derivative..."
# D6up = D5up.derivative()
# D6dn = D5dn.derivative()
# D6up.reduceTree()
# D6dn.reduceTree()
# D6up = D6up.getExpandedExpr()
# D6dn = D6dn.getExpandedExpr()
# D6up.reduceTree()
# D6dn.reduceTree()
# D6up = pt.distributeAllTraces( D6up )
# D6dn = pt.distributeAllTraces( D6dn )
# D6up.reduceTree()
# D6dn.reduceTree()

# print "Calculating seventh derivative..."
# print ">> Taking derivatives..."
# D7up = D6up.derivative()
# D7dn = D6dn.derivative()
# D7up.reduceTree()
# D7dn.reduceTree()
# print ">> Expanding..."
# D7up = D7up.getExpandedExpr()
# D7dn = D7dn.getExpandedExpr()
# D7up.reduceTree()
# D7dn.reduceTree()
# print ">> Distributing traces..."
# D7up = pt.distributeAllTraces( D7up )
# D7dn = pt.distributeAllTraces( D7dn )
# D7up.reduceTree()
# D7dn.reduceTree()

# print "Calculating eighth derivative..."
# print ">> Taking derivatives..."
# D8up = D7up.derivative()
# D8dn = D7dn.derivative()
# D8up.reduceTree()
# D8dn.reduceTree()
# print ">> Expanding..."
# D8up = D8up.getExpandedExpr()
# D8dn = D8dn.getExpandedExpr()
# D8up.reduceTree()
# D8dn.reduceTree()
# print ">> Distributing traces..."
# D8up = pt.distributeAllTraces( D8up )
# D8dn = pt.distributeAllTraces( D8dn )
# D8up.reduceTree()
# D8dn.reduceTree()

print "Rewriting derivatives in KS formalism..."
D1up = pt.rewriteExprInKSFormalism( pt.Sum( [D1up] ) )
D1dn = pt.rewriteExprInKSFormalism( pt.Sum( [D1dn] ) )
D1up.reduceTree()
D1dn.reduceTree()
D2up = pt.rewriteExprInKSFormalism( D2up )
D2dn = pt.rewriteExprInKSFormalism( D2dn )
D2up.reduceTree()
D2dn.reduceTree()
D3up = pt.rewriteExprInKSFormalism( D3up )
D3dn = pt.rewriteExprInKSFormalism( D3dn )
D3up.reduceTree()
D3dn.reduceTree()
D4up = pt.rewriteExprInKSFormalism( D4up )
D4dn = pt.rewriteExprInKSFormalism( D4dn )
D4up.reduceTree()
D4dn.reduceTree()
# D5up = pt.rewriteExprInKSFormalism( D5up )
# D5dn = pt.rewriteExprInKSFormalism( D5dn )
# D5up.reduceTree()
# D5dn.reduceTree()
# D6up = pt.rewriteExprInKSFormalism( D6up )
# D6dn = pt.rewriteExprInKSFormalism( D6dn )
# D6up.reduceTree()
# D6dn.reduceTree()
# D7up = pt.rewriteExprInKSFormalism( D7up )
# D7dn = pt.rewriteExprInKSFormalism( D7dn )
# D7up.reduceTree()
# D7dn.reduceTree()
# D8up = pt.rewriteExprInKSFormalism( D8up )
# D8dn = pt.rewriteExprInKSFormalism( D8dn )
# D8up.reduceTree()
# D8dn.reduceTree()

print "Setting all terms as non-interacting..."
detMup.setAsNoninteracting()
detMdn.setAsNoninteracting()
D1up.setAsNoninteracting()
D1dn.setAsNoninteracting()
D2up.setAsNoninteracting()
D2dn.setAsNoninteracting()
D3up.setAsNoninteracting()
D3dn.setAsNoninteracting()
D4up.setAsNoninteracting()
D4dn.setAsNoninteracting()
# D5up.setAsNoninteracting()
# D5dn.setAsNoninteracting()
# D6up.setAsNoninteracting()
# D6dn.setAsNoninteracting()
# D7up.setAsNoninteracting()
# D7dn.setAsNoninteracting()
# D8up.setAsNoninteracting()
# D8dn.setAsNoninteracting()

print "Calculating spin-up partition function..."
Zup = pt.Sum([detMup, pt.Product([ pt.TermA(), D1up ]), pt.Product([ pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,2), D2up ]) , \
               pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,6), D3up]), \
               pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,24), D4up])])#, \
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,120), D5up]), \
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,720), D6up]), \
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,5040), D7up]),
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,40320), D8up]) ])

Zup = Zup.getExpandedExpr()
Zup.reduceTree()

print "Calculating spin-down partition function..."
Zdn = pt.Sum([detMdn, pt.Product([ pt.TermA(), D1dn ]), pt.Product([ pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,2), D2dn ]), \
               pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,6), D3dn]), \
               pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,24), D4dn])])#, \
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,120), D5dn]), \
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,720), D6dn]), \
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,5040), D7dn]),
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,40320), D8dn]) ])

Zdn.reduceTree()
Zdn.getExpandedExpr()
Zdn.reduceTree()

Z = pt.Product([Zup,Zdn])
Z.reduceTree()

print "Expanding..."
#Z = Z.getExpandedExpr()
Z = ptmt.execParallelDualExpansion( Zup, Zdn )
Z.reduceTree()
Z.reduceTree()

del D1up
del D1dn
del D2up
del D2dn
del D3up
del D3dn
del D4up
del D4dn
# del D5up
# del D5dn
# del D6up
# del D6dn
# del D7up
# del D7dn
# del D8up
# del D8dn
del Zup
del Zdn

print "Truncating odd and high-order terms..."
Z = pt.truncateOddOrders( Z )
Z = pt.truncateAOrder(Z, 4)

print "Indexing matrices..."
Z = pt.indexExpr( Z )

print "Computing path integral..."
Z = pt.pathIntegrateExpression( Z )

print "Expanding..."
Z = ptmt.execParallelSumExpansion( Z )
#Z = Z.getExpandedExpr()
Z.reduceTree()
Z.reduceTree()
print "Computing analytic Fourier transform..."
Z = pt.fourierTransformExpr( Z )

print "Combining like terms..."
Z = pt.combineLikeTerms( Z )
Z = ptmt.execParallelCombineLikeTerms( Z, 50)
Z.simplify()

print "Writing expressions to file..."
f = open( "8thAOrderExpansionExpr.txt", 'w' )
f.write( str( Z ) )
f.close()

f = open( "8thAOrderExpansionInterpreter.txt", 'w' )
f.write( ptout.outputExprToInterpreterFormat( Z ) )
f.close()

print Z
#KK = pt.Sum([Z.terms[-1]])
#KK = pt.fourierTransformExpr( KK )
#print KK
print "Done!"