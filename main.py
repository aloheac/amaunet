import PTSymbolicObjects as pt
import PTMultithreading as ptmt
import cProfile
import sys

Mup = pt.MatrixM( 10, 4, 1, True, "up" )
Mdn = pt.MatrixM( 10, 4, 1, True, "dn" )

detMup = pt.DetM( Mup )
detMdn = pt.DetM( Mdn )

print "Calculating first derivative..."
D1up = detMup.derivative().getExpandedExpr()
D1dn = detMdn.derivative().getExpandedExpr()
D1up.setAsNoninteracting()
D1dn.setAsNoninteracting()
D1up.reduceTree()
D1dn.reduceTree()

print "Calculating second derivative..."
D2up = detMup.derivative().derivative().getExpandedExpr()
D2dn = detMdn.derivative().derivative().getExpandedExpr()
D2up.setAsNoninteracting()
D2dn.setAsNoninteracting()
D2up.reduceTree()
D2dn.reduceTree()

print "Calculating third derivative..."
D3up = detMup.derivative().derivative().derivative().getExpandedExpr()
D3dn = detMdn.derivative().derivative().derivative().getExpandedExpr()
D3up.setAsNoninteracting()
D3dn.setAsNoninteracting()
D3up.reduceTree()
D3dn.reduceTree()

print "Calculating fourth derivative..."
D4up = detMup.derivative().derivative().derivative().derivative().getExpandedExpr()
D4dn = detMdn.derivative().derivative().derivative().derivative().getExpandedExpr()
D4up.setAsNoninteracting()
D4dn.setAsNoninteracting()
D4up.reduceTree()
D4dn.reduceTree()

# print "Calculating fifth derivative..."
# D5up = detMup.derivative().derivative().derivative().derivative().derivative().getExpandedExpr()
# D5dn = detMdn.derivative().derivative().derivative().derivative().derivative().getExpandedExpr()
# D5up.setAsNoninteracting()
# D5dn.setAsNoninteracting()
# D5up.reduceTree()
# D5dn.reduceTree()
# 
# print "Calculating sixth derivative..."
# D6up = detMup.derivative().derivative().derivative().derivative().derivative().derivative().getExpandedExpr()
# D6dn = detMdn.derivative().derivative().derivative().derivative().derivative().derivative().getExpandedExpr()
# D6up.reduceTree()
# D6dn.reduceTree()
# 
# print "Calculating seventh derivative..."
# D7up = D6up.derivative().getExpandedExpr()
# D7dn = D6dn.derivative().getExpandedExpr()
# D7up.reduceTree()
# D7dn.reduceTree()
# 
# print "Calculating eighth derivative..."
# D8up = D7up.derivative().getExpandedExpr()
# D8dn = D7dn.derivative().getExpandedExpr()
# D8up.reduceTree()
# D8dn.reduceTree()
# 
# D6up.setAsNoninteracting()
# D6dn.setAsNoninteracting()
# D7up.setAsNoninteracting()
# D7dn.setAsNoninteracting()
# D8up.setAsNoninteracting()
# D8dn.setAsNoninteracting()

detMup.setAsNoninteracting()
detMdn.setAsNoninteracting()

print "Calculating spin-up partition function..."
#Zup = pt.Sum([detMup, pt.Product([ pt.TermA(), D1up ]), pt.Product([ pt.TermA(), pt.TermA(), D2up ]), pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,6), D3up]), pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,24), D4up]) ])
#Zdn = pt.Sum([detMdn, pt.Product([ pt.TermA(), D1dn ]), pt.Product([ pt.TermA(), pt.TermA(), D2dn ]), pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,6), D3dn]), pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,24), D4dn]) ])

Zup = pt.Sum([detMup, pt.Product([ pt.TermA(), D1up ]), pt.Product([ pt.TermA(), pt.TermA(), D2up ]) , \
               pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,6), D3up]), \
               pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,24), D4up]) ])#, \
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,120), D5up]), \
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,720), D6up]), \
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,5040), D7up]),
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,40320), D8up]) ])
Zup.reduceTree()
Zup = Zup.getExpandedExpr()
Zup.reduceTree()

print "Calculating spin-down partition function..."
Zdn = pt.Sum([detMdn, pt.Product([ pt.TermA(), D1dn ]), pt.Product([ pt.TermA(), pt.TermA(), D2dn ]), \
               pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,6), D3dn]), \
               pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,24), D4dn]) ])#, \
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,120), D5dn]), \
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,720), D6dn]), \
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,5040), D7dn]),
               #pt.Product([pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.TermA(), pt.CoefficientFraction(1,40320), D8dn]) ])
Zdn.reduceTree()
Zdn.getExpandedExpr()
Zdn.reduceTree()

Z = pt.Product([Zup,Zdn])
Z.reduceTree()

print "Expanding stuff..."
#Z = Z.getExpandedExpr()
Z = ptmt.execParallelDualExpansion( Zup, Zdn )
Z.reduceTree()
print "Truncating odd and high-order terms..."
Z = pt.truncateOddOrders( Z )
Z = pt.truncateAOrder(Z, 4)

print "Distributing traces..."
Z = ptmt.execParallelTraceDistribution( Z )

print "Truncating..."
#Z = pt.truncateSingleFlavorTerms( Z )
Z = pt.indexExpr( Z )

print "Generating Gamma expressions..."
Z = pt.generateGammaExpression( Z )

print "Computing analytic Fourier transform..."
Z = pt.fourierTransformExpr( Z )
Z = Z.getExpandedExpr()
Z = Z.getExpandedExpr()
Z.reduceTree()
print Z
