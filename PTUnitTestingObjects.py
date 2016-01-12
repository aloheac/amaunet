"""
Amaunet - High-order Perturbation Theory Expansion for Interacting Quantum Matter
Unit Testing Object Definitions

Andrew C. Loheac, Joaquin E. Drut
Department of Physics and Astronomy
University of North Carolina at Chapel Hill

NOTE: The classes defined in this module are intended only for the unit testing
      of production code. It is not to be part of a release.
"""

import PTSymbolicObjects as pt
from PTSymbolicObjects import MINIMAL_M_DISPLAY

class GenericTestTerm:
    def __init__( self, id, derivativeOrder ):
        self.id = id
        self.derivativeOrder = derivativeOrder
        
    def __str__( self ):
        return "GT" + str( self.id ) + "_d" + str( self.derivativeOrder )
    
    def simplify( self ):
        pass
    
    def derivative( self ):
        return GenericTestTerm( self.id, self.derivativeOrder + 1 )
        
class UnitTest:
    
    testsPassed = 0
    testsFailed = 0
    
    def __init__( self, name, func, correctResult ):
        self.name = name
        self.func = func
        self.correctResult = correctResult
        
    def runMe( self ):
        print "Running unit test '" + str( self.name ) + "'..."
        result = self.func()
        if result == self.correctResult:
            print ">> Unit test PASSED."
            UnitTest.testsPassed += 1
        else:
            print ">> Unit test FAILED. Returned result:"
            print str( result )
            print "   Correct result:"
            print str( self.correctResult )
            UnitTest.testsFailed += 1
            
def executeAllUnitTests():
    # Tests for Sum.
    
    def fA01():
        A = pt.Sum([])
        A.addTerm( GenericTestTerm(1,0) )
        return str( A )
    
    TA01 = UnitTest( "TA01: Sum Basic Construction", fA01, "GT1_d0" )
    TA01.runMe()
    
    def fA02():
        A = pt.Sum([])
        A.addTerm( GenericTestTerm(1,0) )
        A.addTerm( GenericTestTerm(2,0) )
        A.addTerm( GenericTestTerm(3,0) )
        return str( A )
    
    TA02 = UnitTest( "TA02: Sum Construction with 3 Terms", fA02, "GT1_d0 + GT2_d0 + GT3_d0" )
    TA02.runMe()
    
    def fA03():
        A = pt.Sum([])
        A.addTerm( pt.CoefficientFloat(0.0) )
        A.addTerm( GenericTestTerm(1,0) )
        A.addTerm( GenericTestTerm(2,0) )
        A.addTerm( pt.CoefficientFloat(0.0) )
        A.addTerm( GenericTestTerm(3,0) )
        A.addTerm( pt.CoefficientFloat(0.0) )
        A.simplify()
        return str( A )
    
    TA03 = UnitTest( "TA03: Sum simplify()", fA03, "GT1_d0 + GT2_d0 + GT3_d0" )
    TA03.runMe()
    
    def fA04():
        A = pt.Sum([])
        A.addTerm( GenericTestTerm(1,0) )
        A.addTerm( pt.Sum( [ GenericTestTerm(2,0), GenericTestTerm(3,0) ] ) )
        A.addTerm( pt.Sum( [ GenericTestTerm(4,0), pt.Sum( [ GenericTestTerm(5,0), GenericTestTerm(6,0) ] ) ] ) )
        A.reduceTree()
        return str( A )
        
    TA04 = UnitTest( "TA04: Sum reduceTree(), Basic I", fA04, "GT1_d0 + GT2_d0 + GT3_d0 + GT4_d0 + GT5_d0 + GT6_d0" )
    TA04.runMe()
    
    def fA05():
        A = pt.Sum([])
        A.addTerm( GenericTestTerm(1,0) )
        A.addTerm( pt.Sum( [ GenericTestTerm(2,0), GenericTestTerm(3,0) ] ) )
        A.addTerm( pt.Sum( [ pt.Sum( [ GenericTestTerm(4,0), GenericTestTerm(5,0) ] ), pt.Sum( [ pt.Sum( [ GenericTestTerm(6,0), GenericTestTerm(7,0) ] ) ] ), GenericTestTerm(8,0) ] ) )
        A.reduceTree()
        return str( A )
        
    TA05 = UnitTest( "TA05: Sum reduceTree(), Basic II", fA05, "GT1_d0 + GT2_d0 + GT3_d0 + GT4_d0 + GT5_d0 + GT6_d0 + GT7_d0 + GT8_d0" )
    TA05.runMe()
    
    def fA06():
        A = pt.Sum([])
        A.addTerm( GenericTestTerm(1,0) )
        A.addTerm( GenericTestTerm(2,0) )
        A.addTerm( GenericTestTerm(3,0) )
        B = A.derivative()
        return str( B )
        
    TA06 = UnitTest( "TA06: Sum derivative(), Basic I", fA06, "GT1_d1 + GT2_d1 + GT3_d1" )
    TA06.runMe()
    
    def fA07():
        A = pt.Sum([])
        A.addTerm( GenericTestTerm(1,0) )
        A.addTerm( GenericTestTerm(2,0) )
        A.addTerm( GenericTestTerm(3,0) )
        B = A.derivative().derivative()
        return str( B )
        
    TA07 = UnitTest( "TA07: Sum derivative(), Basic II", fA07, "GT1_d2 + GT2_d2 + GT3_d2" )
    TA07.runMe()
    
    def fA08():
        T1 = pt.Product( [ GenericTestTerm(1,0), GenericTestTerm(2,0) ] )
        T2 = pt.Product( [ GenericTestTerm(3,0), GenericTestTerm(4,0) ] )
        A = pt.Sum( [ T1, T2 ] )
        B = A.derivative()
        return str( B )
    
    TA08 = UnitTest( "TA08: Sum derivative(), Products I", fA08, " {GT1_d1}  {GT2_d0} +  {GT1_d0}  {GT2_d1} +  {GT3_d1}  {GT4_d0} +  {GT3_d0}  {GT4_d1}" )
    TA08.runMe()
    
    def fA09():
        T1 = pt.Product( [ GenericTestTerm(1,0), GenericTestTerm(2,0) ] )
        T2 = pt.Product( [ GenericTestTerm(3,0), GenericTestTerm(4,0) ] )
        A = pt.Sum( [ T1, T2 ] )
        B = A.derivative().derivative()
        return str( B )
    
    TA09 = UnitTest( "TA09: Sum derivative(), Products II", fA09, " {GT1_d2}  {GT2_d0} +  {GT1_d1}  {GT2_d1} +  {GT1_d1}  {GT2_d1} +  {GT1_d0}  {GT2_d2} +  {GT3_d2}  {GT4_d0} +  {GT3_d1}  {GT4_d1} +  {GT3_d1}  {GT4_d1} +  {GT3_d0}  {GT4_d2}" )
    TA09.runMe()
    
    def fA10():
        T1 = pt.Product( [ GenericTestTerm(1,0), GenericTestTerm(2,0) ] )
        T2 = pt.Product( [ GenericTestTerm(3,0), GenericTestTerm(4,0) ] )
        A = pt.Sum( [ T1, T2 ] )
        B = A.derivative().derivative().derivative()
        return str( B )
    
    TA10 = UnitTest( "TA10: Sum derivative(), Products III", fA10, " {GT1_d3}  {GT2_d0} +  {GT1_d2}  {GT2_d1} +  {GT1_d2}  {GT2_d1} +  {GT1_d1}  {GT2_d2} +  {GT1_d2}  {GT2_d1} +  {GT1_d1}  {GT2_d2} +  {GT1_d1}  {GT2_d2} +  {GT1_d0}  {GT2_d3} +  {GT3_d3}  {GT4_d0} +  {GT3_d2}  {GT4_d1} +  {GT3_d2}  {GT4_d1} +  {GT3_d1}  {GT4_d2} +  {GT3_d2}  {GT4_d1} +  {GT3_d1}  {GT4_d2} +  {GT3_d1}  {GT4_d2} +  {GT3_d0}  {GT4_d3}" )
    TA10.runMe()
    
    def fA11():
        T1 = pt.Product( [ GenericTestTerm(1,0), GenericTestTerm(2,0) ] )
        T2 = pt.Product( [ GenericTestTerm(3,0), GenericTestTerm(4,0) ] )
        A = pt.Sum( [ T1, T2 ] )
        B = A.derivative().derivative().derivative()
        B.reduceTree()
        return str( B )
    
    TA11 = UnitTest( "TA11: Sum derivative(), reduceTree()", fA11, " {GT1_d3}  {GT2_d0} +  {GT1_d2}  {GT2_d1} +  {GT1_d2}  {GT2_d1} +  {GT1_d1}  {GT2_d2} +  {GT1_d2}  {GT2_d1} +  {GT1_d1}  {GT2_d2} +  {GT1_d1}  {GT2_d2} +  {GT1_d0}  {GT2_d3} +  {GT3_d3}  {GT4_d0} +  {GT3_d2}  {GT4_d1} +  {GT3_d2}  {GT4_d1} +  {GT3_d1}  {GT4_d2} +  {GT3_d2}  {GT4_d1} +  {GT3_d1}  {GT4_d2} +  {GT3_d1}  {GT4_d2} +  {GT3_d0}  {GT4_d3}" )
    TA11.runMe()
    
    def fA12():
        T1 = pt.Product( [ GenericTestTerm(1,0), GenericTestTerm(2,0) ] )
        T2 = pt.Product( [ GenericTestTerm(3,0), GenericTestTerm(4,0) ] )
        A = pt.Sum( [ T1, pt.CoefficientFloat( 0.0 ), T2 ] )
        B = A.derivative().derivative().derivative()
        B.reduceTree()
        B.simplify()
        return str( B )
    
    TA12 = UnitTest( "TA12: Sum derivative(), reduceTree(), simplify()", fA12, " {GT1_d3}  {GT2_d0} +  {GT1_d2}  {GT2_d1} +  {GT1_d2}  {GT2_d1} +  {GT1_d1}  {GT2_d2} +  {GT1_d2}  {GT2_d1} +  {GT1_d1}  {GT2_d2} +  {GT1_d1}  {GT2_d2} +  {GT1_d0}  {GT2_d3} +  {GT3_d3}  {GT4_d0} +  {GT3_d2}  {GT4_d1} +  {GT3_d2}  {GT4_d1} +  {GT3_d1}  {GT4_d2} +  {GT3_d2}  {GT4_d1} +  {GT3_d1}  {GT4_d2} +  {GT3_d1}  {GT4_d2} +  {GT3_d0}  {GT4_d3}" )
    TA12.runMe()
    
    def fA13():
        A = pt.Product( [ GenericTestTerm(1,0), GenericTestTerm(2,0), GenericTestTerm(3,0) ] )
        return str( A )
    
    TA13 = UnitTest( "TA13: Product, Basic Construction", fA13, " {GT1_d0}  {GT2_d0}  {GT3_d0}" )
    TA13.runMe()
    
    def fB01():
        A = pt.Product( [ pt.TermA(), pt.Trace( pt.Sum( [ GenericTestTerm(1,0), GenericTestTerm(2,0) ] ) ) ] )
        B = pt.distributeAllTraces( A )
        return str( B )
    
    TB01 = UnitTest( "TB01: distributeAllTraces(), Distribution I", fB01, " {A}  {Trace[ GT1_d0 ]} +  {A}  {Trace[ GT2_d0 ]}" )
    TB01.runMe()
    
    def fB02():
        A = pt.Trace( pt.Product( [ pt.CoefficientFloat( -1.0 ), pt.MatrixM( 1, 1, 1 ) ] ) )
        B = pt.distributeTrace( A )
        return str( B )
    
    TB02 = UnitTest( "TB02: distributeTrace(), Simplification I", fB02, " {-1.0}  {Trace[ M ]}" )
    TB02.runMe()
    
    def fB03():
        A = pt.Trace( pt.Product( [ pt.CoefficientFloat( -1.0 ), pt.Sum( [ pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) )
        B = pt.distributeTrace( A )
        return str( B )
    
    TB03 = UnitTest( "TB03: distributeTrace(), Simplification II", fB03, " {-1.0}  {Trace[ M ]} +  {-1.0}  {Trace[ M ]}" )
    TB03.runMe()
    
    def fB04():
        A = pt.Trace( pt.Product( [ pt.CoefficientFloat( -1.0 ), pt.Sum( [ pt.Product( [ pt.CoefficientFloat( 3.0 ), pt.MatrixM( 1, 1, 1 ) ] ), pt.Product( [ pt.CoefficientFloat( 5.0 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) ] ) )
        B = pt.distributeTrace( A )
        return str( B )
    
    TB04 = UnitTest( "TB04: distributeTrace(), Simplification III", fB04, " {-1.0}  {3.0}  {Trace[ M ]} +  {-1.0}  {5.0}  {Trace[ M ]}" )
    TB04.runMe()
    
    def fB05():
        A = pt.Product( [ pt.Trace( pt.Product( [ pt.CoefficientFloat( -1.0 ), pt.Product( [ pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) ), pt.Trace( pt.Product( [ pt.CoefficientFloat( -3.0 ), pt.Product( [ pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) ) ] )
        B = pt.distributeAllTraces( A )
        return str( B )
    
    TB05 = UnitTest( "TB05: distributeAllTraces(), Distribution II", fB05, " {-1.0}  {Trace[  {M}  {M} ]}  {-3.0}  {Trace[  {M}  {M} ]}" )
    TB05.runMe()
    
    def fB06():
        A = pt.Sum( [ pt.Trace( pt.Product( [ pt.CoefficientFloat( -1.0 ), pt.Product( [ pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) ), pt.Trace( pt.Product( [ pt.CoefficientFloat( -3.0 ), pt.Product( [ pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) ) ] )
        B = pt.distributeAllTraces( A )
        return str( B )
    
    TB06 = UnitTest( "TB06: distributeAllTraces(), Distribution III", fB06, " {-1.0}  {Trace[  {M}  {M} ]} +  {-3.0}  {Trace[  {M}  {M} ]}" )
    TB06.runMe()
    
    def fB07():
        A = pt.Product( [ pt.TermA(), pt.Sum( [ pt.Trace( pt.Product( [ pt.CoefficientFloat( -1.0 ), pt.Product( [ pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) ), pt.Trace( pt.Product( [ pt.CoefficientFloat( -3.0 ), pt.Product( [ pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) ) ] ) ] )
        B = pt.distributeAllTraces( A )
        return str( B )
    
    TB07 = UnitTest( "TB07: distributeAllTraces(), Distribution IV", fB07, " {A}  {-1.0}  {Trace[  {M}  {M} ]} +  {A}  {-3.0}  {Trace[  {M}  {M} ]}" )
    TB07.runMe()
    
    def fB08():
        A = pt.Product( [ pt.TermA(), pt.Sum( [ pt.Trace( pt.Product( [ pt.CoefficientFloat( -1.0 ), pt.Product( [ pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) ), pt.Trace( pt.Product( [ pt.CoefficientFloat( -3.0 ), pt.Product( [ pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) ) ] ), pt.TermA(), pt.Sum( [ pt.Trace( pt.Product( [ pt.CoefficientFloat( -5.0 ), pt.Product( [ pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) ), pt.Trace( pt.Product( [ pt.CoefficientFloat( -7.0 ), pt.Product( [ pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) ) ] ) ] )
        B = pt.distributeAllTraces( A )
        return str( B )
    
    TB08 = UnitTest( "TB08: distributeAllTraces(), Distribution V", fB08, " {A}  {-1.0}  {Trace[  {M}  {M} ]}  {A}  {-5.0}  {Trace[  {M}  {M} ]} +  {A}  {-1.0}  {Trace[  {M}  {M} ]}  {A}  {-7.0}  {Trace[  {M}  {M} ]} +  {A}  {-3.0}  {Trace[  {M}  {M} ]}  {A}  {-5.0}  {Trace[  {M}  {M} ]} +  {A}  {-3.0}  {Trace[  {M}  {M} ]}  {A}  {-7.0}  {Trace[  {M}  {M} ]}" )
    TB08.runMe()
    
    def fB09():
        A = pt.Product( [ pt.TermA(), pt.Trace( pt.Product( [ pt.CoefficientFloat( -1.0 ), pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ), pt.Trace( pt.Product( [ pt.CoefficientFloat( -3.0 ), pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ) ] )
        B = pt.distributeAllTraces( A )
        return str( B )
    
    TB09 = UnitTest( "TB09: distributeAllTraces(), Multiple Traces I", fB09, " {A}  {-1.0}  {Trace[  {M}  {M} ]}  {-3.0}  {Trace[  {M}  {M} ]}" )
    TB09.runMe()
    
    def fB10():
        M = pt.MatrixM( 10, 4, 1, True, "up" )
        detM = pt.DetM( M )
        D2detM = detM.derivative().derivative()
        A = pt.Product( [pt.TermA(),D2detM] )
        A.reduceTree()
        B = pt.distributeAllTraces( A )
        B.reduceTree()
        return str( B )
    
    TB10 = UnitTest( "TB10: distributeAllTraces(), Matrix and Determinants I", fB10, " {A}  {Det[ M_up ]}  {Trace[  {B_up}  {dM_up / dA} ]}  {Trace[  {B_up}  {dM_up / dA} ]} +  {A}  {Det[ M_up ]}  {-1.0}  {Trace[  {B_up}  {dM_up / dA}  {B_up}  {dM_up / dA} ]}" )
    TB10.runMe()
    
    def fB11():
        M = pt.MatrixM( 10, 4, 1, True, "up" )
        detM = pt.DetM( M )
        D1detM = detM.derivative()
        A = pt.Product( [pt.TermA(),D1detM] )
        A.reduceTree()
        B = pt.distributeAllTraces( A )
        B.reduceTree()
        return str( B )
    
    TB11 = UnitTest( "TB11: distributeAllTraces(), Matrix and Determinants II", fB11, " {A}  {Det[ M_up ]}  {Trace[  {B_up}  {dM_up / dA} ]}" )
    TB11.runMe()
    
    def fC01():
        A = []
        B = pt._elementsAreAllEqual( A )
        return str( B )
    
    TC01 = UnitTest( "TC01: _elementsAreAllEqual(), Trivial I", fC01, "True" )
    TC01.runMe()
    
    def fC02():
        A = [ "up" ]
        B = pt._elementsAreAllEqual( A )
        return str( B )
    
    TC02 = UnitTest( "TC02: _elementsAreAllEqual(), Trivial II", fC02, "True" )
    TC02.runMe()
    
    def fC03():
        A = [ "up", "up", "up", "up" ]
        B = pt._elementsAreAllEqual( A )
        return str( B )
    
    TC03 = UnitTest( "TC03: _elementsAreAllEqual(), Non-trivial I", fC03, "True" )
    TC03.runMe()
    
    def fC04():
        A = [ "up", "up", "dn", "up" ]
        B = pt._elementsAreAllEqual( A )
        return str( B )
    
    TC04 = UnitTest( "TC04: _elementsAreAllEqual(), Non-trivial II", fC04, "False" )
    TC04.runMe()
    
    def fC05():
        A = [ "up", "up", "up", "up", "dn" ]
        B = pt._elementsAreAllEqual( A )
        return str( B )
    
    TC05 = UnitTest( "TC05: _elementsAreAllEqual(), Non-trivial III", fC05, "False" )
    TC05.runMe()
    
    def fC06():
        A = [ "dn,","up", "up", "up", "up", "up" ]
        B = pt._elementsAreAllEqual( A )
        return str( B )
    
    TC06 = UnitTest( "TC06: _elementsAreAllEqual(), Non-trivial IV", fC06, "False" )
    TC06.runMe()
    
    def fC07():
        A = [ "up,","dn", "up", "dn", "dn", "up" ]
        B = pt._elementsAreAllEqual( A )
        return str( B )
    
    TC07 = UnitTest( "TC07: _elementsAreAllEqual(), Non-trivial V", fC07, "False" )
    TC07.runMe()
    
    def fD01():
        pt.MINIMAL_M_DISPLAY = True
        A = pt.MatrixM( 2, 2, 1, True, "up" )
        return str( A )
    
    TD01 = UnitTest( "TD01: MatrixM Construction I", fD01, "M_up")
    TD01.runMe()
    
    def fE01():
        A = pt.IndexedTrace( pt.Product( [ pt.MatrixM(1,1,1), pt.MatrixB(1,1,1) ] ) )
        B = pt.Product( [ pt.TermA(), A ] )
        C = pt._factorInvertedMatrices( B )
        return str( C )
        
    TE01 = UnitTest( "TE01: _factorInvertedMatrices I", fE01, " {A}  {B_(1, 0)}  {Gamma<1>_( 0 1 )[ { M_(0, 1) } ]}" )
    TE01.runMe()
    
    def fE02():
        A = pt.IndexedTrace( pt.Product( [ pt.MatrixM(1,1,1), pt.MatrixB(1,1,1), pt.MatrixM(1,1,1), pt.MatrixB(1,1,1), pt.MatrixM(1,1,1), pt.MatrixB(1,1,1) ] ) )
        B = pt.Product( [ pt.TermA(), A ] )
        C = pt._factorInvertedMatrices( B )
        return str( C )
        
    TE02 = UnitTest( "TE02: _factorInvertedMatrices II", fE02, " {A}  {B_(1, 2)}  {B_(3, 4)}  {B_(5, 0)}  {Gamma<3>_( 0 1 2 3 4 5 )[ { M_(0, 1) }{ M_(2, 3) }{ M_(4, 5) } ]}" )
    TE02.runMe()
    
    def fF01():
        A = pt.Trace( pt.Product( [ pt.MatrixM(1,1,1), pt.MatrixB(1,1,1) ] ) )
        B = pt.Trace( pt.Product( [ pt.MatrixM(1,1,1), pt.MatrixB(1,1,1) ] ) )
        C = pt.Trace( pt.Product( [ pt.MatrixM(1,1,1), pt.MatrixB(1,1,1) ] ) )
        D = pt.Sum( [pt.Product( [ A, B, C ] )] )
        E = pt.indexExpr( D )
        return str( E )
    
    TF01 = UnitTest( "TF01: indexExpr() I", fF01, " {Trace[ { M_(0, 1) }{ B_(1, 0) } ]}  {Trace[ { M_(2, 3) }{ B_(3, 2) } ]}  {Trace[ { M_(4, 5) }{ B_(5, 4) } ]}" )
    TF01.runMe()
    
    def fG01():
        A = pt.Sum( [ GenericTestTerm(1,0), GenericTestTerm(2,0) ] )
        B = pt.Sum( [ GenericTestTerm(3,0), GenericTestTerm(4,0) ] )
        C = pt.Product( [ A, B ] )
        D = C.getExpandedExpr()
        return str( D )
    
    TG01 = UnitTest( "TG01: getExpandedExpr() I", fG01, " {GT1_d0}  {GT3_d0} +  {GT1_d0}  {GT4_d0} +  {GT2_d0}  {GT3_d0} +  {GT2_d0}  {GT4_d0}" )
    TG01.runMe()
    
    def fG02():
        A = pt.Sum( [ GenericTestTerm(1,0), GenericTestTerm(2,0) ] )
        B = pt.Sum( [ GenericTestTerm(3,0), GenericTestTerm(4,0), GenericTestTerm(5,0) ] )
        C = pt.Product( [ GenericTestTerm(0,0), A, B ] )
        D = C.getExpandedExpr()
        D.reduceTree()
        return str( D )
    
    TG02 = UnitTest( "TG02: getExpandedExpr() II", fG02, " {GT0_d0}  {GT1_d0}  {GT3_d0} +  {GT0_d0}  {GT1_d0}  {GT4_d0} +  {GT0_d0}  {GT1_d0}  {GT5_d0} +  {GT0_d0}  {GT2_d0}  {GT3_d0} +  {GT0_d0}  {GT2_d0}  {GT4_d0} +  {GT0_d0}  {GT2_d0}  {GT5_d0}" )
    TG02.runMe()
    
    def fG03():
        A = pt.Sum( [ GenericTestTerm(2,0), GenericTestTerm(3,0) ] )
        B = pt.Product( [GenericTestTerm(1,0), A] )
        C = pt.Sum( [ GenericTestTerm(4,0), GenericTestTerm(5,0), GenericTestTerm(6,0) ] )
        D = pt.Product( [ B, C ] )
        D.reduceTree()
        E = D.getExpandedExpr()
        E.reduceTree()
        return str( E )
    
    TG03 = UnitTest( "TG03: getExpandedExpr() III", fG03, " {GT1_d0}  {GT2_d0}  {GT4_d0} +  {GT1_d0}  {GT2_d0}  {GT5_d0} +  {GT1_d0}  {GT2_d0}  {GT6_d0} +  {GT1_d0}  {GT3_d0}  {GT4_d0} +  {GT1_d0}  {GT3_d0}  {GT5_d0} +  {GT1_d0}  {GT3_d0}  {GT6_d0}" )
    TG03.runMe()
    
    def fG04():
        A = pt.Sum( [ GenericTestTerm(1,0), GenericTestTerm(2,0) ] )
        B = pt.Sum( [ GenericTestTerm(2,0), GenericTestTerm(3,0) ] )
        C = pt.Product( [ GenericTestTerm(4,0), A ] )
        D = pt.Product( [ GenericTestTerm(5,0), B])
        E = pt.Sum( [ GenericTestTerm(6,0), C ] )
        F = pt.Sum( [ GenericTestTerm(7,0), D ] )
        G = pt.Product( [ E, F ] )
        G.reduceTree()
        H = G.getExpandedExpr()
        H.reduceTree()
        return str( H )
    
    TG04 = UnitTest( "TG04: getExpandedExpr() IV", fG04, " {GT6_d0}  {GT7_d0} +  {GT6_d0}  {GT5_d0}  {GT2_d0} +  {GT6_d0}  {GT5_d0}  {GT3_d0} +  {GT4_d0}  {GT1_d0}  {GT7_d0} +  {GT4_d0}  {GT1_d0}  {GT5_d0}  {GT2_d0} +  {GT4_d0}  {GT1_d0}  {GT5_d0}  {GT3_d0} +  {GT4_d0}  {GT2_d0}  {GT7_d0} +  {GT4_d0}  {GT2_d0}  {GT5_d0}  {GT2_d0} +  {GT4_d0}  {GT2_d0}  {GT5_d0}  {GT3_d0}" )
    TG04.runMe()
    
    def fG05():
        A = pt.Sum( [ GenericTestTerm(1,0), GenericTestTerm(2,0) ] )
        B = pt.Sum( [ GenericTestTerm(2,0), GenericTestTerm(3,0) ] )
        C = pt.Product( [ GenericTestTerm(4,0), A ] )
        D = pt.Product( [ GenericTestTerm(5,0), B])
        E = pt.Sum( [ GenericTestTerm(6,0), C ] )
        F = pt.Sum( [ GenericTestTerm(7,0), D ] )
        G = pt.Product( [ E, F ] )
        G.reduceTree()
        
        H = pt.Sum( [ GenericTestTerm(1,0), GenericTestTerm(2,0) ] )
        I = pt.Sum( [ GenericTestTerm(2,0), GenericTestTerm(3,0) ] )
        J = pt.Product( [ GenericTestTerm(4,0), H ] )
        K = pt.Product( [ GenericTestTerm(5,0), I])
        L = pt.Sum( [ GenericTestTerm(6,0), J ] )
        M = pt.Sum( [ GenericTestTerm(7,0), K ] )
        N = pt.Product( [ L, M ] )
        N.reduceTree()
        
        O = pt.Sum( [ G, N ] )
        O.reduceTree()
        P = O.getExpandedExpr()
        P.reduceTree()
        return str( P )
    
    TG05 = UnitTest( "TG05: getExpandedExpr() V", fG05, " {GT6_d0}  {GT7_d0} +  {GT6_d0}  {GT5_d0}  {GT2_d0} +  {GT6_d0}  {GT5_d0}  {GT3_d0} +  {GT4_d0}  {GT1_d0}  {GT7_d0} +  {GT4_d0}  {GT1_d0}  {GT5_d0}  {GT2_d0} +  {GT4_d0}  {GT1_d0}  {GT5_d0}  {GT3_d0} +  {GT4_d0}  {GT2_d0}  {GT7_d0} +  {GT4_d0}  {GT2_d0}  {GT5_d0}  {GT2_d0} +  {GT4_d0}  {GT2_d0}  {GT5_d0}  {GT3_d0} +  {GT6_d0}  {GT7_d0} +  {GT6_d0}  {GT5_d0}  {GT2_d0} +  {GT6_d0}  {GT5_d0}  {GT3_d0} +  {GT4_d0}  {GT1_d0}  {GT7_d0} +  {GT4_d0}  {GT1_d0}  {GT5_d0}  {GT2_d0} +  {GT4_d0}  {GT1_d0}  {GT5_d0}  {GT3_d0} +  {GT4_d0}  {GT2_d0}  {GT7_d0} +  {GT4_d0}  {GT2_d0}  {GT5_d0}  {GT2_d0} +  {GT4_d0}  {GT2_d0}  {GT5_d0}  {GT3_d0}" )
    TG05.runMe()
    
    def fH01():
        A = pt.MatrixM(1,1,1).derivative()
        B = pt.Product( [ A, A ] )
        C = pt.Gamma( B, [ (0, 1), (2, 3) ] )
        D = C.getIntegratedTensorElement( [((1,1),(3,2)), ((4,1),(3,2))] )
        return str( D )
    
    TH01 = UnitTest( "TH01: Gamma, getIntegratedTensorElement(), Basic I", fH01, " {1 / 2}  {T_(1, 3)}  {T_(4, 3)}" )
    TH01.runMe()
    
    def fH02():
        A = pt.MatrixM(1,1,1).derivative()
        B = pt.Product( [ A, A ] )
        C = pt.Gamma( B, [ (0, 1), (2, 3) ] )
        D = C.getIntegratedTensorElement( [((1,1),(3,3)), ((4,1),(3,2))] )
        return str( D )
    
    TH02 = UnitTest( "TH02: Gamma, getIntegratedTensorElement(), Basic II", fH02, "0.0" )
    TH02.runMe()
    
    def fH03():
        A = pt.MatrixM(1,1,1).derivative()
        B = pt.Product( [ A, A ] )
        C = pt.Gamma( B, [ (0, 1), (2, 3) ] )
        D = C.getIntegratedTensorElement( [((1,1),(2,3)), ((4,3),(3,2))] )
        return str( D )
    
    TH03 = UnitTest( "TH03: Gamma, getIntegratedTensorElement(), Basic III", fH03, "0.0" )
    TH03.runMe()
    
    def fH04():
        A = pt.MatrixM(1,1,1).derivative()
        B = pt.Product( [ A, A, A ] )
        C = pt.Gamma( B, [ (0, 1), (2, 3), (3, 4) ] )
        D = C.getIntegratedTensorElement( [((1,1),(2,2)), ((1,2),(2,3)), ((5,2),(2,3))] )
        return str( D )
    
    TH04 = UnitTest( "TH05: Gamma, getIntegratedTensorElement(), Basic IV", fH04, "0.0" )
    TH04.runMe()
    
    def fH05():
        A = pt.MatrixM(1,1,1).derivative()
        B = pt.Product( [ A, A, A, A ] )
        C = pt.Gamma( B, [ (0, 1), (2, 3), (3, 4), (5, 6) ] )
        D = C.getIntegratedTensorElement( [((5,2),(2,3)), ((1,2),(2,3)), ((4,2),(2,3)), ((2,2),(2,3))] )
        return str( D )
    
    TH05 = UnitTest( "TH05: Gamma, getIntegratedTensorElement(), Intermediate I", fH05, " {3 / 8}  {T_(5, 2)}  {T_(1, 2)}  {T_(4, 2)}  {T_(2, 2)}" )
    TH05.runMe()
    
    def fH06():
        A = pt.MatrixM(1,1,1).derivative()
        B = pt.Product( [ A, A, A, A ] )
        C = pt.Gamma( B, [ (0, 1), (2, 3), (3, 4), (5, 6) ] )
        D = C.getIntegratedTensorElement( [((5,3),(2,4)), ((1,3),(2,4)), ((4,2),(2,3)), ((2,2),(2,3))] )
        return str( D )
    
    TH06 = UnitTest( "TH06: Gamma, getIntegratedTensorElement(), Intermediate II", fH06, " {1 / 2}  {1 / 2}  {T_(5, 2)}  {T_(1, 2)}  {T_(4, 2)}  {T_(2, 2)}" )
    TH06.runMe()
    
    def fH07():
        A = pt.MatrixM(1,1,1).derivative()
        B = pt.Product( [ A, A, A, A, A, A ] )
        C = pt.Gamma( B, [ (0, 1), (2, 3), (3, 4), (5, 6), (7, 8), (9, 10) ] )
        D = C.getIntegratedTensorElement( [((5,3),(2,4)), ((1,3),(2,4)), ((4,2),(2,3)), ((2,2),(2,3)), ((7,2),(2,3)), ((2,2),(2,3))] )
        return str( D )
    
    TH07 = UnitTest( "TH07: Gamma, getIntegratedTensorElement(), Intermediate III", fH07, " {3 / 8}  {1 / 2}  {T_(5, 2)}  {T_(1, 2)}  {T_(4, 2)}  {T_(2, 2)}  {T_(7, 2)}  {T_(2, 2)}" )
    TH07.runMe()
    
    def fH08():
        A = pt.MatrixM(1,1,1).derivative()
        B = pt.Product( [ A, A ] )
        C = pt.Gamma( B, [ (0, 1), (2, 3) ] )
        D = C.getIntegratedTensorElement( [((5,3),(1,4)), ((1,3),(2,4))], True )
        return str( D )
    
    TH08 = UnitTest( "TH08: Gamma, getIntegratedTensorElement(), Delta I", fH08, " {1 / 2}  {T_(5, 1)}  {T_(1, 2)}  {Delta<X>( 1, 2 )}" )
    TH08.runMe()
    
    def fH09():
        A = pt.MatrixM(1,1,1).derivative()
        B = pt.Product( [ A, A, A, A ] )
        C = pt.Gamma( B, [ (0, 1), (2, 3), (4, 5), (6, 7) ] )
        D = C.getIntegratedTensorElement( [((5,3),(1,4)), ((1,3),(2,4)), ((5,3),(5,4)), ((1,3),(8,4))], True )
        return str( D )
    
    TH09 = UnitTest( "TH09: Gamma, getIntegratedTensorElement(), Delta II", fH09, " {3 / 8}  {T_(5, 1)}  {T_(1, 2)}  {T_(5, 5)}  {T_(1, 8)}  {Delta<X>( 1, 2 )}  {Delta<X>( 2, 5 )}  {Delta<X>( 5, 8 )}" )
    TH09.runMe()
    
    def fH10():
        A = pt.MatrixM(1,1,1).derivative()
        B = pt.Product( [ A, A, A, A ] )
        C = pt.Gamma( B, [ (0, 1), (2, 3), (4, 5), (6, 7) ] )
        D = C.getIntegratedTensorElement( [((5,5),(1,6)), ((1,5),(2,6)), ((5,3),(5,4)), ((1,3),(8,4))], True )
        return str( D )
    
    TH10 = UnitTest( "TH10: Gamma, getIntegratedTensorElement(), Delta III", fH10, " {1 / 2}  {1 / 2}  {T_(5, 1)}  {T_(1, 2)}  {T_(5, 5)}  {T_(1, 8)}  {Delta<X>( 5, 8 )}  {Delta<X>( 1, 2 )}" )
    TH10.runMe()
    
    def fH11():
        A = pt.MatrixM(1,1,1).derivative()
        B = pt.Product( [ A, A, A, A ] )
        C = pt.Gamma( B, [ (0, 1), (2, 3), (4, 5), (6, 7) ] )
        D = C.getIntegratedTensorElement( [((5,5),(1,7)), ((1,5),(2,6)), ((5,3),(5,4)), ((1,3),(8,4))], True )
        return str( D )
    
    TH11 = UnitTest( "TH11: Gamma, getIntegratedTensorElement(), Delta IV", fH11, "0.0" )
    TH11.runMe()
    
    def fH12():
        A = pt.MatrixM(1,1,1).derivative()
        B = pt.Product( [ A, A, A, A ] )
        C = pt.Gamma( B, [ (0, 1), (2, 3), (4, 5), (6, 7) ] )
        D = C.getIntegratedTensorElement( [((5,5),(1,6)), ((1,5),(2,6)), ((5,5),(5,6)), ((1,3),(8,4))], True )
        return str( D )
    
    TH12 = UnitTest( "TH12: Gamma, getIntegratedTensorElement(), Delta V", fH12, "0.0" )
    TH12.runMe()
    
    def tI01():
        A = pt.calculateAllContractions( 2 )
        return str( A )
    
    TI01 = UnitTest( "TI01: calculateAllContractions(), n = 2", tI01, "[(2,)]" )
    TI01.runMe()
    
    def tI02():
        A = pt.calculateAllContractions( 4 )
        return str( A )
    
    TI02 = UnitTest( "TI02: calculateAllContractions(), n = 4", tI02, "[(4,), (2, 2)]" )
    TI02.runMe()
    
    def tI03():
        A = pt.calculateAllContractions( 6 )
        return str( A )
    
    TI03 = UnitTest( "TI03: calculateAllContractions(), n = 6", tI03, "[(6,), (2, 4), (2, 2, 2)]" )
    TI03.runMe()
    
    def tI04():
        A = pt.calculateAllContractions( 8 )
        return str( A )
    
    TI04 = UnitTest( "TI04: calculateAllContractions(), n = 8", tI04, "[(8,), (2, 6), (2, 2, 4), (2, 2, 2, 2)]" )
    TI04.runMe()
    
    def tI05():
        A = pt.calculateAllContractions( 10 )
        return str( A )
    
    TI05 = UnitTest( "TI05: calculateAllContractions(), n = 10", tI05, "[(10,), (2, 8), (2, 2, 6), (2, 2, 2, 4), (2, 2, 2, 2, 2)]" )
    TI05.runMe()

    def tI06():
        A = pt.getDeltaSignature( (2,) )
        return str( A )
    
    TI06 = UnitTest( "TI06: getDeltaSignature() I", tI06, "[[(0, 1)], []]" )
    TI06.runMe()
    
    def tI07():
        A = pt.getDeltaSignature( (4,) )
        return str( A )
    
    TI07 = UnitTest( "TI07: getDeltaSignature() II", tI07, "[[(0, 1), (1, 2), (2, 3)], []]" )
    TI07.runMe()
    
    def tI08():
        A = pt.getDeltaSignature( (2,2) )
        return str( A )
    
    TI08 = UnitTest( "TI08: getDeltaSignature() III", tI08, "[[(0, 1), (2, 3)], [(1, 2)]]" )
    TI08.runMe()
    
    def tI09():
        A = pt.getDeltaSignature( (6,) )
        return str( A )
    
    TI09 = UnitTest( "TI09: getDeltaSignature() IV", tI09, "[[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)], []]" )
    TI09.runMe()
    
    def tI10():
        A = pt.getDeltaSignature( (4,2) )
        return str( A )
    
    TI10 = UnitTest( "TI10: getDeltaSignature() V", tI10, "[[(0, 1), (1, 2), (2, 3), (2, 3)], [(3, 4)]]" )
    TI10.runMe()
    
    def tI11():
        A = pt.getDeltaSignature( (2,2,2) )
        return str( A )
    
    TI11 = UnitTest( "TI11: getDeltaSignature() VI", tI11, "[[(0, 1), (2, 3), (4, 5)], [(1, 2), (3, 4)]]" )
    TI11.runMe()
    
    def tI12():
        A = pt.generateCoordinateSpacePathIntegral( 2 )
        return str( A )
     
    TI12 = UnitTest( "TI12: generateCoordinateSpacePathIntegral(), n = 2", tI12, " {1 / 2}  { {Delta<XT>( 0, 1 )}}" )
    TI12.runMe()
    
    def tI13():
        A = pt.generateCoordinateSpacePathIntegral( 4 )
        return str( A )
     
    TI13 = UnitTest( "TI13: generateCoordinateSpacePathIntegral(), n = 4", tI13, " {3 / 8}  { {Delta<XT>( 0, 1 )}  {Delta<XT>( 1, 2 )}  {Delta<XT>( 2, 3 )}} +  {1 / 2}  {1 / 2}  { {Delta<XT>( 0, 1 )}  {Delta<XT>( 2, 3 )}  {DeltaBar<XT>( 1, 2 )} +  {Delta<XT>( 0, 1 )}  {Delta<XT>( 3, 2 )}  {DeltaBar<XT>( 1, 3 )} +  {Delta<XT>( 0, 2 )}  {Delta<XT>( 3, 1 )}  {DeltaBar<XT>( 2, 3 )}}" )
    TI13.runMe()
    
    def tI14():
        A = pt.generateCoordinateSpacePathIntegral( 6 )
        return str( A )
     
    TI14 = UnitTest( "TI14: generateCoordinateSpacePathIntegral(), n = 6", tI14, " {5 / 16}  { {Delta<XT>( 0, 1 )}  {Delta<XT>( 1, 2 )}  {Delta<XT>( 2, 3 )}  {Delta<XT>( 3, 4 )}  {Delta<XT>( 4, 5 )}} +  {1 / 2}  {3 / 8}  { {Delta<XT>( 0, 1 )}  {Delta<XT>( 2, 3 )}  {Delta<XT>( 3, 4 )}  {Delta<XT>( 4, 5 )}  {DeltaBar<XT>( 1, 2 )} +  {Delta<XT>( 0, 1 )}  {Delta<XT>( 3, 2 )}  {Delta<XT>( 2, 4 )}  {Delta<XT>( 4, 5 )}  {DeltaBar<XT>( 1, 3 )} +  {Delta<XT>( 0, 1 )}  {Delta<XT>( 4, 2 )}  {Delta<XT>( 2, 3 )}  {Delta<XT>( 3, 5 )}  {DeltaBar<XT>( 1, 4 )} +  {Delta<XT>( 0, 1 )}  {Delta<XT>( 5, 2 )}  {Delta<XT>( 2, 3 )}  {Delta<XT>( 3, 4 )}  {DeltaBar<XT>( 1, 5 )} +  {Delta<XT>( 0, 2 )}  {Delta<XT>( 3, 1 )}  {Delta<XT>( 1, 4 )}  {Delta<XT>( 4, 5 )}  {DeltaBar<XT>( 2, 3 )} +  {Delta<XT>( 0, 2 )}  {Delta<XT>( 4, 1 )}  {Delta<XT>( 1, 3 )}  {Delta<XT>( 3, 5 )}  {DeltaBar<XT>( 2, 4 )} +  {Delta<XT>( 0, 2 )}  {Delta<XT>( 5, 1 )}  {Delta<XT>( 1, 3 )}  {Delta<XT>( 3, 4 )}  {DeltaBar<XT>( 2, 5 )} +  {Delta<XT>( 0, 3 )}  {Delta<XT>( 4, 1 )}  {Delta<XT>( 1, 2 )}  {Delta<XT>( 2, 5 )}  {DeltaBar<XT>( 3, 4 )} +  {Delta<XT>( 0, 3 )}  {Delta<XT>( 5, 1 )}  {Delta<XT>( 1, 2 )}  {Delta<XT>( 2, 4 )}  {DeltaBar<XT>( 3, 5 )} +  {Delta<XT>( 0, 4 )}  {Delta<XT>( 5, 1 )}  {Delta<XT>( 1, 2 )}  {Delta<XT>( 2, 3 )}  {DeltaBar<XT>( 4, 5 )}} +  {1 / 2}  {1 / 2}  {1 / 2}  { {Delta<XT>( 0, 1 )}  {Delta<XT>( 2, 3 )}  {Delta<XT>( 4, 5 )}  {DeltaBar<XT>( 1, 2 )}  {DeltaBar<XT>( 3, 4 )} +  {Delta<XT>( 0, 1 )}  {Delta<XT>( 2, 3 )}  {Delta<XT>( 5, 4 )}  {DeltaBar<XT>( 1, 2 )}  {DeltaBar<XT>( 3, 5 )} +  {Delta<XT>( 0, 1 )}  {Delta<XT>( 2, 4 )}  {Delta<XT>( 5, 3 )}  {DeltaBar<XT>( 1, 2 )}  {DeltaBar<XT>( 4, 5 )} +  {Delta<XT>( 0, 1 )}  {Delta<XT>( 3, 4 )}  {Delta<XT>( 5, 2 )}  {DeltaBar<XT>( 1, 3 )}  {DeltaBar<XT>( 4, 5 )} +  {Delta<XT>( 0, 2 )}  {Delta<XT>( 3, 4 )}  {Delta<XT>( 5, 1 )}  {DeltaBar<XT>( 2, 3 )}  {DeltaBar<XT>( 4, 5 )}}" )
    TI14.runMe()
     
    def tI15():
        A = pt.getDeltaSignature( (2,) )
        B = pt.getIndexPermutations( A , 2 )
        return str( B )
     
    TI15 = UnitTest( "TI15: getIndexPermutations() I", tI15, "[[0, 1]]" )
    TI15.runMe()
     
    def tI16():
        A = pt.getDeltaSignature( (2,2) )
        B = pt.getIndexPermutations( A , 4 )
        return str( B )
     
    TI16 = UnitTest( "TI16: getIndexPermutations() II", tI16, "[[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 3, 1]]" )
    TI16.runMe()
    
    def tI17():
        A = pt.getDeltaSignature( (4,2) )
        B = pt.getIndexPermutations( A , 6 )
        return str( B )
     
    TI17 = UnitTest( "TI17: getIndexPermutations() III", tI17, "[[2, 3, 4, 0, 1, 5], [1, 3, 4, 0, 2, 5], [1, 2, 4, 0, 3, 5], [1, 2, 3, 0, 4, 5], [1, 2, 3, 0, 5, 4], [0, 3, 4, 1, 2, 5], [0, 2, 4, 1, 3, 5], [0, 2, 3, 1, 4, 5], [0, 2, 3, 1, 5, 4], [0, 1, 4, 2, 3, 5], [0, 1, 3, 2, 4, 5], [0, 1, 3, 2, 5, 4], [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 5, 4], [0, 1, 2, 4, 5, 3]]" )
    TI17.runMe()
    
    def tI18():
        A = pt.getDeltaSignature( (2,2,2) )
        B = pt.getIndexPermutations( A , 6 )
        return str( B )
     
    TI18 = UnitTest( "TI18: getIndexPermutations() II", tI18, "[[0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 5, 4], [0, 1, 2, 4, 5, 3], [0, 1, 3, 4, 5, 2], [0, 2, 3, 4, 5, 1]]" )
    TI18.runMe()
    
    print "\n-------------------------------------------------"
    print "Tests PASSED: " + str( UnitTest.testsPassed ) + " , tests FAILED: " + str( UnitTest.testsFailed ) + "."
executeAllUnitTests()