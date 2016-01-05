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
        print str(A)
        B = pt.distributeAllTraces( A )
        return str( B )
    
    TB05 = UnitTest( "TB05: distributeAllTraces(), Distribution II", fB05, " {-1.0}  {Trace[  {M}  {M} ]}  {-3.0}  {Trace[  {M}  {M} ]}" )
    TB05.runMe()
    
    def fB06():
        A = pt.Sum( [ pt.Trace( pt.Product( [ pt.CoefficientFloat( -1.0 ), pt.Product( [ pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) ), pt.Trace( pt.Product( [ pt.CoefficientFloat( -3.0 ), pt.Product( [ pt.MatrixM( 1, 1, 1 ), pt.MatrixM( 1, 1, 1 ) ] ) ] ) ) ] )
        print str(A)
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
    
    print "\n-------------------------------------------------"
    print "Tests PASSED: " + str( UnitTest.testsPassed ) + " , tests FAILED: " + str( UnitTest.testsFailed ) + "."
executeAllUnitTests()