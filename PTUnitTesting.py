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
            
def executeUnitTests():
    def tA01():
        A = pt.MatrixM()
        return str( A )
    
    TA01 = UnitTest( "TA01: MatrixM, Default Construction", tA01, "M" )
    TA01.runMe()
    
    def tA02():
        A = pt.MatrixM( False )
        return str( A )
    
    TA02 = UnitTest( "TA02: MatrixM, Non-interacting", tA02, "M0" )
    TA02.runMe()
    
    def tA03():
        A = pt.MatrixM( True, "up" )
        return str( A )
    
    TA03 = UnitTest( "TA03: MatrixM, Flavor label", tA03, "M_up" )
    TA03.runMe()
    
    def tA04():
        A = pt.MatrixM( False, "up" )
        return str( A )
    
    TA04 = UnitTest( "TA04: MatrixM, Flavor label, non-interacting", tA04, "M0_up" )
    TA04.runMe()
    
    def tA05():
        A = pt.MatrixM( True, "up" )
        B = A.derivative()
        return str( B )
    
    TA05 = UnitTest( "TA05: MatrixM, First derivative (interacting)", tA05, "dM_up / dA" )
    TA05.runMe()
    
    def tA06():
        A = pt.MatrixM( True, "up" )
        B = A.derivative().derivative()
        return str( B )
    
    TA06 = UnitTest( "TA06: MatrixM, Second derivative (interacting)", tA06, "0.0" )
    TA06.runMe()
    
    def tB01():
        A = pt.MatrixB()
        return str( A )
    
    TB01 = UnitTest( "TB01: MatrixB, Default construction", tB01, "B" )
    TB01.runMe()
    
    def tB02():
        A = pt.MatrixB()
        A.setAsNoninteracting()
        return str( A )
    
    TB02 = UnitTest( "TB02: MatrixB, setAsNonInteracting()", tB02, "B0" )
    TB02.runMe()
    
    def tB03():
        A = pt.MatrixB()
        A = A.derivative()
        return str( A )
    
    TB03 = UnitTest( "TB03: MatrixB, setAsNonInteracting()", tB03, " {-1.0}  {B}  {dM / dA}  {B}" )
    TB03.runMe()
    
    def tB04():
        A = pt.MatrixB( True, "up" )
        B = pt.MatrixB( True, "dn" )
        return str( A == B )
    
    TB04 = UnitTest( "TB04: MatrixB, __eq__ I", tB04, "False" )
    TB04.runMe()
    
    def tB05():
        A = pt.MatrixB( True, "up" )
        B = pt.MatrixB( True, "up" )
        return str( A == B )
    
    TB05 = UnitTest( "TB05: MatrixB, __eq__ II", tB05, "True" )
    TB05.runMe()
    
    def tB06():
        A = pt.MatrixB( True, "up" )
        B = pt.MatrixB( True, "up" )
        return str( A == A )
    
    TB06 = UnitTest( "TB06: MatrixB, __eq__ III", tB06, "True" )
    TB06.runMe()
        
    def tB07():
        A = pt.MatrixB( True, "up" )
        A = A.derivative()
        return str( A )
    
    TB07 = UnitTest( "TB07: MatrixB, First derivative (interacting)", tB07, " {-1.0}  {B_up}  {dM_up / dA}  {B_up}" )
    TB07.runMe()
    
    def tB08():
        A = pt.MatrixB( True, "up" )
        A = A.derivative().derivative()
        return str( A )
    
    TB08 = UnitTest( "TB08: MatrixB, Second derivative (interacting)", tB08, " {-1.0}  { {-1.0}  {B_up}  {dM_up / dA}  {B_up}}  {dM_up / dA}  {B_up} +  {-1.0}  {B_up}  {dM_up / dA}  { {-1.0}  {B_up}  {dM_up / dA}  {B_up}}" )
    TB08.runMe()
    
    def tB09():
        A = pt.MatrixB( True, "up" )
        A = A.derivative().derivative().derivative()
        A.reduceTree()
        return str( A )
    
    TB09 = UnitTest( "TB09: MatrixB, Third derivative (interacting)", tB09, "" )
    TB09.runMe()
    
    def tB10():
        A = pt.MatrixB( False, "up" )
        A = A.derivative()
        return str( A )
    
    TB10 = UnitTest( "TB10: MatrixB, First derivative (non-interacting)", tB10, "0.0" )
    TB10.runMe()
    
    def tC01():
        A = pt.MatrixM()
        B = pt.DetM( A )
        return str( B )
    
    TC01 = UnitTest( "TC01: DetM, Default construction", tC01, "Det[ M ]" )
    TC01.runMe()
    
    def tC02():
        A = pt.MatrixM()
        B = pt.DetM( A, True, True )
        return str( B )
    
    TC02 = UnitTest( "TC02: DetM, Inverted", tC02, "1 / Det[ M ]" )
    TC02.runMe()
    
    def tC03():
        A = pt.MatrixM()
        B = pt.DetM( A )
        C = B.derivative()
        return str( C )
    
    TC03 = UnitTest( "TC03: DetM, First derivative", tC03, " {Det[ M ]}  {Trace[  {B}  {dM / dA} ]}" )
    TC03.runMe()
    
    def tC04():
        A = pt.MatrixM()
        B = pt.DetM( A )
        C = B.derivative().derivative()
        C.reduceTree()
        return str( C )
    
    TC04 = UnitTest( "TC04: DetM, Second derivative", tC04, " {Det[ M ]}  {Trace[  {B}  {dM / dA} ]}  {Trace[  {B}  {dM / dA} ]} +  {Det[ M ]}  {Trace[  {-1.0}  {B}  {dM / dA}  {B}  {dM / dA} ]}" )
    TC04.runMe()
    
    def tD01():
        A = pt.MatrixK()
        return str(A)
    
    TD01 = UnitTest( "TD01: MatrixK, Default constructor", tD01, "" )
    TD01.runMe()
    
    def tZ01():
        Mup = pt.MatrixM( False, "up" )
        Mdn = pt.MatrixM( False, "dn" )
        A = pt.Product( [pt.DetM( Mup ), pt.DetM( Mdn )] )
        B = pt.areTermsCommon(A, A)
        return str( B )
    
    TZ01 = UnitTest( "TZ01: areTermsCommon I", tZ01, "True" )
    TZ01.runMe()
    
    def tZ02():
        Mup = pt.MatrixM( False, "up" )
        Mdn = pt.MatrixM( False, "dn" )
        A = pt.Product( [pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn )] )
        B = pt.areTermsCommon(A, A)
        return str( B )
    
    TZ02 = UnitTest( "TZ02: areTermsCommon II", tZ02, "True" )
    TZ02.runMe()
    
    def tZ03():
        Mup = pt.MatrixM( False, "up" )
        Mdn = pt.MatrixM( False, "dn" )
        A = pt.Product( [pt.CoefficientFloat( 1.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn )] )
        B = pt.Product( [pt.CoefficientFloat( 2.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn )] )
        C = pt.areTermsCommon(A, B)
        return str( C )
    
    TZ03 = UnitTest( "TZ03: areTermsCommon III", tZ03, "True" )
    TZ03.runMe()
    
    def tZ04():
        Mup = pt.MatrixM( False, "up" )
        Mdn = pt.MatrixM( False, "dn" )
        A = pt.Product( [pt.CoefficientFloat( 1.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn )] )
        B = pt.Product( [pt.CoefficientFloat( 2.0 ), pt.TermA(), pt.DetM( Mdn ), pt.DetM( Mdn )] )
        C = pt.areTermsCommon(A, B)
        return str( C )
    
    TZ04 = UnitTest( "TZ04: areTermsCommon IV", tZ04, "False" )
    TZ04.runMe()
    
    def tZ05():
        Mup = pt.MatrixM( False, "up" )
        Mdn = pt.MatrixM( False, "dn" )
        A = pt.Product( [pt.CoefficientFloat( 1.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn )] )
        B = pt.Product( [pt.CoefficientFloat( 2.0 ), pt.TermA(), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn )] )
        C = pt.areTermsCommon(A, B)
        return str( C )
    
    TZ05 = UnitTest( "TZ05: areTermsCommon V", tZ05, "False" )
    TZ05.runMe()
    
    def tZ06():
        Mup = pt.MatrixM( False, "up" )
        Mdn = pt.MatrixM( False, "dn" )
        Dup = pt.MatrixK( "up" )
        Dup.fourierTransform()
        Ddn = pt.MatrixK( "dn" )
        Ddn.fourierTransform()
        A = pt.Product( [pt.CoefficientFloat( 1.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Dup ] )
        B = pt.Product( [pt.CoefficientFloat( 2.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Ddn ] )
        C = pt.areTermsCommon(A, B)
        return str( C )
    
    TZ06 = UnitTest( "TZ06: areTermsCommon VI", tZ06, "False" )
    TZ06.runMe()
    
    def tZ07():
        Mup = pt.MatrixM( False, "up" )
        Mdn = pt.MatrixM( False, "dn" )
        Dup = pt.MatrixK( "up" )
        Dup.fourierTransform()
        Ddn = pt.MatrixK( "dn" )
        Ddn.fourierTransform()
        A = pt.Product( [pt.CoefficientFloat( 1.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Ddn, Dup ] )
        B = pt.Product( [pt.CoefficientFloat( 2.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Dup, Ddn ] )
        C = pt.areTermsCommon(A, B)
        return str( C )
    
    TZ07 = UnitTest( "TZ07: areTermsCommon VII", tZ07, "True" )
    TZ07.runMe()
    
    def tZ08():
        Mup = pt.MatrixM( False, "up" )
        Mdn = pt.MatrixM( False, "dn" )
        Dup = pt.MatrixK( "up" )
        Dup.fourierTransform()
        Ddn = pt.MatrixK( "dn" )
        Ddn.fourierTransform()
        F1 = pt.FourierSum( [0,1,0,1], 4 )
        F2 = pt.FourierSum( [0,0,1,1], 4 )
        A = pt.Product( [pt.CoefficientFloat( 1.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Ddn, Dup, F1 ] )
        B = pt.Product( [pt.CoefficientFloat( 2.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Dup, Ddn, F2 ] )
        C = pt.areTermsCommon(A, B)
        return str( C )
    
    TZ08 = UnitTest( "TZ08: areTermsCommon VIII", tZ08, "False" )
    TZ08.runMe()
    
    def tY01():
        Mup = pt.MatrixM( False, "up" )
        Mdn = pt.MatrixM( False, "dn" )
        Dup = pt.MatrixK( "up" )
        Dup.fourierTransform()
        Dup.indices = 0
        Ddn = pt.MatrixK( "dn" )
        Ddn.fourierTransform()
        Ddn.indices = 0
        A = pt.Product( [pt.CoefficientFloat( 1.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Ddn, Dup ] )
        B = pt.Product( [pt.CoefficientFloat( 2.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Dup, Ddn ] )
        C = pt.Sum( [ A, B ] )
        D = pt.combineLikeTerms( C )
        D.simplify()
        return str( D )
    
    TY01 = UnitTest( "TY01: combineLikeTerms() I", tY01, " {A}  {Det[ M_up ]}  {Det[ M_dn ]}  {D_dn_0}  {D_up_0}  {3.0}" )
    TY01.runMe()
    
    def tY02():
        Mup = pt.MatrixM( False, "up" )
        Mdn = pt.MatrixM( False, "dn" )
        Dup = pt.MatrixK( "up" )
        Dup.fourierTransform()
        Dup.indices = 0
        Ddn = pt.MatrixK( "dn" )
        Ddn.fourierTransform()
        Ddn.indices = 0
        A = pt.Product( [pt.CoefficientFloat( 1.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Ddn, Dup ] )
        B = pt.Product( [pt.CoefficientFloat( 2.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Dup, Ddn ] )
        C = pt.Product( [pt.CoefficientFloat( 5.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Dup, Dup ] )
        D = pt.Sum( [ A, B, C ] )
        E = pt.combineLikeTerms( D )
        E.simplify()
        return str( E )
    
    TY02 = UnitTest( "TY02: combineLikeTerms() II", tY02, " {A}  {Det[ M_up ]}  {Det[ M_dn ]}  {D_dn_0}  {D_up_0}  {3.0} +  {A}  {Det[ M_up ]}  {Det[ M_dn ]}  {D_up_0}  {D_up_0}  {5.0}" )
    TY02.runMe()
    
    def tY03():
        Mup = pt.MatrixM( False, "up" )
        Mdn = pt.MatrixM( False, "dn" )
        Dup = pt.MatrixK( "up" )
        Dup.fourierTransform()
        Dup.indices = 0
        Ddn = pt.MatrixK( "dn" )
        Ddn.fourierTransform()
        Ddn.indices = 0
        A1 = pt.Product( [pt.CoefficientFloat( 1.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Ddn, Dup ] )
        A2 = pt.Product( [pt.CoefficientFloat( 9.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Dup, Dup ] )
        A3 = pt.Product( [pt.CoefficientFloat( 2.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Dup, Ddn ] )
        A4 = pt.Product( [pt.CoefficientFloat( 5.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Dup, Dup ] )
        D = pt.Sum( [ A1, A2, A3, A4 ] )
        E = pt.combineLikeTerms( D )
        E.simplify()
        return str( E )
    
    TY03 = UnitTest( "TY02: combineLikeTerms() III", tY03, " {A}  {Det[ M_up ]}  {Det[ M_dn ]}  {D_dn_0}  {D_up_0}  {3.0} +  {A}  {Det[ M_up ]}  {Det[ M_dn ]}  {D_up_0}  {D_up_0}  {14.0}" )
    TY03.runMe()
    
    def tY04():
        Mup = pt.MatrixM( False, "up" )
        Mdn = pt.MatrixM( False, "dn" )
        Dup = pt.MatrixK( "up" )
        Dup.fourierTransform()
        Dup.indices = 0
        Ddn = pt.MatrixK( "dn" )
        Ddn.fourierTransform()
        Ddn.indices = 0
        A1 = pt.Product( [pt.CoefficientFloat( 1.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Ddn, Dup ] )
        A2 = pt.Product( [pt.CoefficientFloat( 9.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Dup, Dup ] )
        A3 = pt.Product( [pt.CoefficientFloat( 2.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Dup, Ddn ] )
        A4 = pt.Product( [pt.CoefficientFloat( 5.0 ), pt.TermA(), pt.DetM( Mup ), pt.DetM( Mdn ), Dup, Dup ] )
        D = pt.Sum( [ A3, A2, A4, A1 ] )
        E = pt.combineLikeTerms( D )
        E.simplify()
        return str( E )
    
    TY04 = UnitTest( "TY02: combineLikeTerms() IV", tY04, " {A}  {Det[ M_up ]}  {Det[ M_dn ]}  {D_up_0}  {D_dn_0}  {3.0} +  {A}  {Det[ M_up ]}  {Det[ M_dn ]}  {D_up_0}  {D_up_0}  {14.0}" )
    TY04.runMe()
    
    def tX01():
        A = [(0,1),(2,3),(0,2)]
        B = pt._constructContractionDict( A )
        return str( B )
        
    TX01 = UnitTest( "TX01: _constructContractionDict() I", tX01, "{0: 0, 1: 0, 2: 0, 3: 0}" )
    TX01.runMe()
    
    def tX02():
        A = [(0,2),(2,3),(0,1)]
        B = pt._constructContractionDict( A )
        return str( B )
        
    TX02 = UnitTest( "TX02: _constructContractionDict() II", tX02, "{0: 0, 1: 0, 2: 0, 3: 0}" )
    TX02.runMe()
    
    def tX03():
        A = [(1,2),(0,1),(0,2)]
        B = pt._constructContractionDict( A )
        return str( B )
        
    TX03 = UnitTest( "TX03: _constructContractionDict() III", tX03, "{0: 0, 1: 0, 2: 0}" )
    TX03.runMe()
    
    def tX04():
        A = [(0,2),(0,3),(1,2)]
        B = pt._constructContractionDict( A )
        return str( B )
        
    TX04 = UnitTest( "TX04: _constructContractionDict() IV", tX04, "{0: 0, 1: 0, 2: 0, 3: 0}" )
    TX04.runMe()
    
    def tX05():
        A = [(1,2),(0,3),(0,2)]
        B = pt._constructContractionDict( A )
        return str( B )
        
    TX05 = UnitTest( "TX05: _constructContractionDict() V", tX05, "{0: 0, 1: 0, 2: 0, 3: 0}" )
    TX05.runMe()
    
    def tX06():
        A = [(0,1),(2,3),(4,5),(6,7),(0,2),(2,4),(4,6)]
        B = pt._constructContractionDict( A )
        return str( B )
        
    TX06 = UnitTest( "TX06: _constructContractionDict() VI", tX06, "{0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0}" )
    TX06.runMe()
    
    def tX07():
        A = [(1,2),(3,4),(0,5),(0,2),(4,6)]
        B = pt._constructContractionDict( A )
        return str( B )
        
    TX07 = UnitTest( "TX07: _constructContractionDict() VII", tX07, "{0: 0, 1: 0, 2: 0, 3: 3, 4: 3, 5: 0, 6: 3}" )
    TX07.runMe()
    
executeUnitTests()
print "\n-------------------------------------------------"
print "Tests PASSED: " + str( UnitTest.testsPassed ) + " , tests FAILED: " + str( UnitTest.testsFailed ) + "."