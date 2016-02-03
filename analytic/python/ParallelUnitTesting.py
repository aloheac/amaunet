import PTSymbolicObjects as pt
import PTMultithreading as ptmt
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
    
    def pA01():
        A = pt.Sum( [GenericTestTerm(0,0), GenericTestTerm(1,0), GenericTestTerm(2,0), GenericTestTerm(3,0)] )
        B = pt.Sum( [GenericTestTerm(4,0), GenericTestTerm(5,0), GenericTestTerm(6,0), GenericTestTerm(7,0)] )
        C = ptmt.execParallelDualExpansion( A, B )
        return str( C )

    PA01 = UnitTest( "PA01: MT execParallelDualExpansion() I", pA01, "" )
    PA01.runMe()

executeUnitTests()
print "\n-------------------------------------------------"
print "Tests PASSED: " + str( UnitTest.testsPassed ) + " , tests FAILED: " + str( UnitTest.testsFailed ) + "."