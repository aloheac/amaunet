/*
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * High-order Perturbation Theory Analytics
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Symbolic Peturbation Theory Expression Objects Source Implementation
 *
 * v. 0.1		02 Feb 2016
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at CHapel Hill
 */

#include "PTSymbolicObjects.h"
#include <sstream>
#include <cstring>
#include <boost/algorithm/string/trim.hpp>

using namespace std;

/*
 * ***********************************************************************
 * EXPRESSION BASE CLASSES : IMPLEMENTATIONS
 * ***********************************************************************
 */

/*
 *  SymbolicTerm
 */

SymbolicTerm::SymbolicTerm() {
	derivativeOrder = 0;
	isInteracting = true;
	flavorLabel = "";
	//indices = new int[ 2 ];
	termID = '0';
}


SymbolicTerm::SymbolicTerm( const SymbolicTerm* st ) {
	cout << "***** Pointer copy constructor called.";
}

SymbolicTerm::~SymbolicTerm() {
	//delete[] indices;
}

bool SymbolicTerm::operator==( const SymbolicTerm &other ) const {
	return ( derivativeOrder == other.derivativeOrder ) and
			( isInteracting == other.isInteracting ) and
			( flavorLabel == other.flavorLabel );
}

ostream& operator<<( ostream& os, const SymbolicTerm &st ) {
	stringstream ss;
	ss << st.to_string();
	os << ss.str();
	return os;
}

Sum SymbolicTerm::getDerivative() {
	return Sum( new CoefficientFloat( 0.0 ) );
}

void SymbolicTerm::simplify() { }

const string SymbolicTerm::to_string() const {
	return "<invalid_term>";
}

void SymbolicTerm::setAsNonInteracting() {
	isInteracting = false;
}

bool SymbolicTerm::isTermInteracting() {
	return isInteracting;
}

char* SymbolicTerm::getFlavorLabel() {
	return flavorLabel;
}

int* SymbolicTerm::getIndices() {
	return indices;
}

void SymbolicTerm::setIndices( int* newIndices ) {
	indices[0] = newIndices[0];
	indices[1] = newIndices[1];
}

char SymbolicTerm::getTermID() {
	return termID;
}

void SymbolicTerm::reduceTree() { }

SymbolicTerm* SymbolicTerm::copy() {
	SymbolicTerm* copy;

	if ( termID == 'S' ) {
		copy = new Sum( dynamic_cast<Sum*>( this ) );
	} else if ( termID == 'P' ) {
		copy = new Product( dynamic_cast<Product*>( this ) );
	} else if ( termID == 'K' ) {
		copy = new MatrixK( dynamic_cast<MatrixK*>( this ) );
	} else if ( termID == 's' ) {
		copy = new MatrixS( dynamic_cast<MatrixS*>( this ) );
	} else if ( termID == 'B' ) {
		copy = new MatrixB( dynamic_cast<MatrixB*>( this ) );
	} else if ( termID == 'A' ) {
		copy = new TermA( dynamic_cast<TermA*>( this ) );
	} else if ( termID == 'T' ) {
		copy = new Trace( dynamic_cast<Trace*>( this ) );
	} else if ( termID == 'L' ) {
		copy = new CoefficientFloat( dynamic_cast<CoefficientFloat*>( this ) );
	} else if ( termID == 'R' ) {
		copy = new CoefficientFraction( dynamic_cast<CoefficientFraction*>( this ) );
	} else if ( termID == 'M' ) {
		copy = dynamic_cast<MatrixM*>( this )->copy();
	} else if ( termID == 'g' ) {
		copy = dynamic_cast<GenericTestTerm*>( this )->copy();
	} else {
		cout << "***ERROR: Invalid SymbolicTerm* copy requested for term ID '" << termID << "'." << endl;
	}

	return copy;
}

/*
 * ***********************************************************************
 * TERM DERIVED CLASSES - SCALAR AND MATRIX OBJECTS : IMPLEMENTATIONS
 * ***********************************************************************
 */

/*
 * GenericTestTerm
 */

GenericTestTerm::GenericTestTerm( int thisId, int thisDerivativeOrder ) : SymbolicTerm() {
	id = thisId;
	derivativeOrder = thisDerivativeOrder;
	termID = 'g';
}

const string GenericTestTerm::to_string() const {
	stringstream ss;
	ss << "GT_" << id << "^" << derivativeOrder;
	return ss.str();
}

Sum GenericTestTerm::getDerivative() {
	return Sum( new GenericTestTerm( id, derivativeOrder + 1 ) );
}

GenericTestTerm* GenericTestTerm::copy() {
	GenericTestTerm* cpy = new GenericTestTerm( 0, 0 );
	cpy->id = id;
	cpy->derivativeOrder = derivativeOrder;
	cpy->termID = 'g';

	return cpy;
}

std::ostream& operator<<( std::ostream& os, const GenericTestTerm &st ) {
	os << st.to_string();
	return os;
}

/*
 * MatrixM
 */

MatrixM::MatrixM() : SymbolicTerm() {
	termID = 'M';
};

MatrixM::MatrixM( char* thisFlavorLabel ) : SymbolicTerm() {
	termID = 'M';
	flavorLabel = thisFlavorLabel;
}

bool MatrixM::operator==( const MatrixM &other ) const {
	return ( derivativeOrder == other.derivativeOrder ) and
				( isInteracting == other.isInteracting ) and
				( flavorLabel == other.flavorLabel ) and
				other.termID == 'M';
}

Sum MatrixM::getDerivative() {
    MatrixM* D = new MatrixM( flavorLabel );
    D->derivativeOrder = derivativeOrder + 1;
    if ( !isInteracting ) {
    	D->setAsNonInteracting();
    }

    vector<SymbolicTerm*> vecD = vector<SymbolicTerm*>();
    vecD.push_back( D );

    return Sum( vecD );
}

const string MatrixM::to_string() const {
	stringstream ss;
	if ( strcmp( flavorLabel, "" ) != 0 ) {
		if ( isInteracting ) {
				if ( derivativeOrder == 0 ) {
					ss << "M_" << flavorLabel;
				} else if ( derivativeOrder == 1 ) {
					ss << "dM_" << flavorLabel << " / dA";
				} else {
					// If placing a matrix that does not vanish for derivativeOrder > 1,
					// changes must be inserted here.
					ss << "0";
				}
		} else {  // M is non-interacting.
			if ( derivativeOrder == 0 ) {
				ss << "M0_" << flavorLabel;
			} else {
				ss << "0";
			}
		}
	} else {
		if ( isInteracting ) {
			if ( derivativeOrder == 0 ) {
				ss << "M";
			} else if ( derivativeOrder == 1 ) {
				ss << "dM / dA";
			} else {
				// If placing a matrix that does not vanish for derivativeOrder > 1,
				// changes must be inserted here.
				ss << "0.0";
			}
		} else {
			if ( derivativeOrder == 0 ) {
				ss << "M0";
			} else {
				ss << "0";
			}
		}
	}

	return ss.str();
}

MatrixM* MatrixM::copy() {
	MatrixM* cpy = new MatrixM();
	cpy->derivativeOrder = derivativeOrder;
	cpy->flavorLabel = flavorLabel;
	cpy->isInteracting = isInteracting;
	cpy->indices[0] = indices[0];
	cpy->indices[1] = indices[1];
	cpy->termID = 'M';

	return cpy;
}

/*
 * MatrixB
 */

MatrixB::MatrixB() : SymbolicTerm() {
	termID = 'B';
}

MatrixB::MatrixB( char* thisFlavorLabel )  : SymbolicTerm() {
	flavorLabel = thisFlavorLabel;
	termID = 'B';
}

MatrixB::MatrixB( MatrixB* b ) : SymbolicTerm() {
	derivativeOrder = b->derivativeOrder;
	isInteracting = b->isInteracting;
	flavorLabel = b->flavorLabel;
	termID = 'B';
}

bool MatrixB::operator==( const MatrixB &other ) const {
	return ( derivativeOrder == other.derivativeOrder ) and
			( isInteracting == other.isInteracting ) and
			( flavorLabel == other.flavorLabel ) and
			( other.termID == 'B' );
}

Sum MatrixB::getDerivative() {
	if ( isInteracting ) {
		MatrixM* dM = new MatrixM( flavorLabel );
		dM->derivativeOrder = derivativeOrder + 1;

		vector<SymbolicTerm*> derivativeTerms = vector<SymbolicTerm*>();
		derivativeTerms.push_back( new CoefficientFloat( -1.0 ) );
		derivativeTerms.push_back( new MatrixB( flavorLabel ) );
		derivativeTerms.push_back( dM );
		derivativeTerms.push_back( new MatrixB( flavorLabel ) );

		vector<SymbolicTerm*> productTerms = vector<SymbolicTerm*>();
		productTerms.push_back( new Product( derivativeTerms ) );

		return Sum( productTerms );
	} else {
		vector<SymbolicTerm*> derivativeTerms = vector<SymbolicTerm*>();
		derivativeTerms.push_back( new CoefficientFloat( 0.0 ) );
		return Sum( derivativeTerms );
	}
}

const std::string MatrixB::to_string() const {
	stringstream ss;
	if ( strcmp( flavorLabel, "" ) == 0 ) {
		if ( isInteracting ) {
			if ( derivativeOrder == 0 ) {
				ss << "B";
			} else if ( derivativeOrder == 1 ) {
				ss << "dB / dA";
			} else {
				ss << "d^" << derivativeOrder << " B / dA^" << derivativeOrder;
			}
		} else {
			if ( derivativeOrder == 0 ) {
				ss << "B0";
			} else {
				ss << "0.0";
			}
		}
	} else {
		if ( isInteracting ) {
			if ( derivativeOrder == 0 ) {
				ss << "B_" << flavorLabel;
			} else if ( derivativeOrder == 1 ) {
				ss << "dB_" << flavorLabel << " / dA";
			} else {
				ss << "d^" << derivativeOrder << " B_" << flavorLabel << " / dA^" << derivativeOrder;
			}
		} else {
			if ( derivativeOrder == 0 ) {
				ss << "B0_" << flavorLabel;
			} else {
				ss << "0.0";
			}
		}
	}

	return ss.str();
}

MatrixB* MatrixB::copy() {
	MatrixB* cpy = new MatrixB();
	cpy->derivativeOrder = derivativeOrder;
	cpy->flavorLabel = flavorLabel;
	cpy->isInteracting = isInteracting;
	cpy->indices[0] = indices[0];
	cpy->indices[1] = indices[1];
	cpy->termID = 'B';

	return cpy;
}

/*
 * MatrixK
 */

MatrixK::MatrixK() : SymbolicTerm() {
	isFourierTransformed = false;
	termID = 'K';
}

MatrixK::MatrixK( char* thisFlavorLabel ) : SymbolicTerm() {
	isFourierTransformed = false;
	flavorLabel = thisFlavorLabel;
	termID = 'K';
}

MatrixK::MatrixK( const MatrixK* k ) : SymbolicTerm() {
	isFourierTransformed = k->isFourierTransformed;
	flavorLabel = k->flavorLabel;
	termID = 'K';
}

bool MatrixK::operator==( const MatrixK &other ) const {
	return ( isFourierTransformed == other.isFourierTransformed ) and ( flavorLabel == other.flavorLabel );
}

const string MatrixK::to_string() const {
	stringstream ss;

	if ( isFourierTransformed ) {
		ss << "D_" << flavorLabel << "_" << indices;
	} else {
		ss << "K_" << flavorLabel << "_" << indices;
	}

	return ss.str();
}

MatrixK* MatrixK::copy() {
	MatrixK* cpy = new MatrixK();
	cpy->derivativeOrder = derivativeOrder;
	cpy->flavorLabel = flavorLabel;
	cpy->isInteracting = isInteracting;
	cpy->indices[0] = indices[0];
	cpy->indices[1] = indices[1];
	cpy->termID = 'K';

	return cpy;
}

void MatrixK::fourierTransform() {
	isFourierTransformed = true;
}

/*
 * MatrixS
 */

MatrixS::MatrixS() : SymbolicTerm() {
	termID = 's';
}

MatrixS::MatrixS( MatrixS* s ) {
	indices[0] = s->indices[0];
	indices[1] = s->indices[1];
	termID = 's';
}

const string MatrixS::to_string() const {
	stringstream ss;
	ss << "S_" << indices;

	return ss.str();
}

MatrixS* MatrixS::copy() {
	MatrixS* cpy = new MatrixS();
	cpy->derivativeOrder = derivativeOrder;
	cpy->flavorLabel = flavorLabel;
	cpy->isInteracting = isInteracting;
	cpy->indices[0] = indices[0];
	cpy->indices[1] = indices[1];
	cpy->termID = 'B';

	return cpy;
}

/*
 * TermA
 */

TermA::TermA() : SymbolicTerm() {
	termID = 'A';
}

TermA::TermA( const TermA* A ) : SymbolicTerm() {
	termID = 'A';
}

const std::string TermA::to_string() const {
	return "A";
}

TermA* TermA::copy() {
	return new TermA();
}

bool TermA::operator==( const TermA &other ) const {
	return termID == other.termID;  // Really should always be true if other is valid.
}

/*
 * CoefficientFloat;
 */

CoefficientFloat::CoefficientFloat( double val ) : SymbolicTerm() {
	value = val;
	termID = 'L';
}

CoefficientFloat::CoefficientFloat( const CoefficientFloat* f ) : SymbolicTerm() {
	value = f-> value;
	termID = 'L';
}

const string CoefficientFloat::to_string() const {
	stringstream ss;
	ss << value;
	return ss.str();
}

CoefficientFloat* CoefficientFloat::copy() {
	CoefficientFloat* cpy = new CoefficientFloat( value );
	return cpy;
}

Sum CoefficientFloat::getDerivative() {
	return Sum();
}

double CoefficientFloat::eval() {
	return value;
}

/*
 * CoefficientFraction;
 */

CoefficientFraction::CoefficientFraction( int n, int d ) : SymbolicTerm() {
	num = n;
	den = d;
	termID = 'R';
}

CoefficientFraction::CoefficientFraction( CoefficientFraction* f ) : SymbolicTerm() {
	num = f->num;
	den = f->den;
	termID = 'R';
}

const string CoefficientFraction::to_string() const {
	stringstream ss;
	ss << num << " / " << den;
	return ss.str();
}

CoefficientFraction* CoefficientFraction::copy() {
	CoefficientFraction* cpy = new CoefficientFraction( num, den );
	return cpy;
}

Sum CoefficientFraction::getDerivative() {
	return Sum(); // TODO
}

double CoefficientFraction::eval() {
	return num / den;
}

/*
 * ***********************************************************************
 * TERM DERIVED CLASSES - EXPRESSION CLASSES
 * ***********************************************************************
 */

/*
 * Sum
 */

Sum::Sum() {
	terms = vector<SymbolicTerm*>();
	termID = 'S';
}

Sum::Sum( std::vector<SymbolicTerm*> thisTerms ) {
	terms = thisTerms;
	termID = 'S';
}

Sum::Sum( SymbolicTerm* term ) {
	terms = vector<SymbolicTerm*>();
	termID = 'S';
	terms.push_back( term );
}

Sum::Sum( const Sum* s ) {
	// TODO
}

Sum::~Sum() {
//	for ( int i = 0; i < terms.size(); i++ ) {
//		if ( terms[ i ] != 0 ) {
//			cout << "destroy " << terms[ i ]->to_string() << endl;
//			delete terms[ i ];
//		}
//	}
}

const std::string Sum::to_string() const {
	stringstream ss;
	bool SPLIT_SUMS_BY_LINE = false;
	unsigned int counter = 0;
	for ( vector<SymbolicTerm*>::const_iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		ss << (*iter)->to_string();
		if ( counter != (terms.size() - 1) ) {
			if ( SPLIT_SUMS_BY_LINE ) ss << "\n\n";
			ss << " + ";
		}
		counter++;
	}

	return ss.str();
}

Sum* Sum::copy() {
	Sum* cpy = new Sum();
	for ( vector<SymbolicTerm*>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		cpy->addTerm( (*iter)->copy() );
	}

	return cpy;
}

void Sum::simplify() {
	vector<SymbolicTerm*> simplifiedTerms = vector<SymbolicTerm*>();

	for ( vector<SymbolicTerm*>::iterator iter = terms.begin(); iter != terms.end(); ) {
		(*iter)->simplify();
		unpackTrivialExpression( (*iter) );
		string trimmedRepresentation = (*iter)->to_string();  // To be trimmed on next statement.
		boost::trim( trimmedRepresentation );
		if ( trimmedRepresentation == "0" or trimmedRepresentation == "{0}" or isZeroTrace( *iter ) ) {
			//delete *iter;
			iter = terms.erase( iter );
		} else {
			++iter;
		}
	}

	if ( terms.size() == 0 ) {
		CoefficientFloat* zero = new CoefficientFloat( 0.0 );
		terms.push_back( zero );
	}

}

void Sum::reduceTree() {
	for ( vector<SymbolicTerm*>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		if ( (*iter)->getTermID() == 'S' ) {
			(*iter)->reduceTree();
			Sum* iteratingSum = dynamic_cast<Sum*>(*iter); // TODO: Add exception check here.
			for ( vector<SymbolicTerm*>::iterator iter_term = iteratingSum->getIteratorBegin(); iter_term != iteratingSum->getIteratorEnd(); ++iter_term ) {
				unpackTrivialExpression( *iter_term );
			}

		} else if ( (*iter)->getTermID() == 'P' or (*iter)->getTermID() == 'T' ) {
			(*iter)->reduceTree();
			unpackTrivialExpression( *iter );
		}  else {
			unpackTrivialExpression( *iter );
		}
	}
}


Sum Sum::getDerivative() {
	reduceTree();
	Sum D = Sum();
	Sum termDerivative;
	for ( vector<SymbolicTerm*>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		termDerivative = (*iter)->getDerivative();
		cout << "in getDerivative: " << (*iter)->to_string() << endl;
		D.addTerm( termDerivative.copy() );
	}

	D.simplify();

	return D;
}

Sum* Sum::getExpandedExpr() {
	Sum* expandedSum = new Sum();
	for ( vector<SymbolicTerm*>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		if ( (*iter)->getTermID() == 'S' ) {
			Sum* s = dynamic_cast<Sum*>( *iter );
			expandedSum->addTerm( s->getExpandedExpr() );
		} else if ( (*iter)->getTermID() == 'P' ) {
			Product* prod = dynamic_cast<Product*>( *iter );
			expandedSum->addTerm( prod->getExpandedExpr() );
		}
	}

	return expandedSum;
}

void Sum::addTerm( SymbolicTerm* t ) {
	terms.push_back( t );
}

int Sum::getNumberOfTerms() {
	return terms.size();
}

void Sum::setAsNonInteracting() {
	for ( vector<SymbolicTerm*>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		(*iter)->setAsNonInteracting();
	}
}

bool Sum::operator==( const Sum &other ) const {
	return false; // TODO
}

vector<SymbolicTerm*>::iterator Sum::getIteratorBegin() {
	return terms.begin();
}

vector<SymbolicTerm*>::iterator Sum::getIteratorEnd() {
	return terms.end();
}

/*
 * Product
 */

Product::Product() : SymbolicTerm() {
	terms = vector<SymbolicTerm*>();
	termID = 'P';
}

Product::Product( vector<SymbolicTerm*> t ) : SymbolicTerm() {
	terms = t;
	termID = 'P';
}

Product::Product( SymbolicTerm* term ) : SymbolicTerm() {
	terms = vector<SymbolicTerm*>();
	terms.push_back( term );
	termID = 'P';

}

Product::Product( const Product* p ) {

}

Product::~Product() {
//	for ( int i = 0; i < terms.size(); i++ ) {
//		if ( terms[ i ] != 0 ) {
//			delete terms[ i ];
//		}
//	}
}

const string Product::to_string() const {
	stringstream ss;

	ss << " ";
	for ( vector<SymbolicTerm*>::const_iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		ss << "{" << (*iter)->to_string() << "}";
		if ( iter != terms.end()-- ) {
			ss << " ";
		}
	}

	return ss.str();
}

Product* Product::copy() {
	Product* cpy = new Product();
	for ( vector<SymbolicTerm*>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		cpy->addTerm( (*iter)->copy() );
	}

	return cpy;
}

void Product::simplify() {
	for ( vector<SymbolicTerm*>::iterator iter =  terms.begin(); iter != terms.end(); ++iter ) {
		(*iter)->simplify();
		unpackTrivialExpression( *iter );

		string trimmedRepresentation = (*iter)->to_string();  // To be trimmed on next statement.
		boost::trim( trimmedRepresentation );

		if ( trimmedRepresentation == "0" or trimmedRepresentation == "-0" or isZeroTrace( *iter ) ) {
			CoefficientFloat* zero = new CoefficientFloat( 0.0 );
			vector<SymbolicTerm*> zeroVector = vector<SymbolicTerm*>();
			zeroVector.push_back( zero );
			// delete terms?
			terms = zeroVector;
			return;
		} else if ( trimmedRepresentation == "1" ) {
			cout << "delete called in product simplify()" << endl;
			delete *iter;
			iter = terms.erase( iter );
		}
	}
}

void Product::reduceTree() {
	for ( vector<SymbolicTerm*>::iterator iter =  terms.begin(); iter != terms.end(); ++iter ) {
		if ( (*iter)->getTermID() == 'P' ) {
			(*iter)->reduceTree();
			Product* castProduct = dynamic_cast<Product*>( *iter );
			for ( vector<SymbolicTerm*>::iterator term_iter = castProduct->getIteratorBegin(); term_iter != castProduct->getIteratorEnd(); ++term_iter ) {
				unpackTrivialExpression( *term_iter );
			}
		} else if ( (*iter)->getTermID() == 'S' or (*iter)->getTermID() == 'T' ) {
			(*iter)->reduceTree();
			unpackTrivialExpression( *iter );
		} else {
			unpackTrivialExpression( *iter );
		}
	}
}

Sum Product::getDerivative() {
	reduceTree();
	Sum D = Sum();
	for ( unsigned int i = 0; i < terms.size(); i++ ) {
		Product* derivativeTerm = this->copy();
		derivativeTerm->terms[ i ] = dynamic_cast<SymbolicTerm*>( derivativeTerm->terms[ i ]->getDerivative().copy() );
		D.addTerm( derivativeTerm );
	}

	D.simplify();
	return D;
}
//def getExpandedExpr( self ):
//        checkGarbageCollection()
//        if len( self.terms ) == 0 or len( self.terms ) == 1:
//            # No expansion to be done; return the original object.
//            return self
//
//        elif len( self.terms ) == 2:
//            if isinstance( unpackTrivialExpr( self.terms[0] ), Sum ):  # Case 1: First term is a Sum, or both terms are a Sum.
//                expandedSum = Sum([])
//                factorA = self.terms[1]
//                if isinstance( factorA, Product ) and factorA.containsSum():
//                    factorA = factorA.getExpandedExpr()
//                for i in range( 0, len( self.terms[0].terms ) ):
//                    factorB = self.terms[0].terms[i]
//                    if isinstance( factorB, Product ) and factorB.containsSum():
//                        factorB = factorB.getExpandedExpr()
//                    expandedSum.terms.append( Product( [ deepcopy( factorB ), deepcopy( factorA ) ] ).getExpandedExpr() )
//                return expandedSum
//
//            elif isinstance( unpackTrivialExpr( self.terms[1] ), Sum ):  # Case 2: Second term is a Sum.
//                expandedSum = Sum([])
//                factorA = self.terms[0]
//                if isinstance( factorA, Product ) and factorA.containsSum():
//                    factorA = factorA.getExpandedExpr()
//                for i in range( 0, len( self.terms[1].terms ) ):
//                    factorB = self.terms[1].terms[i]
//                    if isinstance( factorB, Product ) and factorB.containsSum():
//                        factorB = factorB.getExpandedExpr()
//                    expandedSum.terms.append( Product( [ deepcopy( factorA ), deepcopy( factorB ) ] ).getExpandedExpr() )
//                return expandedSum
//            else:  # Case 3: Neither term is a sum, so no expansion is necessary.
//                return self
//
//        else:  # Recursively expand the product.
//            return Product( [ Product( self.terms[0:2] ).getExpandedExpr(), Product( self.terms[2:] ).getExpandedExpr() ] ).getExpandedExpr()
Sum* Product::getExpandedExpr() {
//	if ( terms.size() == 0 or terms.size() == 1 ) {
//		vector<SymbolicTerm*> trivialSum = vector<SymbolicTerm*>( this );
//		return Sum( trivialSum );
//	} else if ( terms.size() == 2 ) {
//		unpackTrivialExpression( terms[0] );
//		if ( terms[0]->getTermID() == 'S') {  // Case I: First term is a Sum, or both terms are a Sum.
//			Sum expandedSum = Sum();
//			Sum expandedFactorA;
//			SymbolicTerm* factorA = terms[1];
//			Sum* castSum = dynamic_cast<Sum*>( terms[0] );
//
//			if ( factorA->getTermID() == 'P' ) {
//				Product* castFactorA = dynamic_cast<Product*>( factorA );
//				if ( castFactorA->containsSum() ) {
//					expandedFactorA = castFactorA->getExpandedExpr();
//				}
//			}
//
//			for ( vector<SymbolicTerm*>::iterator iter = castSum->getIteratorBegin(); iter != castSum->getIteratorEnd(); ++iter ) {
//				Sum expandedFactorB;
//				SymbolicTerm* factorB = (*iter);
//
//				if ( factorB->getTermID() == 'P' ) {
//					Product* castFactorB = dynamic_cast<Product*>( factorB );
//					if ( castFactorB->containsSum() ) {
//						expandedFactorB = castFactorB->getExpandedExpr();
//					}
//				}
//
//				Product* finalFactorA, finalFactorB;
//				if ( expandedFactorA == 0 ) {  // If expandedFactorA remains null, then factorA did not require expansion.
//					vector<SymbolicTerm*> trivial = vector<SymbolicTerm*>();
//					trivial.push_back( factorA);
//					expandedFactorA = Sum( trivial ).copy();
//				}
//
//				if ( expandedFactorB == 0 ) { // If expandedFactor B remains null, then factorB did not require expansion.
//					vector<SymbolicTerm*> trivial = vector<SymbolicTerm*>();
//					trivial.push_back( factorB );
//					expandedFactorB = Sum( trivial ).copy();
//				}
//
//				// expandedSum.terms.append( Product( [ deepcopy( factorA ), deepcopy( factorB ) ] ).getExpandedExpr() )
//				vector<SymbolicTerm*> productTermsToExpand = vector<SymbolicTerm*>();
//				productTermsToExpand.push_back( expandedFactorA );
//				productTermsToExpand.push_back( expandedFactorB );
//				Product* productToExpand = new Product( productTermsToExpand );
//				vector<SymbolicTerm*> trivial = vector<SymbolicTerm*>();
//				trivial.push_back( )
//				expandedSum.addTerm(  )
//			}
//
//		}
//	}
	return new Sum();
}

void Product::addTerm( SymbolicTerm* t ) {
	terms.push_back( t );
}

int Product::getNumberOfTerms() {
	return terms.size();
}

void Product::setAsNonInteracting() {
	for ( vector<SymbolicTerm*>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		(*iter)->setAsNonInteracting();
	}
}

bool Product::operator==( const Sum &other ) const {
	return false;
}

bool Product::containsSum() {
	for ( vector<SymbolicTerm*>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		if ( (*iter)->getTermID() == 'S' ) return true;
	}
	return false;
}

vector<SymbolicTerm*>::const_iterator Product::getIteratorBegin() const {
	return terms.begin();
}

vector<SymbolicTerm*>::const_iterator Product::getIteratorEnd() const {
	return terms.end();
}

vector<SymbolicTerm*>::iterator Product::getIteratorBegin() {
	return terms.begin();
}

vector<SymbolicTerm*>::iterator Product::getIteratorEnd() {
	return terms.end();
}

/*
 * Trace
 */

Trace::Trace( SymbolicTerm* expr ) : SymbolicTerm() {

}

Trace::Trace( const Trace* ) : SymbolicTerm() {

}

Trace::~Trace() {

}

const string Trace::to_string() const {
	stringstream ss;
	ss << "Trace[";
	ss << "]";
	return ss.str();
}

void Trace::simplify() {

}

void Trace::reduceTree() {

}

Sum Trace::getDerivative(){
	return Sum(); // TODO
}

void Trace::setAsNonInteracting() {
	isInteracting = false;
}

bool Trace::operator==( const Trace &other ) const {
	return false; // TODO
}

void Trace::rewriteInKSFormalism() {

}

/*
 * FourierSum
 */

void unpackTrivialExpression( SymbolicTerm*& st ) {
	if ( st->getTermID() == 'P' ) {
		Product* castProduct = dynamic_cast<Product*>( st );
		if ( castProduct->terms.size() == 1 ) {
			st = castProduct->terms[ 0 ];
		}
	} else if ( st->getTermID() == 'S' ) {
		Sum* castSum = dynamic_cast<Sum*>( st );
		if ( castSum->terms.size() == 1) {
			st = castSum->terms[ 0 ];
		}
	}
}

bool isZeroTrace( SymbolicTerm* tr ) {
	return false; // TODO
}

std::ostream& operator<<( std::ostream& os, const TermA &obj ) {
	os << obj.to_string();
	return os;
}

std::ostream& operator<<( std::ostream& os, const MatrixM &obj ) {
	os << obj.to_string();
	return os;
}
