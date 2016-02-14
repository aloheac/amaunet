/*
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * High-order Perturbation Theory Analytics
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Symbolic Peturbation Theory Expression Objects Source Implementation
 *
 * v. 0.1		09 Feb 2016
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
	indices[ 0 ] = 0;
	indices[ 1 ] = 0;
	termID = '0';
}


SymbolicTerm::SymbolicTerm( const SymbolicTerm* st ) {

}

SymbolicTerm::~SymbolicTerm() { }

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

	// TODO: Shouldn't need dynamic_cast calls here -- polymorphism should take care of it.
	if ( termID == 'S' ) {
		copy = dynamic_cast<Sum*>( this )->copy();
	} else if ( termID == 'P' ) {
		copy = dynamic_cast<Product*>( this )->copy();
	} else if ( termID == 'K' ) {
		copy = dynamic_cast<MatrixK*>( this )->copy();
	} else if ( termID == 's' ) {
		copy = dynamic_cast<MatrixS*>( this )->copy();
	} else if ( termID == 'B' ) {
		copy = dynamic_cast<MatrixB*>( this )->copy();
	} else if ( termID == 'A' ) {
		copy = dynamic_cast<TermA*>( this )->copy();
	} else if ( termID == 'T' ) {
		copy = dynamic_cast<Trace*>( this )->copy();
	} else if ( termID == 'L' ) {
		copy = dynamic_cast<CoefficientFloat*>( this )->copy();
	} else if ( termID == 'R' ) {
		copy = dynamic_cast<CoefficientFraction*>( this )->copy();
	} else if ( termID == 'M' ) {
		copy = dynamic_cast<MatrixM*>( this )->copy();
	} else if ( termID == 'D' ) {
		copy = dynamic_cast<DetM*>( this )->copy();
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
				ss << "0";
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

		return Sum( new Product( derivativeTerms ) );
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
				ss << "0";
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
				ss << "0";
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
		ss << "D_" << flavorLabel << "_( " << indices[0] << ", " << indices[1] << " )";
	} else {
		ss << "K_" << flavorLabel << "_( " << indices[0] << ", " << indices[1] << " )";
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

const string MatrixS::to_string() const {
	stringstream ss;
	ss << "S_(" << indices[0] << ", " << indices[1] << ")";

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
 * DetM
 */

DetM::DetM( MatrixM thisMatrix ) : SymbolicTerm() {
	termID = 'D';
	isInverted = false;
	matrix = thisMatrix;
}

DetM::DetM( MatrixM thisMatrix, bool inverted ) : SymbolicTerm() {
	termID = 'D';
	isInverted = inverted;
	matrix = thisMatrix;
}

bool DetM::isDetInverted() {
	return isInverted;
}

void DetM::setAsNonInteracting() {
	isInteracting = false;
	matrix.setAsNonInteracting();
}

const std::string DetM::to_string() const {
	stringstream ss;
	ss << "Det[ " << matrix.to_string() << " ]";
	return ss.str();
}

bool DetM::operator==( const DetM &other ) const {
	return ( isInverted == other.isInverted and matrix == other.matrix );
}

Sum DetM::getDerivative() {
	//return Product( [ DetM( self.M, True, False ), Trace( Product( [ MatrixB( True, self.M.flavorLabel ), self.M.derivative() ] ) ) ] )
	// TODO: Add exception handling for an inverted or non-interacting determinant.
	Product* derivative = new Product();
	derivative->addTerm( new DetM( matrix ) );

	Product* traceExpr = new Product();
	traceExpr->addTerm( new MatrixB( matrix.getFlavorLabel() ) );
	traceExpr->addTerm( matrix.getDerivative().copy() );

	Trace* tr = new Trace( traceExpr );
	derivative->addTerm( tr );

	return Sum( derivative );
}

DetM* DetM::copy() {
	DetM* cpy = new DetM( matrix, isInverted );
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

Sum::Sum( const Sum &s ) {
	terms = vector<SymbolicTerm*>();
	termID = 'S';
	for ( vector<SymbolicTerm*>::const_iterator iter = s.terms.begin(); iter != s.terms.end(); ++iter ) {
			addTerm( (*iter)->copy() );
	}
}

Sum::~Sum() {
	for ( int i = 0; i < terms.size(); i++ ) {
			delete terms[ i ];
	}
}

Sum& Sum::operator=( const Sum &rhs ) {
	terms = vector<SymbolicTerm*>();
	termID = 'S';
	for ( vector<SymbolicTerm*>::const_iterator iter = rhs.terms.begin(); iter != rhs.terms.end(); ++iter ) {
			addTerm( (*iter)->copy() );
	}

	return *this;
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
	cpy->termID = 'S';
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
	vector<SymbolicTerm*> reducedExpression;
	for ( vector<SymbolicTerm*>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		if ( (*iter)->getTermID() == 'S' ) {
			(*iter)->reduceTree();
			Sum* iteratingSum = dynamic_cast<Sum*>(*iter); // TODO: Add exception check here.
			for ( vector<SymbolicTerm*>::iterator iter_term = iteratingSum->getIteratorBegin(); iter_term != iteratingSum->getIteratorEnd(); ++iter_term ) {
				unpackTrivialExpression( *iter_term );
				reducedExpression.push_back( *iter_term );
			}

		} else if ( (*iter)->getTermID() == 'P' or (*iter)->getTermID() == 'T' ) {
			(*iter)->reduceTree();
			unpackTrivialExpression( *iter );
			reducedExpression.push_back( *iter );
		}  else {
			unpackTrivialExpression( *iter );
			reducedExpression.push_back( *iter );
		}
	}

	terms = reducedExpression;
}


Sum Sum::getDerivative() {
	reduceTree();
	Sum D;
	SymbolicTerm* termDerivative;
	for ( vector<SymbolicTerm*>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		termDerivative = (*iter)->getDerivative().copy();
		D.addTerm( termDerivative );
	}

	D.simplify();

	return D;
}

Sum Sum::getExpandedExpr() {
	Sum expandedSum;
	for ( vector<SymbolicTerm*>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		if ( (*iter)->getTermID() == 'S' ) {
			Sum* s = dynamic_cast<Sum*>( *iter );
			expandedSum.addTerm( s->getExpandedExpr().copy() );
		} else if ( (*iter)->getTermID() == 'P' ) {
			Product* prod = dynamic_cast<Product*>( *iter );
			expandedSum.addTerm( prod->getExpandedExpr().copy() );
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

Product::~Product() {
	for ( int i = 0; i < terms.size(); i++ ) {
			delete terms[ i ];
	}
}

Product& Product::operator=( const Product &rhs ) {
	terms = vector<SymbolicTerm*>();
	termID = 'P';
	for ( vector<SymbolicTerm*>::const_iterator iter = rhs.terms.begin(); iter != rhs.terms.end(); ++iter ) {
		addTerm( (*iter)->copy() );
	}

	return *this;
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
	termID = 'P';
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
			delete *iter;
			iter = terms.erase( iter );
		}
	}
}

void Product::reduceTree() {
	vector<SymbolicTerm*> reducedExpression;
	for ( vector<SymbolicTerm*>::iterator iter =  terms.begin(); iter != terms.end(); ++iter ) {
		if ( (*iter)->getTermID() == 'P' ) {
			(*iter)->reduceTree();
			Product* castProduct = dynamic_cast<Product*>( *iter );
			for ( vector<SymbolicTerm*>::iterator term_iter = castProduct->getIteratorBegin(); term_iter != castProduct->getIteratorEnd(); ++term_iter ) {
				unpackTrivialExpression( *term_iter );
				reducedExpression.push_back( *term_iter );
			}
		} else if ( (*iter)->getTermID() == 'S' or (*iter)->getTermID() == 'T' ) {
			(*iter)->reduceTree();
			unpackTrivialExpression( *iter );
			reducedExpression.push_back( *iter );
		} else {
			unpackTrivialExpression( *iter );
			reducedExpression.push_back( *iter );
		}
	}

	terms = reducedExpression;
}

Sum Product::getDerivative() {
	reduceTree();
	Sum D;
	Product* derivativeTerm;
	for ( unsigned int i = 0; i < terms.size(); i++ ) {
		derivativeTerm = copy();
		derivativeTerm->terms[ i ] = dynamic_cast<SymbolicTerm*>( terms[ i ]->getDerivative().copy() );
		D.addTerm( derivativeTerm );
	}

	D.simplify();
	return D;
}

Sum Product::getExpandedExpr() {
	if ( terms.size() == 0 or terms.size() == 1 ) {

		// No expansion to be done for expressions with zero or one term.
		// Return a copy of this instance.
		return Sum( copy() );

	} else if ( terms.size() == 2 ) {

		SymbolicTerm* firstTerm = terms[0]->copy();
		SymbolicTerm* secondTerm = terms[1]->copy();

		unpackTrivialExpression( firstTerm );
		unpackTrivialExpression( secondTerm );

		if ( firstTerm->getTermID() == 'S' ) {  // Case I: First term is a Sum, or both terms are a Sum.

			Sum expandedSum;
			Sum* currentSum = dynamic_cast<Sum*>( firstTerm );

			SymbolicTerm* factorA = secondTerm;
			SymbolicTerm* factorB;
			Product *productA, *productB;

			if ( factorA->getTermID() == 'P' ) {
				productA = dynamic_cast<Product*>( factorA );
				if ( productA->containsSum() ) {
					factorA = productA->getExpandedExpr().copy();  // copy() returns Sum*.
				}
			}

			for ( vector<SymbolicTerm*>::iterator iter = currentSum->terms.begin(); iter != currentSum->terms.end(); ++iter ) {
				factorB = *iter;
				unpackTrivialExpression( factorB );  // Need?
				if ( factorB->getTermID() == 'P' ) {
					productB = dynamic_cast<Product*>( factorB );
					if ( productB->containsSum() ) {
						factorB = productB->getExpandedExpr().copy();  // copy() returns Sum*.
						unpackTrivialExpression( factorB );  // Need?
					}
				}

				vector<SymbolicTerm*> copiedFactors;
				copiedFactors.push_back( factorB->copy() );
				copiedFactors.push_back( factorA->copy() );

				SymbolicTerm* expandedProduct = Product( copiedFactors ).getExpandedExpr().copy();
				unpackTrivialExpression( expandedProduct );
				expandedSum.addTerm( expandedProduct );
			}

			return expandedSum;

		} else if ( secondTerm->getTermID() == 'S' ) {  // Case II: Second term is a sum.

			Sum expandedSum;
			Sum* currentSum = dynamic_cast<Sum*>( secondTerm );

			SymbolicTerm* factorA = firstTerm;
			SymbolicTerm* factorB;
			Product *productA, *productB;

			if ( factorA->getTermID() == 'P' ) {
				productA = dynamic_cast<Product*>( factorA );
				if ( productA->containsSum() ) {
					factorA = productA->getExpandedExpr().copy();  // copy() returns Sum*.
				}
			}

			for ( vector<SymbolicTerm*>::iterator iter = currentSum->terms.begin(); iter != currentSum->terms.end(); ++iter ) {
				factorB = *iter;
				unpackTrivialExpression( factorB );  // Need?
				if ( factorB->getTermID() == 'P' ) {
					productB = dynamic_cast<Product*>( factorB );
					if ( productB->containsSum() ) {
						factorB = productB->getExpandedExpr().copy();  // copy() returns Sum*.
						unpackTrivialExpression( factorB );  // Need?
					}
				}

				vector<SymbolicTerm*> copiedFactors;
				copiedFactors.push_back( factorA->copy() );
				copiedFactors.push_back( factorB->copy() );

				SymbolicTerm* expandedProduct = Product( copiedFactors ).getExpandedExpr().copy();
				unpackTrivialExpression( expandedProduct );
				expandedSum.addTerm( expandedProduct );
			}

			return expandedSum;

		} else {  // Case III: Neither term is a sum, so no expansion is necessary.
			return Sum( copy() );
		}

	} else {  // Recursively expand the product.

		Product headProduct, tailProduct, completeProduct;
		SymbolicTerm *expandedHeadProduct, *expandedTailProduct;

		headProduct.addTerm( terms[0]->copy() );  // Don't need copy here?
		headProduct.addTerm( terms[1]->copy() );

		for ( int i = 2; i < terms.size(); i++ ) {
			tailProduct.addTerm( terms[ i ]->copy() );
		}

		expandedHeadProduct = headProduct.getExpandedExpr().copy();
		expandedTailProduct = tailProduct.getExpandedExpr().copy();

		unpackTrivialExpression( expandedHeadProduct );  // Don't need these?
		unpackTrivialExpression( expandedTailProduct );

		completeProduct.addTerm( expandedHeadProduct );
		completeProduct.addTerm( expandedTailProduct );

		return completeProduct.getExpandedExpr();
	}
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

Trace::Trace( SymbolicTerm* thisExpr ) : SymbolicTerm() {
	expr = thisExpr;
	termID = 'T';
}

Trace::~Trace() {
	delete expr;
}

const string Trace::to_string() const {
	stringstream ss;
	ss << "Trace[ " << expr->to_string() << " ]";
	return ss.str();
}

Trace* Trace::copy() {
	return new Trace( expr->copy() );
}

void Trace::simplify() {
	expr->simplify();
}

void Trace::reduceTree() {
	expr->reduceTree();
}

Sum Trace::getDerivative(){
	return Sum( new Trace( expr->getDerivative().copy() ) );
}

void Trace::setAsNonInteracting() {
	isInteracting = false;
	expr->setAsNonInteracting();
}

bool Trace::operator==( const Trace &other ) const {
	return *expr == other;
}

void Trace::rewriteInKSFormalism() {

}

/*
 * Delta
 */

Delta::Delta( int a, int b ) {
	indices[0] = a;
	indices[1] = b;
	isBar = false;
}

Delta::Delta( int a, int b, bool type ) {
	indices[0] = a;
	indices[1] = b;
	isBar = type;
}

const std::string Delta::to_string() const {
	stringstream ss;
	ss << "Delta";
	if ( isBar ) ss << "Bar";
	ss << "( " << indices[0] << ", " << indices[1] << " )";
	return ss.str();
}

bool Delta::operator==( const Delta &other ) const {
	return isBar == other.isBar and indices[0] == other.indices[0] and indices[1] == other.indices[1];
}

bool Delta::isDeltaBar() {
	return isBar;
}

/*
 * FourierSum
 */

bool unpackTrivialExpression( SymbolicTerm*& st ) {  // TODO: Separate out into helper.
	if ( st->getTermID() == 'P' ) {
		SymbolicTerm* tmp;
		Product* castProduct = dynamic_cast<Product*>( st );
		if ( castProduct->terms.size() == 1 ) {
			tmp = castProduct->terms[ 0 ]->copy();
			delete st;
			st = tmp;

			if ( unpackTrivialExpression( st ) ) {
				unpackTrivialExpression( st );
			}

			return true;
		}
	} else if ( st->getTermID() == 'S' ) {
		SymbolicTerm* tmp;
		Sum* castSum = dynamic_cast<Sum*>( st );
		if ( castSum->terms.size() == 1) {
			tmp = castSum->terms[ 0 ]->copy();
			delete st;
			st = tmp;

			if ( unpackTrivialExpression( st ) ) {
				unpackTrivialExpression( st );
			}
			return true;
		}
	}

	return false;
}

bool isZeroTrace( SymbolicTerm* tr ) {
	if ( tr->getTermID() == 'T' ) {
		Trace* castTrace = dynamic_cast<Trace*>( tr );
		if ( castTrace->expr->getTermID() == 'S' ) {
			Sum* castExpr = dynamic_cast<Sum*>( castTrace->expr );
			if ( castExpr->getNumberOfTerms() == 0 ) return true;
		} else if ( castTrace->expr->getTermID() == 'P' ) {
			Product* castExpr = dynamic_cast<Product*>( castTrace->expr );
			if ( castExpr->getNumberOfTerms() == 0 ) return true;
		}

		string trimmedRepresentation = tr->to_string();  // To be trimmed on next statement.
		boost::trim( trimmedRepresentation );

		if ( trimmedRepresentation == "0" or trimmedRepresentation == "{0}" ) return true;
	}

	return false;
}

std::ostream& operator<<( std::ostream& os, const TermA &obj ) {
	os << obj.to_string();
	return os;
}

std::ostream& operator<<( std::ostream& os, const MatrixM &obj ) {
	os << obj.to_string();
	return os;
}
