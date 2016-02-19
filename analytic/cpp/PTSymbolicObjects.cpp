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
 * University of North Carolina at Chapel Hill
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
	return Sum( SymbolicTermPtr( new CoefficientFloat( 0.0 ) ) );
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

const char* SymbolicTerm::getFlavorLabel() {
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

SymbolicTermPtr SymbolicTerm::copy() {
	SymbolicTermPtr cpy( new SymbolicTerm() );
	cpy->derivativeOrder = derivativeOrder;
	cpy->isInteracting = isInteracting;
	cpy->flavorLabel = flavorLabel;
	cpy->indices[ 0 ] = indices[ 0 ];
	cpy->indices[ 1 ] = indices[ 1 ];
	termID = '0';

	return cpy;
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
	return Sum( SymbolicTermPtr( new GenericTestTerm( id, derivativeOrder + 1 ) ) );
}

SymbolicTermPtr GenericTestTerm::copy() {
	GenericTestTermPtr cpy( new GenericTestTerm( 0, 0 ) );
	cpy->id = id;
	cpy->derivativeOrder = derivativeOrder;
	cpy->termID = 'g';

	return SymbolicTermPtr( static_pointer_cast<SymbolicTerm>( cpy ) );
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

MatrixM::MatrixM( const char* thisFlavorLabel ) : SymbolicTerm() {
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
    MatrixMPtr D( new MatrixM( flavorLabel ) );
    D->derivativeOrder = derivativeOrder + 1;
    if ( !isInteracting ) {
    	D->setAsNonInteracting();
    }

    vector<SymbolicTermPtr> vecD;
    vecD.push_back( static_pointer_cast<SymbolicTerm>( D ) );

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

SymbolicTermPtr MatrixM::copy() {
	MatrixMPtr cpy( new MatrixM() );
	cpy->derivativeOrder = derivativeOrder;
	cpy->flavorLabel = flavorLabel;
	cpy->isInteracting = isInteracting;
	cpy->indices[0] = indices[0];
	cpy->indices[1] = indices[1];
	cpy->termID = 'M';

	return static_pointer_cast<SymbolicTerm>( cpy );
}

/*
 * MatrixB
 */

MatrixB::MatrixB() : SymbolicTerm() {
	termID = 'B';
}

MatrixB::MatrixB( const char* thisFlavorLabel )  : SymbolicTerm() {
	flavorLabel = thisFlavorLabel;
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
		MatrixMPtr dM( new MatrixM( flavorLabel ) );
		dM->derivativeOrder = derivativeOrder + 1;

		vector<SymbolicTermPtr> derivativeTerms;
		derivativeTerms.push_back( SymbolicTermPtr( new CoefficientFloat( -1.0 ) ) );
		derivativeTerms.push_back( SymbolicTermPtr( new MatrixB( flavorLabel ) ) );
		derivativeTerms.push_back( static_pointer_cast<SymbolicTerm>( dM ) );
		derivativeTerms.push_back( SymbolicTermPtr( new MatrixB( flavorLabel ) ) );

		return Sum( SymbolicTermPtr( new Product( derivativeTerms ) ) );
	} else {
		vector<SymbolicTermPtr> derivativeTerms;
		derivativeTerms.push_back( SymbolicTermPtr( new CoefficientFloat( 0.0 ) ) );
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

SymbolicTermPtr MatrixB::copy() {
	MatrixBPtr cpy( new MatrixB() );
	cpy->derivativeOrder = derivativeOrder;
	cpy->flavorLabel = flavorLabel;
	cpy->isInteracting = isInteracting;
	cpy->indices[0] = indices[0];
	cpy->indices[1] = indices[1];
	cpy->termID = 'B';

	return static_pointer_cast<SymbolicTerm>( cpy );
}

/*
 * MatrixK
 */

MatrixK::MatrixK() : SymbolicTerm() {
	isFourierTransformed = false;
	termID = 'K';
}

MatrixK::MatrixK( const char* thisFlavorLabel ) : SymbolicTerm() {
	isFourierTransformed = false;
	flavorLabel = thisFlavorLabel;
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

SymbolicTermPtr MatrixK::copy() {
	MatrixKPtr cpy( new MatrixK() );
	cpy->derivativeOrder = derivativeOrder;
	cpy->flavorLabel = flavorLabel;
	cpy->isInteracting = isInteracting;
	cpy->indices[0] = indices[0];
	cpy->indices[1] = indices[1];
	cpy->termID = 'K';

	return static_pointer_cast<SymbolicTerm>( cpy );
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

SymbolicTermPtr MatrixS::copy() {
	MatrixSPtr cpy( new MatrixS() );
	cpy->derivativeOrder = derivativeOrder;
	cpy->flavorLabel = flavorLabel;
	cpy->isInteracting = isInteracting;
	cpy->indices[0] = indices[0];
	cpy->indices[1] = indices[1];
	cpy->termID = 'B';

	return static_pointer_cast<SymbolicTerm>( cpy );
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
	// TODO: Add exception handling for an inverted or non-interacting determinant.
	ProductPtr derivative( new Product() );
	derivative->addTerm( SymbolicTermPtr( new DetM( matrix ) ) );

	ProductPtr traceExpr( new Product() );
	traceExpr->addTerm( SymbolicTermPtr( new MatrixB( matrix.getFlavorLabel() ) ) );
	traceExpr->addTerm( SymbolicTermPtr( matrix.getDerivative().copy() ) );

	SymbolicTermPtr tr = TracePtr( new Trace( static_pointer_cast<SymbolicTerm>( traceExpr ) ) );
	derivative->addTerm( tr );

	return Sum( static_pointer_cast<SymbolicTerm>( derivative ) );
}

SymbolicTermPtr DetM::copy() {
	SymbolicTermPtr cpy( new DetM( matrix, isInverted ) );
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

SymbolicTermPtr TermA::copy() {
	return SymbolicTermPtr( new TermA() );
}

bool TermA::operator==( const TermA &other ) const {
	return termID == other.termID;  // Really should always be true if other is a valid TermA.
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

SymbolicTermPtr CoefficientFloat::copy() {
	SymbolicTermPtr cpy( new CoefficientFloat( value ) );
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

SymbolicTermPtr CoefficientFraction::copy() {
	SymbolicTermPtr cpy( new CoefficientFraction( num, den ) );
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
	terms = vector<SymbolicTermPtr>();
	termID = 'S';
}

Sum::Sum( std::vector<SymbolicTermPtr> thisTerms ) {
	terms = thisTerms;
	termID = 'S';
}

Sum::Sum( SymbolicTermPtr term ) {
	terms = vector<SymbolicTermPtr>();
	termID = 'S';
	terms.push_back( term );
}

Sum::Sum( const Sum &s ) {
	terms = vector<SymbolicTermPtr>();
	termID = 'S';
	for ( vector<SymbolicTermPtr>::const_iterator iter = s.terms.begin(); iter != s.terms.end(); ++iter ) {
			addTerm( (*iter)->copy() );
	}
}

Sum::~Sum() {
	for ( int i = 0; i < terms.size(); i++ ) {
			terms[ i ].reset();
	}
}

Sum& Sum::operator=( const Sum &rhs ) {
	terms = vector<SymbolicTermPtr>();
	termID = 'S';
	for ( vector<SymbolicTermPtr>::const_iterator iter = rhs.terms.begin(); iter != rhs.terms.end(); ++iter ) {
			addTerm( (*iter)->copy() );
	}

	return *this;
}

const std::string Sum::to_string() const {
	stringstream ss;
	bool SPLIT_SUMS_BY_LINE = false;
	unsigned int counter = 0;
	for ( vector<SymbolicTermPtr>::const_iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		ss << (*iter)->to_string();
		if ( counter != (terms.size() - 1) ) {
			if ( SPLIT_SUMS_BY_LINE ) ss << "\n\n";
			ss << " + ";
		}
		counter++;
	}

	return ss.str();
}

SymbolicTermPtr Sum::copy() {
	SumPtr cpy( new Sum() );
	for ( vector<SymbolicTermPtr>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		cpy->addTerm( (*iter)->copy() );
	}
	cpy->termID = 'S';

	return static_pointer_cast<SymbolicTerm>( cpy );
}

void Sum::simplify() {
	vector<SymbolicTermPtr> simplifiedTerms;

	for ( vector<SymbolicTermPtr>::iterator iter = terms.begin(); iter != terms.end(); ) {
		(*iter)->simplify();
		unpackTrivialExpression( (*iter) );
		string trimmedRepresentation = (*iter)->to_string();  // To be trimmed on next statement.
		boost::trim( trimmedRepresentation );
		if ( trimmedRepresentation == "0" or trimmedRepresentation == "{0}" or isZeroTrace( *iter ) ) {
			(*iter).reset();  // Verify.
			iter = terms.erase( iter );
		} else {
			++iter;
		}
	}

	if ( terms.size() == 0 ) {
		SymbolicTermPtr zero( new CoefficientFloat( 0.0 ) );
		terms.push_back( zero );
	}

}

void Sum::reduceTree() {
	vector<SymbolicTermPtr> reducedExpression;
	for ( vector<SymbolicTermPtr>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		unpackTrivialExpression( *iter );
		if ( (*iter)->getTermID() == 'S' ) {
			(*iter)->reduceTree();
			SumPtr iteratingSum = dynamic_pointer_cast<Sum>(*iter); // TODO: Add exception check here.
			for ( vector<SymbolicTermPtr>::iterator iter_term = iteratingSum->getIteratorBegin(); iter_term != iteratingSum->getIteratorEnd(); ++iter_term ) {
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
	SymbolicTermPtr termDerivative( nullptr );
	for ( vector<SymbolicTermPtr>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		termDerivative = (*iter)->getDerivative().copy();
		D.addTerm( termDerivative );
	}

	D.simplify();

	return D;
}

Sum Sum::getExpandedExpr() {
	Sum expandedSum;
	for ( vector<SymbolicTermPtr>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		if ( (*iter)->getTermID() == 'S' ) {
			SumPtr s = dynamic_pointer_cast<Sum>( *iter );
			expandedSum.addTerm( s->getExpandedExpr().copy() );
		} else if ( (*iter)->getTermID() == 'P' ) {
			ProductPtr prod = dynamic_pointer_cast<Product>( *iter );
			expandedSum.addTerm( prod->getExpandedExpr().copy() );
		}
	}

	return expandedSum;
}

void Sum::addTerm( SymbolicTermPtr t ) {
	terms.push_back( t );
}

int Sum::getNumberOfTerms() {
	return (int)terms.size();
}

void Sum::setAsNonInteracting() {
	for ( vector<SymbolicTermPtr>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		(*iter)->setAsNonInteracting();
	}
}

bool Sum::operator==( const Sum &other ) const {
	return false; // TODO
}

vector<SymbolicTermPtr>::iterator Sum::getIteratorBegin() {
	return terms.begin();
}

vector<SymbolicTermPtr>::iterator Sum::getIteratorEnd() {
	return terms.end();
}

vector<SymbolicTermPtr>::const_iterator Sum::getIteratorBegin() const {
	return terms.begin();
}

vector<SymbolicTermPtr>::const_iterator Sum::getIteratorEnd() const {
	return terms.end();
}

/*
 * Product
 */

Product::Product() : SymbolicTerm() {
	terms = vector<SymbolicTermPtr>();
	termID = 'P';
}

Product::Product( vector<SymbolicTermPtr> t ) : SymbolicTerm() {
	terms = t;
	termID = 'P';
}

Product::Product( SymbolicTermPtr term ) : SymbolicTerm() {
	terms = vector<SymbolicTermPtr>();
	terms.push_back( term );
	termID = 'P';

}

Product::~Product() {
	for ( int i = 0; i < terms.size(); i++ ) {
			terms[ i ].reset();
	}
}

Product& Product::operator=( const Product &rhs ) {
	terms = vector<SymbolicTermPtr>();
	termID = 'P';
	for ( vector<SymbolicTermPtr>::const_iterator iter = rhs.terms.begin(); iter != rhs.terms.end(); ++iter ) {
		addTerm( (*iter)->copy() );
	}

	return *this;
}

const string Product::to_string() const {
	stringstream ss;

	ss << " ";
	for ( vector<SymbolicTermPtr>::const_iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		ss << "{" << (*iter)->to_string() << "}";
		if ( iter != terms.end()-- ) {
			ss << " ";
		}
	}

	return ss.str();
}

SymbolicTermPtr Product::copy() {
	ProductPtr cpy( new Product() );
	for ( vector<SymbolicTermPtr>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		cpy->addTerm( (*iter)->copy() );
	}
	termID = 'P';
	return static_pointer_cast<SymbolicTerm>( cpy );
}

void Product::simplify() {
	for ( vector<SymbolicTermPtr>::iterator iter =  terms.begin(); iter != terms.end(); ++iter ) {
		(*iter)->simplify();
		unpackTrivialExpression( *iter );

		string trimmedRepresentation = (*iter)->to_string();  // To be trimmed on next statement.
		boost::trim( trimmedRepresentation );

		if ( trimmedRepresentation == "0" or trimmedRepresentation == "-0" or isZeroTrace( *iter ) ) {
			CoefficientFloatPtr zero( new CoefficientFloat( 0.0 ) );
			vector<SymbolicTermPtr> zeroVector;
			zeroVector.push_back( zero );
			// delete terms?
			terms = zeroVector;
			return;
		} else if ( trimmedRepresentation == "1" ) {
			(*iter).reset(); // Verify.
			iter = terms.erase( iter );
		}
	}
}

void Product::reduceTree() {
	vector<SymbolicTermPtr> reducedExpression;
	for ( vector<SymbolicTermPtr>::iterator iter =  terms.begin(); iter != terms.end(); ++iter ) {
		unpackTrivialExpression( *iter );
		if ( (*iter)->getTermID() == 'P' ) {
			(*iter)->reduceTree();
			ProductPtr castProduct = dynamic_pointer_cast<Product>( *iter );
			for ( vector<SymbolicTermPtr>::iterator term_iter = castProduct->getIteratorBegin(); term_iter != castProduct->getIteratorEnd(); ++term_iter ) {
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
	ProductPtr derivativeTerm;
	for ( unsigned int i = 0; i < terms.size(); i++ ) {
		derivativeTerm = ProductPtr( static_pointer_cast<Product>( copy() ) );
		derivativeTerm->terms[ i ] = dynamic_pointer_cast<SymbolicTerm>( terms[ i ]->getDerivative().copy() );
		D.addTerm( static_pointer_cast<SymbolicTerm>( derivativeTerm ) );
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

		SymbolicTermPtr firstTerm = terms[0]->copy();
		SymbolicTermPtr secondTerm = terms[1]->copy();

		unpackTrivialExpression( firstTerm );
		unpackTrivialExpression( secondTerm );

		if ( firstTerm->getTermID() == 'S' ) {  // Case I: First term is a Sum, or both terms are a Sum.

			Sum expandedSum;
			SumPtr currentSum = dynamic_pointer_cast<Sum>( firstTerm );

			SymbolicTermPtr factorA = secondTerm;
			SymbolicTermPtr factorB;
			ProductPtr productA, productB;

			if ( factorA->getTermID() == 'P' ) {
				productA = dynamic_pointer_cast<Product>( factorA );
				if ( productA->containsSum() ) {
					factorA = productA->getExpandedExpr().copy();  // copy() returns Sum*.
				}
			}

			for ( vector<SymbolicTermPtr>::iterator iter = currentSum->terms.begin(); iter != currentSum->terms.end(); ++iter ) {
				factorB = *iter;
				unpackTrivialExpression( factorB );  // Need?
				if ( factorB->getTermID() == 'P' ) {
					productB = dynamic_pointer_cast<Product>( factorB );
					if ( productB->containsSum() ) {
						factorB = productB->getExpandedExpr().copy();  // copy() returns Sum*.
						unpackTrivialExpression( factorB );  // Need?
					}
				}

				vector<SymbolicTermPtr> copiedFactors;
				copiedFactors.push_back( factorB->copy() );
				copiedFactors.push_back( factorA->copy() );

				SymbolicTermPtr expandedProduct = Product( copiedFactors ).getExpandedExpr().copy();
				unpackTrivialExpression( expandedProduct );
				expandedSum.addTerm( expandedProduct );
			}

			return expandedSum;

		} else if ( secondTerm->getTermID() == 'S' ) {  // Case II: Second term is a sum.

			Sum expandedSum;
			SumPtr currentSum = dynamic_pointer_cast<Sum>( secondTerm );

			SymbolicTermPtr factorA = firstTerm;
			SymbolicTermPtr factorB;
			ProductPtr productA, productB;

			if ( factorA->getTermID() == 'P' ) {
				productA = dynamic_pointer_cast<Product>( factorA );
				if ( productA->containsSum() ) {
					factorA = productA->getExpandedExpr().copy();  // copy() returns Sum*.
				}
			}

			for ( vector<SymbolicTermPtr>::iterator iter = currentSum->terms.begin(); iter != currentSum->terms.end(); ++iter ) {
				factorB = *iter;
				unpackTrivialExpression( factorB );  // Need?
				if ( factorB->getTermID() == 'P' ) {
					productB = dynamic_pointer_cast<Product>( factorB );
					if ( productB->containsSum() ) {
						factorB = productB->getExpandedExpr().copy();  // copy() returns Sum*.
						unpackTrivialExpression( factorB );  // Need?
					}
				}

				vector<SymbolicTermPtr> copiedFactors;
				copiedFactors.push_back( factorA->copy() );
				copiedFactors.push_back( factorB->copy() );

				SymbolicTermPtr expandedProduct = Product( copiedFactors ).getExpandedExpr().copy();
				unpackTrivialExpression( expandedProduct );
				expandedSum.addTerm( expandedProduct );
			}

			return expandedSum;

		} else {  // Case III: Neither term is a sum, so no expansion is necessary.
			return Sum( copy() );
		}

	} else {  // Recursively expand the product.

		Product headProduct, tailProduct, completeProduct;
		SymbolicTermPtr expandedHeadProduct, expandedTailProduct;

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

void Product::addTerm( SymbolicTermPtr t ) {
	terms.push_back( t );
}

int Product::getNumberOfTerms() {
	return (int)terms.size();
}

void Product::setAsNonInteracting() {
	for ( vector<SymbolicTermPtr>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		(*iter)->setAsNonInteracting();
	}
}

bool Product::operator==( const Sum &other ) const {
	return false;
}

bool Product::containsSum() {
	for ( vector<SymbolicTermPtr>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		if ( (*iter)->getTermID() == 'S' ) return true;
	}
	return false;
}

vector<SymbolicTermPtr>::const_iterator Product::getIteratorBegin() const {
	return terms.begin();
}

vector<SymbolicTermPtr>::const_iterator Product::getIteratorEnd() const {
	return terms.end();
}

vector<SymbolicTermPtr>::iterator Product::getIteratorBegin() {
	return terms.begin();
}

vector<SymbolicTermPtr>::iterator Product::getIteratorEnd() {
	return terms.end();
}

/*
 * Trace
 */

Trace::Trace( SymbolicTermPtr thisExpr ) : SymbolicTerm() {
	expr = thisExpr;
	termID = 'T';
}

Trace::~Trace() {
	expr.reset();  // Verify.
}

const string Trace::to_string() const {
	stringstream ss;
	ss << "Trace[ " << expr->to_string() << " ]";
	return ss.str();
}

SymbolicTermPtr Trace::copy() {
	return SymbolicTermPtr( new Trace( expr->copy() ) );
}

void Trace::simplify() {
	expr->simplify();
}

void Trace::reduceTree() {
	expr->reduceTree();
}

Sum Trace::getDerivative(){
	return Sum( SymbolicTermPtr( new Trace( expr->getDerivative().copy() ) ) );
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

bool unpackTrivialExpression( SymbolicTermPtr& st ) {  // TODO: Separate out into helper.
	if ( st->getTermID() == 'P' ) {
		SymbolicTermPtr tmp;
		ProductPtr castProduct = dynamic_pointer_cast<Product>( st );
		if ( castProduct->terms.size() == 1 ) {
			tmp = castProduct->terms[ 0 ]->copy();
			st.reset();  // Verify.
			st = tmp;

			if ( unpackTrivialExpression( st ) ) {
				unpackTrivialExpression( st );
			}

			return true;
		}
	} else if ( st->getTermID() == 'S' ) {
		SymbolicTermPtr tmp;
		SumPtr castSum = dynamic_pointer_cast<Sum>( st );
		if ( castSum->terms.size() == 1) {
			tmp = castSum->terms[ 0 ]->copy();
			st.reset();  // Verify.
			st = tmp;

			if ( unpackTrivialExpression( st ) ) {
				unpackTrivialExpression( st );
			}

			return true;
		}
	}

	return false;
}

bool isZeroTrace( SymbolicTermPtr tr ) {
	if ( tr->getTermID() == 'T' ) {
		TracePtr castTrace = dynamic_pointer_cast<Trace>( tr );
		if ( castTrace->expr->getTermID() == 'S' ) {
			SumPtr castExpr = dynamic_pointer_cast<Sum>( castTrace->expr );
			if ( castExpr->getNumberOfTerms() == 0 ) return true;
		} else if ( castTrace->expr->getTermID() == 'P' ) {
			ProductPtr castExpr = dynamic_pointer_cast<Product>( castTrace->expr );
			if ( castExpr->getNumberOfTerms() == 0 ) return true;
		}

		string trimmedRepresentation = tr->to_string();  // To be trimmed on next statement.
		boost::trim( trimmedRepresentation );

		if ( trimmedRepresentation == "0" or trimmedRepresentation == "{0}" ) return true;
	}

	return false;
}

Sum distributeTrace( SymbolicTermPtr tr ) {
	if ( tr->getTermID() != 'T' ) {
		cout << "***ERROR: An object other than a trace was passed to distributeTrace()." << endl;
		return Sum();  // TODO: Throw exception.
	}

	// The passed SymbolicTermPtr is a known instance of Trace; cast it.
	TracePtr castTrace = static_pointer_cast<Trace>( tr );

	// Simple case: a single product inside a trace has been passed in.
	if ( castTrace->expr->getTermID() == 'P' ) {
		// Consider relatively more complicated cases, where the argument to the trace is an expression like
		// a ( b c + d e ); we must expand the product and reduce the representation.

		// Argument to the trace is a known instance of a Product; cast it.
		ProductPtr castExpr = static_pointer_cast<Product>( castTrace->expr );

		SymbolicTermPtr expandedArgument( castExpr->copy() );  // Set scope; expression will be expanded next if necessary.
		if ( castExpr->containsSum() ) {
			expandedArgument = castExpr->getExpandedExpr().copy();
			expandedArgument->reduceTree();
		}

		// If the expansion resulted in the product turning into a sum, evaluate it appropriately through a recursive
		// call.
		if ( expandedArgument->getTermID() == 'S' ) {
			return distributeTrace( Trace( expandedArgument ).copy() );

		} else if ( expandedArgument->getTermID() == 'P' ){  // Otherwise, if the expression is a simple, single product, proceed.
			Product simplifiedProduct;
			vector<SymbolicTermPtr> matrixTerms;
			ProductPtr castExpandedArgument = static_pointer_cast<Product>( expandedArgument );
			for ( vector<SymbolicTermPtr>::iterator iter = castExpandedArgument->getIteratorBegin(); iter != castExpandedArgument->getIteratorEnd(); ++iter ) {
				if ( (*iter)->getTermID() == 'L' or (*iter)->getTermID() == 'R') {
					simplifiedProduct.addTerm( (*iter)->copy() );
				} else if ( (*iter)->getTermID() == 'B' or (*iter)->getTermID() == 'M' ) {
					matrixTerms.push_back( (*iter)->copy() );
				}
			}

			SymbolicTermPtr productOfMatrixTerms = static_pointer_cast<SymbolicTerm>( Product( matrixTerms ).copy() );
			unpackTrivialExpression( productOfMatrixTerms );
			simplifiedProduct.addTerm( Trace( productOfMatrixTerms ).copy() );
			return Sum( simplifiedProduct.copy() );

		} else {  // Expanded argument is somehow not a Sum or Product; invalid expression!
			cout << "***ERROR: Something bad happened in distributeTraces()." << endl;
			return Sum();
		}

	// Typical case: a sum of terms inside a trace has been passed in to be evaluated.
	} else if ( castTrace->expr->getTermID() == 'S' ){
		Sum distributedSum;
		SumPtr castExpr = static_pointer_cast<Sum>( castTrace->expr );

		for ( vector<SymbolicTermPtr>::iterator iter = castExpr->getIteratorBegin(); iter != castExpr->getIteratorEnd(); ++iter ) {
			distributedSum.addTerm( distributeTrace( Trace( (*iter)->copy() ).copy() ).copy() );  // TODO: Expensive!
		}
		distributedSum.reduceTree();
		return distributedSum;
	} else {  // Something trivial has been passed in; just return it.
		return Sum( tr );  // May need to call unpackTrivialExpression on the returned object.
	}
}

Sum distributeAllTraces( SymbolicTermPtr expr ) {
	// First, before doing anything, check if the expression can be distributed as a while. However, only do so if
	// necessary since an expansion is computationally intensive.
	if ( expr->getTermID() == 'P' ) {
		ProductPtr castProduct = static_pointer_cast<Product>( expr );
		if ( castProduct->containsSum() ) {
			expr = static_pointer_cast<SymbolicTerm>( castProduct->getExpandedExpr().copy() );
		}
	}

	// Reduce the expression tree.
	expr->reduceTree();

	SumPtr distributedExpr( new Sum() );
	if ( expr->getTermID() == 'S' ) {
		SumPtr castSum = static_pointer_cast<Sum>( expr );
		for ( vector<SymbolicTermPtr>::iterator iter = castSum->getIteratorBegin(); iter != castSum->getIteratorEnd(); ++iter ) {  // Loop over terms.
			distributedExpr->addTerm( distributeAllTraces( *iter ).copy() );
		}
	} else if ( expr->getTermID() == 'P' ) {
		ProductPtr castProduct = static_pointer_cast<Product>( expr );
		Product distributedProduct;
		for ( vector<SymbolicTermPtr>::iterator iter = castProduct->getIteratorBegin(); iter != castProduct->getIteratorEnd(); ++iter  ) {
			if ( (*iter)->getTermID() == 'T' ) {
				distributedProduct.addTerm( distributeTrace( *iter ).copy() );
			} else {
				distributedProduct.addTerm( *iter );
			}
		}

		SymbolicTermPtr expandedAndDistributedProduct;
		if ( distributedProduct.containsSum() ) {
			expandedAndDistributedProduct = SymbolicTermPtr( distributedProduct.getExpandedExpr().copy() );
		} else {
			expandedAndDistributedProduct = SymbolicTermPtr( distributedProduct.copy() );
		}

		distributedExpr->addTerm( expandedAndDistributedProduct );
	} else if ( expr->getTermID() == 'T' ) {
		distributedExpr->addTerm( distributeTrace( expr ).copy() ); // need expr.copy() as argument?
	}

	distributedExpr->reduceTree();
	return *distributedExpr;
}

/*
 * INPUT REDIRECTION OPERATOR OVERLOADS
 */

std::ostream& operator<<( std::ostream& os, const TermA &obj ) {
	os << obj.to_string();
	return os;
}

std::ostream& operator<<( std::ostream& os, const MatrixM &obj ) {
	os << obj.to_string();
	return os;
}
