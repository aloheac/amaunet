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

#include <sstream>
#include <cstring>
#include <boost/algorithm/string/trim.hpp>
#include "PTSymbolicObjects.h"
#include "PathIntegration.h"
#include <iostream>
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
	termID = TermTypes::INVALID_TERM ;
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

void SymbolicTerm::setIndices( int i, int j ) {
	indices[0] = i;
	indices[1] = j;
}

TermTypes SymbolicTerm::getTermID() {
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
	termID = TermTypes::INVALID_TERM ;

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
	termID = TermTypes::GENERIC_TEST_TERM ;
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
	cpy->termID = TermTypes::GENERIC_TEST_TERM;

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
	termID = TermTypes::MATRIX_M;
};

MatrixM::MatrixM( const char* thisFlavorLabel ) : SymbolicTerm() {
	termID = TermTypes::MATRIX_M;
	flavorLabel = thisFlavorLabel;
}

bool MatrixM::operator==( const MatrixM &other ) const {
	return ( derivativeOrder == other.derivativeOrder ) and
				( isInteracting == other.isInteracting ) and
				( flavorLabel == other.flavorLabel ) and
				other.termID == TermTypes::MATRIX_M;
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
	cpy->termID = TermTypes::MATRIX_M;;

	return static_pointer_cast<SymbolicTerm>( cpy );
}

/*
 * MatrixB
 */

MatrixB::MatrixB() : SymbolicTerm() {
	termID = TermTypes::MATRIX_B;
}

MatrixB::MatrixB( const char* thisFlavorLabel )  : SymbolicTerm() {
	flavorLabel = thisFlavorLabel;
	termID = TermTypes::MATRIX_B;
}

bool MatrixB::operator==( const MatrixB &other ) const {
	return ( derivativeOrder == other.derivativeOrder ) and
			( isInteracting == other.isInteracting ) and
			( flavorLabel == other.flavorLabel ) and
			( other.termID == TermTypes::MATRIX_B );
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
	cpy->termID = TermTypes::MATRIX_B;

	return static_pointer_cast<SymbolicTerm>( cpy );
}

/*
 * MatrixK
 */

MatrixK::MatrixK() : SymbolicTerm() {
	isFourierTransformed = false;
	termID = TermTypes::MATRIX_K;
}

MatrixK::MatrixK( const char* thisFlavorLabel ) : SymbolicTerm() {
	isFourierTransformed = false;
	flavorLabel = thisFlavorLabel;
	termID = TermTypes::MATRIX_K;
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
	cpy->termID = TermTypes::MATRIX_K;

	return static_pointer_cast<SymbolicTerm>( cpy );
}

void MatrixK::fourierTransform() {
	isFourierTransformed = true;
}

/*
 * MatrixS
 */

MatrixS::MatrixS() : SymbolicTerm() {
	termID = TermTypes::MATRIX_S;
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
	cpy->termID = TermTypes::MATRIX_S;

	return static_pointer_cast<SymbolicTerm>( cpy );
}

/*
 * DetM
 */

DetM::DetM( MatrixM thisMatrix ) : SymbolicTerm() {
	termID = TermTypes::DET_M;
	isInverted = false;
	matrix = thisMatrix;
}

DetM::DetM( MatrixM thisMatrix, bool inverted ) : SymbolicTerm() {
	termID = TermTypes::DET_M;
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
	termID = TermTypes::TERM_A;
}

TermA::TermA( const TermA* A ) : SymbolicTerm() {
	termID = TermTypes::TERM_A;
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
	termID = TermTypes::COEFFICIENT_FLOAT;
}

CoefficientFloat::CoefficientFloat( const CoefficientFloat* f ) : SymbolicTerm() {
	value = f-> value;
	termID = TermTypes::COEFFICIENT_FLOAT;
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
	return Sum( SymbolicTermPtr( new CoefficientFloat( 0.0 ) ) );
}

double CoefficientFloat::eval() {
	return value;
}

/*
 * CoefficientFraction;
 */

CoefficientFraction::CoefficientFraction() : SymbolicTerm() {
	num = 0;
	den = 1;
	termID = TermTypes::COEFFICIENT_FRACTION;
}

CoefficientFraction::CoefficientFraction( int n, int d ) : SymbolicTerm() {
	num = n;
	den = d;
	termID = TermTypes::COEFFICIENT_FRACTION;
}

CoefficientFraction::CoefficientFraction( CoefficientFraction* f ) : SymbolicTerm() {
	num = f->num;
	den = f->den;
	termID = TermTypes::COEFFICIENT_FRACTION;
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
	return Sum( SymbolicTermPtr( new CoefficientFloat( 0.0 ) ) );
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
	termID = TermTypes::SUM;
}

Sum::Sum( std::vector<SymbolicTermPtr> thisTerms ) {
	terms = thisTerms;
	termID = TermTypes::SUM;
}

Sum::Sum( SymbolicTermPtr term ) {
	terms = vector<SymbolicTermPtr>();
	termID = TermTypes::SUM;
	terms.push_back( term );
}

Sum::Sum( const Sum &s ) {
	terms = vector<SymbolicTermPtr>();
	termID = TermTypes::SUM;
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
	termID = TermTypes::SUM;
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
	cpy->termID = TermTypes::SUM;

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
			iter = terms.erase( iter );  // Note that we shouldn't need to check for a break condition as we do for
		} else {                         // Product::simplify() since the for loop doesn't advance the iterator; we
			++iter;                      // do it ourselves.
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
		if ( (*iter)->getTermID() == TermTypes::SUM ) {
			(*iter)->reduceTree();
			SumPtr iteratingSum = dynamic_pointer_cast<Sum>(*iter); // TODO: Add exception check here.
			for ( vector<SymbolicTermPtr>::iterator iter_term = iteratingSum->getIteratorBegin(); iter_term != iteratingSum->getIteratorEnd(); ++iter_term ) {
				unpackTrivialExpression( *iter_term );
				reducedExpression.push_back( *iter_term );
			}

		} else if ( (*iter)->getTermID() == TermTypes::PRODUCT or (*iter)->getTermID() == TermTypes::TRACE ) {
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
		if ( (*iter)->getTermID() == TermTypes::SUM ) {
			SumPtr s = dynamic_pointer_cast<Sum>( *iter );
			expandedSum.addTerm( s->getExpandedExpr().copy() );
		} else if ( (*iter)->getTermID() == TermTypes::PRODUCT ) {
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
	termID = TermTypes::PRODUCT;
}

Product::Product( vector<SymbolicTermPtr> t ) : SymbolicTerm() {
	terms = t;
	termID = TermTypes::PRODUCT;
}

Product::Product( SymbolicTermPtr term ) : SymbolicTerm() {
	terms = vector<SymbolicTermPtr>();
	terms.push_back( term );
	termID = TermTypes::PRODUCT;

}

Product::~Product() {
	for ( int i = 0; i < terms.size(); i++ ) {
			terms[ i ].reset();
	}
}

Product& Product::operator=( const Product &rhs ) {
	terms = vector<SymbolicTermPtr>();
	termID = TermTypes::PRODUCT;
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
	termID = TermTypes::PRODUCT;
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

			// Note that we could have the case where the last element of the vector could be removed, and the call
			// vector::erase() returns the iterator terms.end(). The for loop then advances the iterator, which
			// becomes invalid and results in a seg fault on the next recursive call. Therefore, check if the
			// iterator has reached the end of the data structure, and if so, break out of the loop.
			if ( iter == terms.end() ) break;
		}
	}
}

void Product::reduceTree() {
	vector<SymbolicTermPtr> reducedExpression;
	for ( vector<SymbolicTermPtr>::iterator iter =  terms.begin(); iter != terms.end(); ++iter ) {
		unpackTrivialExpression( *iter );
		if ( (*iter)->getTermID() == TermTypes::PRODUCT ) {
			(*iter)->reduceTree();
			ProductPtr castProduct = dynamic_pointer_cast<Product>( *iter );
			for ( vector<SymbolicTermPtr>::iterator term_iter = castProduct->getIteratorBegin(); term_iter != castProduct->getIteratorEnd(); ++term_iter ) {
				unpackTrivialExpression( *term_iter );
				reducedExpression.push_back( *term_iter );
			}
		} else if ( (*iter)->getTermID() == TermTypes::SUM or (*iter)->getTermID() == TermTypes::TRACE ) {
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

		if ( firstTerm->getTermID() == TermTypes::SUM ) {  // Case I: First term is a Sum, or both terms are a Sum.

			Sum expandedSum;
			SumPtr currentSum = dynamic_pointer_cast<Sum>( firstTerm );

			SymbolicTermPtr factorA = secondTerm;
			SymbolicTermPtr factorB;
			ProductPtr productA, productB;

			if ( factorA->getTermID() == TermTypes::PRODUCT ) {
				productA = dynamic_pointer_cast<Product>( factorA );
				if ( productA->containsSum() ) {
					factorA = productA->getExpandedExpr().copy();  // copy() returns Sum*.
				}
			}

			for ( vector<SymbolicTermPtr>::iterator iter = currentSum->terms.begin(); iter != currentSum->terms.end(); ++iter ) {
				factorB = *iter;
				unpackTrivialExpression( factorB );  // Need?
				if ( factorB->getTermID() == TermTypes::PRODUCT ) {
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

		} else if ( secondTerm->getTermID() == TermTypes::SUM ) {  // Case II: Second term is a sum.

			Sum expandedSum;
			SumPtr currentSum = dynamic_pointer_cast<Sum>( secondTerm );

			SymbolicTermPtr factorA = firstTerm;
			SymbolicTermPtr factorB;
			ProductPtr productA, productB;

			if ( factorA->getTermID() == TermTypes::PRODUCT ) {
				productA = dynamic_pointer_cast<Product>( factorA );
				if ( productA->containsSum() ) {
					factorA = productA->getExpandedExpr().copy();  // copy() returns Sum*.
				}
			}

			for ( vector<SymbolicTermPtr>::iterator iter = currentSum->terms.begin(); iter != currentSum->terms.end(); ++iter ) {
				factorB = *iter;
				unpackTrivialExpression( factorB );  // Need?
				if ( factorB->getTermID() == TermTypes::PRODUCT ) {
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
		if ( (*iter)->getTermID() == TermTypes::SUM ) return true;
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
	termID = TermTypes::TRACE;
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

void Trace::rewriteInKSFormalism() {  // TODO: Verify first derivative is taken for MatrixM exactly.
	if ( expr->getTermID() != TermTypes::PRODUCT ) {
		cout << "***ERROR: Argument of a Trace must be a single Product before rewriting in KS formalism. Distribute trace operators first." << endl;
	}

	ProductPtr castExpr = static_pointer_cast<Product>( expr );
	ProductPtr newExpr( new Product() );
	for ( vector<SymbolicTermPtr>::iterator iter = castExpr->getIteratorBegin(); iter != castExpr->getIteratorEnd(); ++iter ) {
		if ( (*iter)->getTermID() == TermTypes::MATRIX_B ) {
			newExpr->addTerm( SymbolicTermPtr( new MatrixK( (*iter)->getFlavorLabel() ) ) );
		} else if( (*iter)->getTermID() == TermTypes::MATRIX_M ) {
			newExpr->addTerm( SymbolicTermPtr( new MatrixS() ) );
		} else {
			cout << "***ERROR: Unexpected term encountered in Trace argument when rewriting in KS formalism." << endl;
		}
	}

	expr.reset();  // Delete current expression.
	expr = newExpr;
}

/*
 * Delta
 */

Delta::Delta( int a, int b ) {
	indices[0] = a;
	indices[1] = b;
	isBar = false;
	termID = TermTypes::DELTA;
}

Delta::Delta( int a, int b, bool type ) {
	indices[0] = a;
	indices[1] = b;
	isBar = type;
	termID = TermTypes::DELTA;
}

const std::string Delta::to_string() const {
	stringstream ss;
	ss << "Delta";
	if ( isBar ) ss << "Bar";
	ss << "( " << indices[0] << ", " << indices[1] << " )";
	return ss.str();
}

SymbolicTermPtr Delta::copy() {
	return SymbolicTermPtr( new Delta( indices[0], indices[1] ) );
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

FourierSum::FourierSum( vector<IndexContraction> i, int orderInK ) {
	indices = i;
	order = orderInK;
	termID = TermTypes::FOURIER_SUM;
}

FourierSum::~FourierSum() {
	indices.clear();
}

const string FourierSum::to_string() const {
	stringstream ss;
	ss << "FourierSum[";
	for ( vector<IndexContraction>::const_iterator indexPair = indices.begin(); indexPair != indices.end(); ++indexPair ) {
		ss << " ( " << indexPair->i << ", " << indexPair->j << " ) ";
	}
	ss << "]";

	return ss.str();
}

bool FourierSum::operator==( const FourierSum &other ) const {
	vector<IndexContraction> lhs( indices );
	vector<IndexContraction> rhs( other.indices );

	if ( lhs.size() != rhs.size() ) return false;

	sort( lhs.begin(), lhs.end() );
	sort( rhs.begin(), rhs.end() );

	for ( int i = 0; i < lhs.size(); i++ ) {
		if ( lhs[ i ].i != rhs[ i ].i or lhs[ i ].j != rhs[ i ].j ) return false;
	}

	return true;
}

/*
 * ***********************************************************************
 * GENERIC HELPER FUNCTIONS
 * ***********************************************************************
 */

bool unpackTrivialExpression( SymbolicTermPtr& st ) {  // TODO: Separate out into helper.
	if ( st->getTermID() == TermTypes::PRODUCT ) {
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
	} else if ( st->getTermID() == TermTypes::SUM ) {
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
	if ( tr->getTermID() == TermTypes::TRACE ) {
		TracePtr castTrace = dynamic_pointer_cast<Trace>( tr );
		if ( castTrace->expr->getTermID() == TermTypes::SUM ) {
			SumPtr castExpr = dynamic_pointer_cast<Sum>( castTrace->expr );
			if ( castExpr->getNumberOfTerms() == 0 ) return true;
		} else if ( castTrace->expr->getTermID() == TermTypes::PRODUCT ) {
			ProductPtr castExpr = dynamic_pointer_cast<Product>( castTrace->expr );
			if ( castExpr->getNumberOfTerms() == 0 ) return true;
		}

		string trimmedRepresentation = tr->to_string();  // To be trimmed on next statement.
		boost::trim( trimmedRepresentation );

		if ( trimmedRepresentation == "0" or trimmedRepresentation == "{0}" ) return true;
	}

	return false;
}

int getTerminatedContraction( map<int, int> contractedIndexMapping, int index ) {
	if ( contractedIndexMapping.count( index ) == 0 or contractedIndexMapping[ index ] == index  ) {
		return index;
	} else {
		return getTerminatedContraction( contractedIndexMapping, contractedIndexMapping[ index ] );
	}
}

map<int, int> constructContractionDictionary( DeltaContractionSet contractions ) {
	contractions.orderContractionIndices();
	contractions.sortContractions();

	map<int, int> contractedIndexMapping;

	for ( vector<IndexContraction>::iterator indexPair = contractions.getIteratorBegin(); indexPair != contractions.getIteratorEnd(); ++indexPair ) {
		// Recall that always j >= i.
		if ( contractedIndexMapping.count( indexPair->j ) == 0 ) {  // j is not in the dictionary, add the mapping.
			contractedIndexMapping[ indexPair->j ] = indexPair->i;
		} else if ( contractedIndexMapping[ indexPair->j ] > indexPair->i ) {
			contractedIndexMapping[ indexPair->j ] = indexPair->i;
		}

		if ( contractedIndexMapping.count( indexPair->i ) == 0 ) {  // i has not been referenced before, add the mapping
			contractedIndexMapping[ indexPair->i ] = indexPair->i;  // to itself. If j is contracted with another index
		}															// less then i, this mapping will be updated next.

		if ( contractedIndexMapping[ indexPair->j ] < indexPair->i ) {
			contractedIndexMapping[ indexPair->i ] = contractedIndexMapping[ indexPair->j ];
		}

		// Iterate over current contents of the index mapping dictionary and terminate all contractions.
		for ( map<int, int>::iterator key = contractedIndexMapping.begin(); key != contractedIndexMapping.end(); ++key ) {
			contractedIndexMapping[ key->first ] = getTerminatedContraction( contractedIndexMapping, key->second );
		}
	}

	return contractedIndexMapping;
}

Sum distributeTrace( SymbolicTermPtr tr ) {
	if ( tr->getTermID() != TermTypes::TRACE ) {
		cout << "***ERROR: An object other than a trace was passed to distributeTrace()." << endl;
		return Sum();  // TODO: Throw exception.
	}

	// The passed SymbolicTermPtr is a known instance of Trace; cast it.
	TracePtr castTrace = static_pointer_cast<Trace>( tr );

	// Simple case: a single product inside a trace has been passed in.
	if ( castTrace->expr->getTermID() == TermTypes::PRODUCT ) {
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
		if ( expandedArgument->getTermID() == TermTypes::SUM ) {
			return distributeTrace( Trace( expandedArgument ).copy() );

		} else if ( expandedArgument->getTermID() == TermTypes::PRODUCT ){  // Otherwise, if the expression is a simple, single product, proceed.
			Product simplifiedProduct;
			vector<SymbolicTermPtr> matrixTerms;
			ProductPtr castExpandedArgument = static_pointer_cast<Product>( expandedArgument );
			for ( vector<SymbolicTermPtr>::iterator iter = castExpandedArgument->getIteratorBegin(); iter != castExpandedArgument->getIteratorEnd(); ++iter ) {
				if ( (*iter)->getTermID() == TermTypes::COEFFICIENT_FLOAT or (*iter)->getTermID() == TermTypes::COEFFICIENT_FRACTION) {
					simplifiedProduct.addTerm( (*iter)->copy() );
				} else if ( (*iter)->getTermID() == TermTypes::MATRIX_B or (*iter)->getTermID() == TermTypes::MATRIX_M ) {
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
	} else if ( castTrace->expr->getTermID() == TermTypes::SUM ){
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
	if ( expr->getTermID() == TermTypes::PRODUCT ) {
		ProductPtr castProduct = static_pointer_cast<Product>( expr );
		if ( castProduct->containsSum() ) {
			expr = static_pointer_cast<SymbolicTerm>( castProduct->getExpandedExpr().copy() );
		}
	}

	// Reduce the expression tree.
	expr->reduceTree();

	SumPtr distributedExpr( new Sum() );
	if ( expr->getTermID() == TermTypes::SUM ) {
		SumPtr castSum = static_pointer_cast<Sum>( expr );
		for ( vector<SymbolicTermPtr>::iterator iter = castSum->getIteratorBegin(); iter != castSum->getIteratorEnd(); ++iter ) {  // Loop over terms.
			distributedExpr->addTerm( distributeAllTraces( *iter ).copy() );
		}
	} else if ( expr->getTermID() == TermTypes::PRODUCT ) {
		ProductPtr castProduct = static_pointer_cast<Product>( expr );
		Product distributedProduct;
		for ( vector<SymbolicTermPtr>::iterator iter = castProduct->getIteratorBegin(); iter != castProduct->getIteratorEnd(); ++iter  ) {
			if ( (*iter)->getTermID() == TermTypes::TRACE ) {
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
	} else if ( expr->getTermID() == TermTypes::TRACE ) {
		distributedExpr->addTerm( distributeTrace( expr ).copy() ); // need expr.copy() as argument?
	}

	distributedExpr->reduceTree();
	return *distributedExpr;
}

Sum truncateAOrder( SymbolicTermPtr expr, int highestOrder ) {
	if ( expr->getTermID() != TermTypes::SUM ) {
		return Sum();  // TODO: Throw exception.
	}

	Sum truncatedSum;
	SumPtr castExpr = static_pointer_cast<Sum>( expr );

	for ( vector<SymbolicTermPtr>::iterator term = castExpr->getIteratorBegin(); term != castExpr->getIteratorEnd(); ++term ) {
		int orderInA = 0;

		if ( (*term)->getTermID() == TermTypes::PRODUCT ) {
			ProductPtr castTerm = static_pointer_cast<Product>( *term );

			for ( vector<SymbolicTermPtr>::iterator factor = castTerm->getIteratorBegin(); factor != castTerm->getIteratorEnd(); ++factor ) {
				if ( (*factor)->getTermID() == TermTypes::TERM_A ) orderInA++;
			}
		}

		if ( orderInA <= highestOrder ) truncatedSum.addTerm( (*term)->copy() );
	}

	return truncatedSum;
}

Sum truncateOddOrders( SymbolicTermPtr expr ) {
	if ( expr->getTermID() != TermTypes::SUM ) {
		return Sum();  // TODO: Throw exception.
	}

	Sum truncatedSum;
	SumPtr castExpr = static_pointer_cast<Sum>( expr );

	for ( vector<SymbolicTermPtr>::iterator term = castExpr->getIteratorBegin(); term != castExpr->getIteratorEnd(); ++term ) {
		if ( (*term)->getTermID() != TermTypes::PRODUCT ) {
			return Sum();  // TODO: Throw exception.
		}

		ProductPtr castProduct = static_pointer_cast<Product>( *term );
		int orderInA = 0;

		for ( vector<SymbolicTermPtr>::iterator factor = castProduct->getIteratorBegin(); factor != castProduct->getIteratorEnd(); ++factor ) {
			if ( (*factor)->getTermID() == TermTypes::TERM_A ) orderInA++;
		}

		if ( orderInA % 2 == 0 ) {  // If the order of the term in the parameter A is odd, erase the vector element.
			truncatedSum.addTerm( (*term)->copy() );
		}
	}

	return truncatedSum;
}

void rewriteSumInKSFormalism( SymbolicTermPtr expr ) {
	if ( expr->getTermID() != TermTypes::SUM ) {
		// TODO: Throw exception.
		cout << "***ERROR: Expression passed to rewriteSumInKSFormalism() must be an instance of a Sum." << endl;
		return;
	}

	SumPtr castExpr = static_pointer_cast<Sum>( expr );
	for ( vector<SymbolicTermPtr>::iterator iter = castExpr->getIteratorBegin(); iter != castExpr->getIteratorEnd(); ++iter ) {
		if ( (*iter)->getTermID() == TermTypes::PRODUCT ) {
			ProductPtr castTerm = static_pointer_cast<Product>( *iter );
			for ( vector<SymbolicTermPtr>::iterator factor = castTerm->getIteratorBegin(); factor != castTerm->getIteratorEnd(); ++factor ) {
				if ( (*factor)->getTermID() == TermTypes::TRACE) {
					TracePtr castTrace = static_pointer_cast<Trace>( *factor );
					castTrace->rewriteInKSFormalism();
				}
			}
		} else if ( (*iter)->getTermID() == TermTypes::TRACE ) {
			TracePtr castTrace = static_pointer_cast<Trace>( *iter );
			castTrace->rewriteInKSFormalism();
		}
	}
}

void indexExpression( SymbolicTermPtr expr ) {
	if ( expr->getTermID() != TermTypes::SUM ) {
		cout << "***ERROR: indexExpression() expects a Sum as the passed expression." << endl;
		return;
	}

	SumPtr castExpr = static_pointer_cast<Sum>( expr );
	for ( vector<SymbolicTermPtr>::iterator term = castExpr->getIteratorBegin(); term != castExpr->getIteratorEnd(); ++term ) {
		if ( (*term)->getTermID() == TermTypes::PRODUCT ) {
			ProductPtr castTerm = static_pointer_cast<Product>( *term );
			unsigned int nextIndex = 0;

			for ( vector<SymbolicTermPtr>::iterator factor = castTerm->getIteratorBegin(); factor != castTerm->getIteratorEnd(); ++factor ) {
				if ( (*factor)->getTermID() == TermTypes::TRACE ) {
					TracePtr castTrace = static_pointer_cast<Trace>( *factor );

					if ( castTrace->expr->getTermID() != TermTypes::PRODUCT ) {  // TODO: Consider case where a single matrix is the argument to the Trace.
						cout << "***ERROR: Expression passed to indexExpression() contains a badly formed trace. Be sure to distribute all traces first." << endl;
						return;
					}

					ProductPtr castTraceExpr = static_pointer_cast<Product>( castTrace->expr );

					for ( int i = 0; i < castTraceExpr->getNumberOfTerms() - 1; i++ ) {
						castTraceExpr->terms[ i ]->setIndices( nextIndex + i, nextIndex + i + 1 );
					}
					castTraceExpr->terms[ castTraceExpr->getNumberOfTerms() - 1 ]->setIndices( nextIndex + castTraceExpr->getNumberOfTerms() - 1, nextIndex);
					nextIndex += castTraceExpr->getNumberOfTerms() + 1;

					*factor = SymbolicTermPtr( castTraceExpr );
				}
			}

		}
	}
}

Sum fourierTransformExpression( SymbolicTermPtr expr ) {
	Sum transformedExpression;

	if ( expr->getTermID() != TermTypes::SUM) {
		return Sum(); // TODO: Raise exception.
	}

	SumPtr castExpr = static_pointer_cast<Sum>( expr );
	for ( vector<SymbolicTermPtr>::iterator term = castExpr->getIteratorBegin(); term != castExpr->getIteratorEnd(); ++term ) {
		if ( (*term)->getTermID() != TermTypes::PRODUCT ) {
			return Sum();  // TODO: Raise exception.
		}

		ProductPtr castTerm = static_pointer_cast<Product>( *term );

		Product transformedProduct;
		int orderInK = 0;
		map<int, int> deltaMapping;
		vector<IndexContraction> fourierIndices;
		vector<IndexContraction> indexPairsToBeContracted;

		for ( vector<SymbolicTermPtr>::iterator factor = castTerm->getIteratorBegin(); factor != castTerm->getIteratorEnd(); ++factor ) {
			if ( (*factor)->getTermID() == TermTypes::MATRIX_K ) {
				MatrixKPtr castFactor = static_pointer_cast<MatrixK>( *factor );
				castFactor->fourierTransform();

				transformedProduct.addTerm( castFactor->copy() );
				fourierIndices.push_back( IndexContraction( castFactor->getIndices()[0], castFactor->getIndices()[1] ) );

				orderInK++;
			} else if ( (*factor)->getTermID() == TermTypes::DELTA ) {
				DeltaPtr castFactor = static_pointer_cast<Delta>( *factor );
				if ( not castFactor->isDeltaBar() ) {
					indexPairsToBeContracted.push_back( IndexContraction( castFactor->getIndices()[0], castFactor->getIndices()[1] ) );
				}
			} else {
				transformedProduct.addTerm( (*factor)->copy() );
			}
		}
	}
}

/*
 * ***********************************************************************
 * INPUT REDIRECTION OPERATOR OVERLOADS
 * ***********************************************************************
 */

std::ostream& operator<<( std::ostream& os, const TermA &obj ) {
	os << obj.to_string();
	return os;
}

std::ostream& operator<<( std::ostream& os, const MatrixM &obj ) {
	os << obj.to_string();
	return os;
}

std::ostream& operator<<( std::ostream& os, const TermTypes &obj ) {
	os << static_cast<char>( obj );
	return os;
}

std::ostream& operator<<( std::ostream& os, const map<int, int> &obj ) {
	os << "[";
	for ( map<int, int>::const_iterator pair = obj.begin(); pair != obj.end(); ++pair ) {
		os << " " << pair->first << " : " << pair->second << " ";
	}
	os << "]";
}