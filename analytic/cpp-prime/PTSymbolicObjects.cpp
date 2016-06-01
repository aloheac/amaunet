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
#include <set>
#include <boost/algorithm/string/trim.hpp>
#include "PTSymbolicObjects.h"
#include "PathIntegration.h"

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
	flavorLabel = "";
	indices[ 0 ] = 0;
	indices[ 1 ] = 0;
	termID = TermTypes::INVALID_TERM ;
}

SymbolicTerm::SymbolicTerm( const SymbolicTerm &other ) : flavorLabel( other.flavorLabel ), termID( other.termID ) {
	indices[ 0 ] = other.indices[ 0 ];
	indices[ 1 ] = other.indices[ 1 ];
}

SymbolicTerm::~SymbolicTerm() { }

SymbolicTerm& SymbolicTerm::operator=( const SymbolicTerm& rhs ) {
	flavorLabel = rhs.flavorLabel;
	termID = rhs.termID;
	indices[ 0 ] = rhs.indices[ 0 ];
	indices[ 1 ] = rhs.indices[ 1 ];

	return *this;
}

bool SymbolicTerm::operator==( const SymbolicTerm &other ) const {
	return ( flavorLabel == other.flavorLabel );
}

ostream& operator<<( ostream& os, const SymbolicTerm &st ) {
	stringstream ss;
	ss << st.to_string();
	os << ss.str();
	return os;
}

void SymbolicTerm::simplify() { }

const string SymbolicTerm::to_string() const {
	return "<invalid_term>";
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

MatrixK::MatrixK( const MatrixK &other ) : SymbolicTerm( other ) {
	isFourierTransformed = other.isFourierTransformed;
}

MatrixK& MatrixK::operator=( const MatrixK &rhs ) {
	SymbolicTerm::operator=( rhs );
	isFourierTransformed = rhs.isFourierTransformed;

	return *this;
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
	cpy->flavorLabel = flavorLabel;
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
	cpy->flavorLabel = flavorLabel;
	cpy->indices[0] = indices[0];
	cpy->indices[1] = indices[1];
	cpy->termID = TermTypes::MATRIX_S;

	return static_pointer_cast<SymbolicTerm>( cpy );
}

/*
 * TermA
 */

TermA::TermA() : SymbolicTerm() {
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
 * TermE
 */

TermE::TermE( int thisOrder ) : order( thisOrder ) { }

TermE::TermE( int thisOrder, const char* thisFlavorLabel ) : order( thisOrder ) {
	flavorLabel = thisFlavorLabel;
}

const std::string TermE::to_string() const {
	stringstream ss;
	ss << "E" << order;

	if ( flavorLabel != "" ) ss << "_" << flavorLabel;

	return ss.str();
}

bool TermE::operator==( const TermE &other ) const {
	return ( order == other.order ) and ( flavorLabel == other.flavorLabel );
}

SymbolicTermPtr TermE::copy() {
	return SymbolicTermPtr( new TermE( order, flavorLabel ) );
}

unsigned int TermE::getOrder() {
	return order;
}

SymbolicTermPtr TermE::getFullExpression() {
	Product expr;

	int sign = 1;
	if ( order % 2 == 0 ) sign = -1;

	expr.addTerm( CoefficientFractionPtr( new CoefficientFraction( sign, order ) ) );

	Product traceArgument;
	for ( int i = 0; i < order; i++ ) {
		traceArgument.addTerm( MatrixKPtr( new MatrixK( flavorLabel ) ) );
		traceArgument.addTerm( MatrixSPtr( new MatrixS() ) );
	}

	expr.addTerm( TracePtr( new Trace( traceArgument.copy() ) ) );

	return expr.copy();
}

/*
 * CoefficientFloat;
 */

CoefficientFloat::CoefficientFloat( double val ) : SymbolicTerm() {
	value = val;
	termID = TermTypes::COEFFICIENT_FLOAT;
}

const string CoefficientFloat::to_string() const {
	stringstream ss;
	if ( abs( value ) < 1E-10 ) {
		ss << "0";
	} else {
		ss << value;
	}

	return ss.str();
}

CoefficientFloat CoefficientFloat::operator*( const CoefficientFloat& obj ) const {
	return CoefficientFloat( value * obj.value );
}

CoefficientFloat CoefficientFloat::operator+( const CoefficientFloat& obj ) const {
	return CoefficientFloat( value + obj.value );
}

CoefficientFraction CoefficientFloat::operator*( const CoefficientFraction& obj ) const {
	return obj * (*this);
}

CoefficientFraction CoefficientFloat::operator+( const CoefficientFraction& obj ) const {
	return obj + (*this);
}

CoefficientFloat& CoefficientFloat::operator*=( const CoefficientFloat& obj ) {
	this->value *= obj.value;
	return *this;
}

CoefficientFloat& CoefficientFloat::operator+=( const CoefficientFloat& obj ) {
	this->value += obj.value;
	return *this;
}

SymbolicTermPtr CoefficientFloat::copy() {
	SymbolicTermPtr cpy( new CoefficientFloat( value ) );
	return cpy;
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

CoefficientFraction::CoefficientFraction( double n, double d ) : SymbolicTerm() {
	num = n;
	den = d;
	termID = TermTypes::COEFFICIENT_FRACTION;
}

const string CoefficientFraction::to_string() const {
	stringstream ss;
	ss << num << " / " << den;
	return ss.str();
}

CoefficientFraction CoefficientFraction::operator*( const CoefficientFraction& obj ) const {
	CoefficientFraction product( num * obj.num, den * obj.den );
	product.reduce();
	return product;
}

CoefficientFraction CoefficientFraction::operator+( const CoefficientFraction& obj ) const {
	CoefficientFraction sum( num * obj.den + den * obj.num, den * obj.den );
	sum.reduce();
	return sum;
}

CoefficientFraction CoefficientFraction::operator*( const CoefficientFloat& obj ) const {
	CoefficientFraction product( num * obj.value, den );
	product.reduce();
	return product;
}

CoefficientFraction CoefficientFraction::operator+( const CoefficientFloat& obj ) const {
	CoefficientFraction sum( num + obj.value * den, den );
	sum.reduce();
	return sum;
}

CoefficientFraction& CoefficientFraction::operator*=( const CoefficientFraction& obj ) {
	CoefficientFraction product = (*this) * obj;
	this->num = product.num;
	this->den = product.den;

	return (*this);
}

CoefficientFraction& CoefficientFraction::operator+=( const CoefficientFraction& obj ) {
	CoefficientFraction sum = (*this) + obj;
	this->num = sum.num;
	this->den = sum.den;

	return (*this);
}

SymbolicTermPtr CoefficientFraction::copy() {
	SymbolicTermPtr cpy( new CoefficientFraction( num, den ) );
	return cpy;
}

double CoefficientFraction::eval() {
	return num / den;
}


void CoefficientFraction::reduce() {
	if ( abs( floor( num ) - num ) == 0 and abs( floor( den ) - den ) == 0 ) {
		double thisGCD = (double)gcd( (unsigned long)abs( num ), (unsigned long)abs( den ) );
		num /= thisGCD;
		den /= thisGCD;
	}
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
		if ( trimmedRepresentation == "0" or trimmedRepresentation == "{0}" or trimmedRepresentation == "0 / 1" or trimmedRepresentation == "{0 / 1}" or isZeroTrace( *iter ) ) {
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

Sum Sum::getExpandedExpr() {
	Sum expandedSum;
	for ( vector<SymbolicTermPtr>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		if ( (*iter)->getTermID() == TermTypes::SUM ) {
			SumPtr s = dynamic_pointer_cast<Sum>( *iter );
			expandedSum.addTerm( s->getExpandedExpr().copy() );
		} else if ( (*iter)->getTermID() == TermTypes::PRODUCT ) {
			ProductPtr prod = dynamic_pointer_cast<Product>( *iter );
			expandedSum.addTerm( prod->getExpandedExpr().copy() );
		} else {
			expandedSum.addTerm( (*iter)->copy() );
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

void Sum::clear() {
	terms.clear();
}

void Sum::combineCoefficients() {
	for ( vector<SymbolicTermPtr>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		if ( (*iter)->getTermID() == TermTypes::PRODUCT ) {
			ProductPtr castTerm = static_pointer_cast<Product>( *iter );
			castTerm->combineCoefficients();
		}
	}
}

bool Sum::operator==( const Sum &other ) const {
	return false; // TODO
}

void Sum::reduceFourierSumIndices() {
	for ( vector<SymbolicTermPtr>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		if ( (*iter)->getTermID() == TermTypes::PRODUCT ) {
			ProductPtr castTerm = static_pointer_cast<Product>( *iter );
			castTerm->reduceFourierSumIndices();
		}
	}
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
	for ( vector<SymbolicTermPtr>::iterator iter =  terms.begin(); iter != terms.end(); ) {
		(*iter)->simplify();
		unpackTrivialExpression( *iter );

		string trimmedRepresentation = (*iter)->to_string();  // To be trimmed on next statement.
		boost::trim( trimmedRepresentation );

		if ( trimmedRepresentation == "0" or trimmedRepresentation == "-0" or trimmedRepresentation == "0 / 1" or isZeroTrace( *iter ) ) {
			CoefficientFloatPtr zero( new CoefficientFloat( 0.0 ) );
			vector<SymbolicTermPtr> zeroVector;
			zeroVector.push_back( zero );
			// delete terms?
			terms = zeroVector;
			return;
		} else if ( trimmedRepresentation == "1" or trimmedRepresentation == "1 / 1" ) {
			if ( getNumberOfTerms() != 1 ) {  // We don't want to delete the last remaining factor if it is one.
				(*iter).reset(); // Verify.
				iter = terms.erase(iter);

				// Note that we could have the case where the last element of the vector could be removed, and the call
				// vector::erase() returns the iterator terms.end(). The for loop then advances the iterator, which
				// becomes invalid and results in a seg fault on the next recursive call. Therefore, check if the
				// iterator has reached the end of the data structure, and if so, break out of the loop.
				if (iter == terms.end()) break;
			} else {
				++iter;
			}
		} else {
			++iter;
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

bool Product::operator==( const Sum &other ) const {
	return false;
}

bool Product::containsSum() {
	for ( vector<SymbolicTermPtr>::iterator iter = terms.begin(); iter != terms.end(); ++iter ) {
		if ( (*iter)->getTermID() == TermTypes::SUM ) return true;
	}
	return false;
}

void Product::zero() {
	vector<SymbolicTermPtr> zero;
	zero.push_back( CoefficientFloatPtr( new CoefficientFloat( 0.0 ) ) );
	terms = zero;
}

void Product::clear() {
	terms.clear();
}

void Product::reduceFourierSumIndices() {
	for ( vector<SymbolicTermPtr>::iterator factor = terms.begin(); factor != terms.end(); ++factor ) {
		if ( (*factor)->getTermID() == TermTypes::FOURIER_SUM ) {
			FourierSumPtr castFactor = static_pointer_cast<FourierSum>( *factor );
			castFactor->reduceDummyIndices();
		}
	}
}

void Product::combineCoefficients() {
	CoefficientFraction runningCoefficient( 1, 1 );
	for ( vector<SymbolicTermPtr>::iterator factor = terms.begin(); factor != terms.end(); ) {
		if ( (*factor)->getTermID() == TermTypes::COEFFICIENT_FRACTION ) {
			runningCoefficient *= *(static_pointer_cast<CoefficientFraction>( *factor ));
			factor = terms.erase( factor );
		} else if ( (*factor)->getTermID() == TermTypes::COEFFICIENT_FLOAT ) {
			runningCoefficient = runningCoefficient * ( *(static_pointer_cast<CoefficientFloat>( *factor )) );
			factor = terms.erase( factor );
		} else {
			++factor;
		}
	}

	terms.push_back( runningCoefficient.copy() );
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

bool Trace::operator==( const Trace &other ) const {
	return *expr == other;
}

bool Trace::operator<( const Trace &other ) const {
	int thisNumTerms = 0;
	if ( expr->getTermID() == TermTypes::PRODUCT ) {
		ProductPtr castProduct = static_pointer_cast<Product>( expr );
		thisNumTerms = castProduct->getNumberOfTerms();
	} else if ( expr->getTermID() == TermTypes::SUM ) {
		cout << "***WARNING: Traces were not distributed prior to a call to Trace::operator<. Results cannot be trusted." << endl;
		SumPtr castSum = static_pointer_cast<Sum>( expr );
		thisNumTerms = castSum->getNumberOfTerms();
	}

	int otherNumTerms = 0;
	if ( other.expr->getTermID() == TermTypes::PRODUCT ) {
		ProductPtr castProduct = static_pointer_cast<Product>( other.expr );
		otherNumTerms = castProduct->getNumberOfTerms();
	} else if ( other.expr->getTermID() == TermTypes::SUM ) {
		cout << "***WARNING: Traces were not distributed prior to a call to Trace::operator<. Results cannot be trusted." << endl;
		SumPtr castSum = static_pointer_cast<Sum>( expr );
		otherNumTerms = castSum->getNumberOfTerms();
	}

	return thisNumTerms < otherNumTerms;
}

bool Trace::operator>( const Trace &other ) const {
	int thisNumTerms = 0;
	if ( expr->getTermID() == TermTypes::PRODUCT ) {
		ProductPtr castProduct = static_pointer_cast<Product>( expr );
		thisNumTerms = castProduct->getNumberOfTerms();
	} else if ( expr->getTermID() == TermTypes::SUM ) {
		cout << "***WARNING: Traces were not distributed prior to a call to Trace::operator>. Results cannot be trusted." << endl;
		SumPtr castSum = static_pointer_cast<Sum>( expr );
		thisNumTerms = castSum->getNumberOfTerms();
	}

	int otherNumTerms = 0;
	if ( other.expr->getTermID() == TermTypes::PRODUCT ) {
		ProductPtr castProduct = static_pointer_cast<Product>( other.expr );
		otherNumTerms = castProduct->getNumberOfTerms();
	} else if ( other.expr->getTermID() == TermTypes::SUM ) {
		cout << "***WARNING: Traces were not distributed prior to a call to Trace::operator>. Results cannot be trusted." << endl;
		SumPtr castSum = static_pointer_cast<Sum>( expr );
		otherNumTerms = castSum->getNumberOfTerms();
	}

	return thisNumTerms > otherNumTerms;
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

SymbolicTermPtr FourierSum::copy() {
	return FourierSumPtr( new FourierSum( indices, order ) );
}

void FourierSum::reduceDummyIndices() {
	// Reduce all single loop contractions (a, a) to (0, 0).
	for ( vector<IndexContraction>::iterator indexPair = indices.begin(); indexPair != indices.end(); ++indexPair ) {
		if ( indexPair->i == indexPair->j ) {
			indexPair->i = 0;
			indexPair->j = 0;
		}
	}

	set<int> presentIndices;
	for ( vector<IndexContraction>::iterator indexPair = indices.begin(); indexPair != indices.end(); ++indexPair ) {
		presentIndices.insert( indexPair->i );
		presentIndices.insert( indexPair->j );
	}

	map<int, int> reducedIndexMapping;
	int mappedIndex = 0;
	for ( set<int>::iterator dummyIndex = presentIndices.begin(); dummyIndex != presentIndices.end(); ++dummyIndex ) {
		reducedIndexMapping[ *dummyIndex ] = mappedIndex;
		mappedIndex++;
	}

	for ( vector<IndexContraction>::iterator indexPair = indices.begin(); indexPair != indices.end(); ++indexPair ) {
		indexPair->i = reducedIndexMapping[ indexPair->i ];
		indexPair->j = reducedIndexMapping[ indexPair->j ];
	}
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

int getProductAOrder( SymbolicTermPtr prod ) {
	if ( prod->getTermID() != TermTypes::PRODUCT ) {
		cout << "***ERROR: An expression other than a Product was passed to getProductAOrder()." << endl;
		return false;  // TODO: Throw exception.
	}

	ProductPtr castTerm = static_pointer_cast<Product>( prod );

	int orderInA = 0;


	for ( vector<SymbolicTermPtr>::iterator factor = castTerm->getIteratorBegin(); factor != castTerm->getIteratorEnd(); ++factor ) {
		if ( (*factor)->getTermID() == TermTypes::TERM_A ) orderInA++;
	}

	return orderInA;
}

int getTerminatedContraction( map<int, int> contractedIndexMapping, int index ) {
	if ( contractedIndexMapping.count( index ) == 0 or contractedIndexMapping[ index ] == index  ) {
		return index;
	} else {
		return getTerminatedContraction( contractedIndexMapping, contractedIndexMapping[ index ] );
	}
}

vector< set<int> > groupContractions( vector< set<int> > groups ) {
	for ( vector< set<int> >::iterator group = groups.begin(); group != groups.end(); ++group ) {
		for ( set<int>::iterator index = group->begin(); index != group->end(); ++index ) {
			for ( vector< set<int> >::iterator comparedGroup = groups.begin(); comparedGroup != groups.end(); ++comparedGroup ) {

				if ( group != comparedGroup and comparedGroup->count( *index ) > 0 ) {
					for ( set<int>::iterator indexToAppend = comparedGroup->begin(); indexToAppend != comparedGroup->end(); ++indexToAppend ) {
						group->insert( *indexToAppend );
					}

					groups.erase( comparedGroup );
					return groupContractions( groups );
				}

			}
		}
	}

	return groups;
}

map<int, int> constructContractionDictionary( DeltaContractionSet contractions ) {
	contractions.orderContractionIndices();
	contractions.sortContractions();

	vector< set<int> > contractionIndexSet;
	for ( vector<IndexContraction>::iterator indexPair = contractions.getIteratorBegin(); indexPair != contractions.getIteratorEnd(); ++indexPair ) {
		set<int> nextIndexPair;
		nextIndexPair.insert( indexPair->i );
		nextIndexPair.insert( indexPair->j );
		contractionIndexSet.push_back( nextIndexPair );
	}

	contractionIndexSet = groupContractions( contractionIndexSet );

	map<int, int> contractionDictionary;
	for ( vector< set<int> >::iterator group = contractionIndexSet.begin(); group != contractionIndexSet.end(); ++group ) {
		int smallestGroupIndex = *(min_element( group->begin(), group->end() ));

		for( set<int>::iterator index = group->begin(); index != group->end(); ++index ) {
			contractionDictionary[ *index ] = smallestGroupIndex;
		}
	}

	return contractionDictionary;
}

bool areTermsCommon( SymbolicTermPtr termA, SymbolicTermPtr termB ) {
	if ( termA->getTermID() != TermTypes::PRODUCT or termB->getTermID() != TermTypes::PRODUCT ) {
		cout << "***ERROR: An expression that is not a Product was passed to areTermsCommon()." << endl;
		return false;  // TODO: Raise exception.
	}

	// Check if order in A in both terms are the same.
	if ( getProductAOrder( termA ) != getProductAOrder( termB ) ) return false;

	// Check if FourierSum signatures for both terms are the same.
	SymbolicTermPtr termAFourierSum, termBFourierSum;
	ProductPtr castTermA = static_pointer_cast<Product>( termA );
	ProductPtr castTermB = static_pointer_cast<Product>( termB );

	for ( vector<SymbolicTermPtr>::iterator factor = castTermA->getIteratorBegin(); factor != castTermA->getIteratorEnd(); ++factor ) {
		if ( (*factor)->getTermID() == TermTypes::FOURIER_SUM ) {
			termAFourierSum = *factor;
			break;
		}
	}

	for ( vector<SymbolicTermPtr>::iterator factor = castTermB->getIteratorBegin(); factor != castTermB->getIteratorEnd(); ++factor ) {
		if ( (*factor)->getTermID() == TermTypes::FOURIER_SUM ) {
			termBFourierSum = *factor;
			break;
		}
	}

	// Make sure we have valid pointers for FourierSum objects.
	if ( termAFourierSum == nullptr or termBFourierSum == nullptr or termAFourierSum->getTermID() != TermTypes::FOURIER_SUM or termBFourierSum->getTermID() != TermTypes::FOURIER_SUM ) return false;

	FourierSumPtr castTermAFourierSum = static_pointer_cast<FourierSum>( termAFourierSum );
	FourierSumPtr castTermBFourierSum = static_pointer_cast<FourierSum>( termBFourierSum );

	if ( not( *castTermAFourierSum == *castTermBFourierSum ) ) return false;

	return true;
}

Sum truncateAOrder( SymbolicTermPtr expr, int highestOrder ) {
	if ( expr->getTermID() != TermTypes::SUM ) {
		cout << "***ERROR: An expression other than a Sum was passed to truncateAOrders()." << endl;
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
		cout << "***ERROR: An expression other than a Sum was passed to truncateOddOrders()." << endl;
		return Sum();  // TODO: Throw exception.
	}

	Sum truncatedSum;
	SumPtr castExpr = static_pointer_cast<Sum>( expr );

	for ( vector<SymbolicTermPtr>::iterator term = castExpr->getIteratorBegin(); term != castExpr->getIteratorEnd(); ++term ) {
		if ( (*term)->getTermID() != TermTypes::PRODUCT ) {
			if ( (*term)->to_string() != "0" and (*term)->to_string() != "1"  and (*term)->to_string() != "1 / 0 "  and (*term)->to_string() != "1 / 1" ) {
				cout << "***WARNING: (WB1) A term other then a product, zero, or one was encountered when truncating odd orders. The solution may still be correct, but should be inspected." << endl;
			}

			(*term) = Product( *term ).copy();
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
					nextIndex += castTraceExpr->getNumberOfTerms();

					*factor = SymbolicTermPtr( castTraceExpr );
				}
			}

		}
	}
}

Sum fourierTransformExpression( SymbolicTermPtr expr ) {
	Sum transformedExpression;

	if ( expr->getTermID() != TermTypes::SUM) {
		cout << "***ERROR: An expression other than a Sum was passed to fourierTransformExpression()." << endl;
		return Sum(); // TODO: Raise exception.
	}

	SumPtr castExpr = static_pointer_cast<Sum>( expr );
	for ( vector<SymbolicTermPtr>::iterator term = castExpr->getIteratorBegin(); term != castExpr->getIteratorEnd(); ++term ) {
		if ( (*term)->getTermID() != TermTypes::PRODUCT ) {
			if ( (*term)->to_string() != "0" and (*term)->to_string() != "1"  and (*term)->to_string() != "1 / 0 "  and (*term)->to_string() != "1 / 1" ) {
				cout << "***WARNING: (WC1) A term other then a product, zero, or one was encountered when computing symbolic Fourier transform. The solution may still be correct, but should be inspected." << endl;
			}

			(*term) = Product( *term ).copy();
		}

		ProductPtr castTerm = static_pointer_cast<Product>( *term );

		Product transformedProduct;
		int orderInK = 0;
		vector<IndexContraction> fourierIndices;
		DeltaContractionSet indexPairsToBeContracted;

		for ( vector<SymbolicTermPtr>::iterator factor = castTerm->getIteratorBegin(); factor != castTerm->getIteratorEnd(); ++factor ) {
			if ( (*factor)->getTermID() == TermTypes::MATRIX_K ) {
				MatrixKPtr castFactor = static_pointer_cast<MatrixK>( *factor );
				castFactor->fourierTransform();

				// TODO: Update MatrixK indices with momentum-space index.

				transformedProduct.addTerm( castFactor->copy() );
				fourierIndices.push_back( IndexContraction( castFactor->getIndices()[0], castFactor->getIndices()[1] ) );

				orderInK++;
			} else if ( (*factor)->getTermID() == TermTypes::DELTA ) {
				DeltaPtr castFactor = static_pointer_cast<Delta>( *factor );
				if ( not castFactor->isDeltaBar() ) {
					indexPairsToBeContracted.addContraction( IndexContraction( castFactor->getIndices()[0], castFactor->getIndices()[1] ) );
				}
			} else {
				transformedProduct.addTerm( (*factor)->copy() );
			}
		}

		map<int, int> indexDictionary = constructContractionDictionary( indexPairsToBeContracted );
		if ( orderInK > 0 ) {
			for ( map<int, int>::iterator pair = indexDictionary.begin(); pair != indexDictionary.end(); ++pair ) {
				for ( vector<IndexContraction>::iterator contraction = fourierIndices.begin(); contraction != fourierIndices.end(); ++contraction ) {
					contraction->i = indexDictionary[ contraction->i ];
					contraction->j = indexDictionary[ contraction->j ];
				}
			}

			transformedProduct.addTerm( FourierSumPtr( new FourierSum( fourierIndices, orderInK ) ) );
		}

		transformedExpression.addTerm( transformedProduct.copy() );
	}

	return transformedExpression;
}

Sum combineLikeTerms( Sum &expr ) {
	Sum reducedSum;

	expr.simplify();

	for ( vector<SymbolicTermPtr>::iterator term = expr.getIteratorBegin(); term != expr.getIteratorEnd(); ++term ) {
		CoefficientFraction runningLikeTermsCoefficient( 0, 1 );

		CoefficientFraction termATotalCoefficient( 1, 1 );
		Product combinedFactorA;

		if ( (*term)->getTermID() != TermTypes::PRODUCT ) {
			// Under the current perturbation theory formalism, an expanded expression in standard form should always
			// be a Sum of Product instances, except for the case where the expression is identically zero. We consider
			// and account for that case here. The expression was mathematically simplified before entering the loop.
			if ( (*term)->to_string() != "0" and (*term)->to_string() != "1"  and (*term)->to_string() != "1 / 0 "  and (*term)->to_string() != "1 / 1" ) {
				cout << "***WARNING: (WA1) A term other then a product, zero, or one was encountered when combining like terms. The solution may still be correct, but should be inspected." << endl;
			}

			(*term) = Product( *term ).copy();
		}

		ProductPtr castTerm = static_pointer_cast<Product>( *term );

		for ( vector<SymbolicTermPtr>::iterator factor = castTerm->getIteratorBegin(); factor != castTerm->getIteratorEnd(); ++factor ) {
			if ( (*factor)->getTermID() == TermTypes::COEFFICIENT_FLOAT ) {
				CoefficientFloatPtr castFactor = static_pointer_cast<CoefficientFloat>( *factor );
				termATotalCoefficient *= CoefficientFraction( castFactor->eval(), 1 );
			} else if ( (*factor)->getTermID() == TermTypes::COEFFICIENT_FRACTION ) {
				CoefficientFractionPtr castFactor = static_pointer_cast<CoefficientFraction>( *factor );
				termATotalCoefficient *= (*castFactor);
			} else {
				combinedFactorA.addTerm( (*factor)->copy() );
			}
		}

		runningLikeTermsCoefficient = runningLikeTermsCoefficient + termATotalCoefficient;

		for ( vector<SymbolicTermPtr>::iterator secondTerm = expr.getIteratorBegin(); secondTerm != expr.getIteratorEnd(); ++secondTerm ) {

			if ( (*secondTerm)->getTermID() != TermTypes::PRODUCT ) {
				// Under the current perturbation theory formalism, an expanded expression in standard form should always
				// be a Sum of Product instances, except for the case where the expression is identically zero. We consider
				// and account for that case here. The expression was mathematically simplified before entering the loop.
				if ( (*secondTerm)->to_string() != "0" and (*secondTerm)->to_string() != "1"  and (*secondTerm)->to_string() != "1 / 0 "  and (*secondTerm)->to_string() != "1 / 1" ) {
					cout << "***WARNING: (WA2) A term other then a product, zero, or one was encountered when combining like terms. The solution may still be correct, but should be inspected." << endl;
				}

				(*secondTerm) = Product( *secondTerm ).copy();
			}

			ProductPtr castSecondTerm = static_pointer_cast<Product>( *secondTerm );

			if ( term != secondTerm and areTermsCommon( *term, *secondTerm ) ) {
				CoefficientFraction termBCoefficient( 1, 1 );

				for ( vector<SymbolicTermPtr>::iterator secondFactor = castSecondTerm->getIteratorBegin(); secondFactor != castSecondTerm->getIteratorEnd(); ++secondFactor ) {
					if ( (*secondFactor)->getTermID() == TermTypes::COEFFICIENT_FLOAT ) {
						CoefficientFloatPtr castSecondFactor = static_pointer_cast<CoefficientFloat>( *secondFactor );
						termBCoefficient *= CoefficientFraction( castSecondFactor->eval(), 1 );
					} else if ( (*secondFactor)->getTermID() == TermTypes::COEFFICIENT_FRACTION ) {
						CoefficientFractionPtr castSecondFactor = static_pointer_cast<CoefficientFraction>( *secondFactor );
						termBCoefficient *= (*castSecondFactor);
					}
				}

				runningLikeTermsCoefficient += termBCoefficient;
				(*castSecondTerm).zero();
			}
		}

		combinedFactorA.addTerm( runningLikeTermsCoefficient.copy() );
		reducedSum.addTerm( combinedFactorA.copy() );
	}

	return reducedSum;
}

Sum combineLikeTerms( Sum &expr, int groupSize ) {
	expr.simplify();  // Remove any zero terms in the sum before proceeding.

	Sum simplifiedExpression;

	if ( expr.getNumberOfTerms() <= groupSize ) {
		simplifiedExpression = combineLikeTerms( expr );
		simplifiedExpression.simplify();
		return simplifiedExpression;
	} else {
		int i = 0;
		while ( i < expr.getNumberOfTerms() ) {
			cout << ">> Combining like terms for batch range " << i << " to " << i + groupSize << " of " << expr.getNumberOfTerms() << "." << endl;

			vector<SymbolicTermPtr>::iterator endingIterator;
			if ( i + groupSize > expr.getNumberOfTerms() ) {
				endingIterator = expr.getIteratorEnd();
			} else {
				endingIterator = expr.getIteratorBegin();
				endingIterator += i + groupSize;
			}

			vector<SymbolicTermPtr> nextGroupOfTerms( expr.getIteratorBegin() += i, endingIterator );
			Sum nextGroup( nextGroupOfTerms );
			simplifiedExpression.addTerm( combineLikeTerms( nextGroup, groupSize ).copy() );
			i += groupSize;
		}

		simplifiedExpression.reduceTree();
		simplifiedExpression.simplify();
		return combineLikeTerms( simplifiedExpression, groupSize );
	}
}

Sum generateExponentialSeries( int order, Product x ) {
	Sum series;

	series.addTerm( CoefficientFloatPtr( new CoefficientFloat( 1 ) ) );

	for ( int i = 0; i < order; i++ ) {
		Product nextTerm;

		nextTerm.addTerm( CoefficientFractionPtr( new CoefficientFraction( 1, factorial( i + 1 ) ) ) );

		for ( int j = 0; j < i + 1; j++ ) {
			nextTerm.addTerm( x.copy() );
		}

		series.addTerm( nextTerm.copy() );
	}

	series.reduceTree();
	series.simplify();

	return series;
}

Sum generateDeterminantExpansion( int order, const char* flavorLabel,  bool insertFullE ) {
	Product expansion;

	for ( int i = 1; i < order + 1; i++ ) {
		Product nextExpansion;

		for ( int j = 0; j < i; j++ ) {
			nextExpansion.addTerm( TermAPtr( new TermA() ) );
		}

		if ( insertFullE ) {
			nextExpansion.addTerm( TermE( i, flavorLabel ).getFullExpression() );
		} else {
			nextExpansion.addTerm( TermE( i, flavorLabel ).copy() );
		}

		Sum series = generateExponentialSeries( (int)ceil( order / i ), nextExpansion );
		expansion.addTerm( series.copy() );
	}

	expansion.reduceTree();
	return Sum( expansion.copy() );
}

unsigned long gcd( unsigned long a, unsigned long b ) {
	if ( a == 0 ) {
		return b;
	} else {
		return gcd( b % a, a );
	}
}

int factorial( int n ) {
	if ( n == 0 or n == 1 ) {
		return 1;
	} else {
		return n * factorial( n - 1 );
	}
}

Sum sortTracesByOrder( Sum &expr ) {
	Sum sortedExpression;

	for ( vector<SymbolicTermPtr>::iterator term = expr.getIteratorBegin(); term != expr.getIteratorEnd(); ++term ) {
		if ( (*term)->getTermID() == TermTypes::PRODUCT ) {
			ProductPtr castProduct = static_pointer_cast<Product>( *term );

			vector<Trace> traces;
			Product sortedProduct;
			for ( vector<SymbolicTermPtr>::iterator factor = castProduct->getIteratorBegin(); factor != castProduct->getIteratorEnd(); ++factor ) {
				if ( (*factor)->getTermID() == TermTypes::TRACE ) {
					TracePtr castTrace = static_pointer_cast<Trace>( *factor );
					traces.push_back( *castTrace );
				} else {
					sortedProduct.addTerm( (*factor)->copy() );
				}
			}

			sort( traces.begin(), traces.end() );

			for ( vector<Trace>::iterator tr = traces.begin(); tr != traces.end(); ++tr ) {
				sortedProduct.addTerm( tr->copy() );
			}

			sortedExpression.addTerm( sortedProduct.copy() );

		} else {
			sortedExpression.addTerm( (*term)->copy() );
		}
	}

	return sortedExpression;
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
	return os;
}