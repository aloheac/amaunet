#include "PTSymbolicObjects.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

using namespace std;

/*
 * *****************************************************************************
 * OBJECT CONSTRUCTORS
 * ******************************************************************************
 */

/*
 * Primary constructor for TermA.
 * @param initOrder Order (power) of A.
 */
TermA::TermA( int initOrder ) : Factor() {
	order = initOrder;
	stringRepresentation = "{ A^" + pt_util::str( initOrder ) + " }";
}

/*
 * Primary constructor for TermD. Note that since D is in momentum space and
 * is diagonal, D is indexed by a single index; hence initIndices should have
 * a length of one.
 * @param initIndices Unique momentum-space index for this matrix.
 */
TermD::TermD( vector<int> initIndices ) : Factor() {
	indices = initIndices;
	stringRepresentation = "{ D_" + pt_util::str( initIndices ) + " }";
}

Expression::Expression( int initNumTerms, Factor* initTerms ) {
	numTerms = initNumTerms;
	terms = initTerms;
}

Expression::Expression() { }

/*
 * Default constructor for ExpressionInterpreter. Takes no arguments.
 */
ExpressionInterpreter::ExpressionInterpreter() { }

/*
 * Default constructor for an empty Sum.
 */
Sum::Sum() {
	products = vector<Product>();
}

/*
 * Default constructor for an empty Product.
 */
Product::Product() {
	factors = vector<Factor>();
}

Factor::Factor() {
	indices = vector<int>();
	numIndices = 0;
	stringRepresentation = "{ <invalid_factor> }";
}

Factor::Factor( vector<int> initIndices ) {
	indices = initIndices;
	numIndices = indices.size();
	stringRepresentation = "{ <invalid_factor> }";
}

FourierSum::FourierSum( vector<int> initIndices ) : Factor( initIndices ) {

}

CoefficientFloat::CoefficientFloat( double initValue ) : Factor() {
	value = initValue;
}

/*
 * *****************************************************************************
 * CLASS METHOD IMPLEMENTATIONS
 * ******************************************************************************
 */

std::string Factor::to_string() const {
	return stringRepresentation;
}

void Factor::setFlavorLabel( string label ) {
	flavorLabel = label;
}

string Factor::getFlavorLabel() {
	return flavorLabel;
}

std::string TermA::to_string() const {
	return "A^" + pt_util::str( order );
}

std::string TermD::to_string() const {
	return "D_" + pt_util::str( indices );
}

std::string FourierSum::to_string() const {
	std::string s = "FourierSum[ ";
	for ( int i = 0; i < numIndices; i++ ) {
		s += pt_util::str( indices[i] );
	}

	return s;
}

double Factor::eval( vector<int> indices, pt_util::PTSystemParameters params ) {
	return 0.0;
}

double TermA::eval( vector<int> indices, pt_util::PTSystemParameters params ) {
	return std::sqrt( std::exp( params.TAU * params.BARE_COUPLING ) - 1.0 );
}

double TermD::eval( vector<int> indices, pt_util::PTSystemParameters params ) {
	return 0.0;
}

double CoefficientFloat::eval( vector<int> indices, pt_util::PTSystemParameters params ) {
	return value;
}

string CoefficientFloat::to_string() const {
	return pt_util::str( value );
}

void Sum::addProduct( Product p ) {
	products.push_back( p );
}

vector<Product>::const_iterator Sum::getIterator() const {
	return products.begin();
}

void Product::addFactor( Factor f ) {
	factors.push_back( f );
}

vector<Factor>::const_iterator Product::getIterator() const {
	return factors.begin();
}

std::ostream& operator<<( std::ostream& os, const Factor& f ) {
	os << f.to_string();
	return os;
}

std::ostream& operator<<( std::ostream& os, const Product& prod ) {
	os << prod.to_string();
	return os;
}

ostream& operator<<( ostream& os, const Sum& s ) {
	os << s.to_string();
	return os;
}

string Product::to_string() const {
	stringstream ss;
	for ( vector<Factor>::const_iterator iter = getIterator(); iter != factors.end(); ++iter ) {
		ss << (*iter) << " ";
	}

	return ss.str();
}

string Sum::to_string() const {
	stringstream ss;
	for( vector<Product>::const_iterator iter = getIterator(); iter != products.end(); ++iter ) {
		ss << (*iter) << "\n";
	}

	return ss.str();
}

Sum ExpressionInterpreter::parseExpression( string expr ) {
	using namespace boost;

	typedef tokenizer< char_separator<char> > tokenizer;
	char_separator<char> sepLine( ";" );
	char_separator<char> sepTerm( "/" );
	char_separator<char> sepIndex( "," );
	tokenizer tokLine( expr, sepLine );
	Sum loadedExpression = Sum();

	for ( tokenizer::iterator begLine = tokLine.begin(); begLine != tokLine.end(); ++begLine  ) {
		tokenizer tokTerm( *begLine, sepTerm );
		Product nextProduct = Product();
		for ( tokenizer::iterator begTerm = tokTerm.begin(); begTerm != tokTerm.end(); ++begTerm ) {
			if ( (*begTerm).length() > 0 ) {

				if ( (*begTerm)[0] == 'A' ) {
					int order = atoi( (*begTerm).substr( 2, (*begTerm).size() - 1 ).c_str() );  // Extract order n from syntax "A,n".
					nextProduct.addFactor( TermA( order ) );

				} else if ( (*begTerm)[0] == 'F' ) {
					string indexList = (*begTerm).substr( 2, (*begTerm).size() - 1 ); // List of indices as a string.
																					  // Extracted from syntax "D,i,j,k,...".
					vector<string> tokenizedIndexList = vector<string>();
					split( tokenizedIndexList, indexList, is_any_of( "," ) ); // tokenizedIndexList is passed by reference,
																			  // tokens are added to this container.

					vector<int> indices = vector<int>(); // Vector to hold indices in and
														 // pass toTerm D constructor.

					// Cast string tokens in tokenizedIndexList to int and place in indices.
					for ( vector<string>::iterator iter = tokenizedIndexList.begin(); iter < tokenizedIndexList.end(); ++iter ) {
						// TODO: Add exception handler here.
						indices.push_back( atoi( (*iter).c_str() ) );
					}

					nextProduct.addFactor( FourierSum( indices ) );

				} else if ( (*begTerm)[0] == 'D' ) {
					// TODO: Add exception handler here.

					// See comments above for term type 'F' regarding extracting arguments from
					// the string for this term using the boost library method split().
					string args = (*begTerm).substr( 2, (*begTerm).size() - 1 ).c_str();
					vector<string> tokenizedArgs = vector<string>();
					split( tokenizedArgs, args, is_any_of( "," ) );

					// For a TermD object, the syntax should provide no more then three elements;
					// the flavor label and the corresponding spatial indices. If this isn't satisfied,
					// raise an exception. TODO: Raise the exception.
					if ( tokenizedArgs.size() != 3 ) {
						cout << "***ERROR: (ParseException): Bad syntax formatting encountered for factor of type TermD (found argument list size of " << tokenizedArgs.size() << ")." << endl;
					}

					vector<int> indexVector = vector<int>();
					indexVector.push_back( atoi( tokenizedArgs[1].c_str() ) );
					indexVector.push_back( atoi( tokenizedArgs[2].c_str() ) );
					TermD newFactor = TermD( indexVector );
					newFactor.setFlavorLabel( tokenizedArgs[0] );

					nextProduct.addFactor( newFactor );

				} else if ( (*begTerm)[0] == 'C' ) {

					double value = atof( (*begTerm).substr( 2, (*begTerm).size() - 1 ).c_str() );
					nextProduct.addFactor( CoefficientFloat( value ) );
				} else {
					cout << "***ERROR: (ParseExpression): Invalid token encountered: " << (*begTerm) << "." << endl;
				}
			} else {
				cout << "***ERROR: (ParseExpression) Zero-length token encountered." << endl;
			}
		}

		loadedExpression.addProduct( nextProduct );
	}

	return loadedExpression;

}




