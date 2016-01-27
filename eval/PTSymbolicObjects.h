/*
 * PTSymbolicObjects.h
 *
 *  Created on: Jan 22, 2016
 *      Author: loheac
 */

#ifndef PTSYMBOLICOBJECTS_H_
#define PTSYMBOLICOBJECTS_H_

#include <string>
#include <vector>
#include <ostream>
#include "PTEvalUtils.h"

/*
 * ****************************************************************************
 * BASE AND DERIVED CLASS DEFINITIONS FOR EXPRESSON OBJECTS
 * ****************************************************************************
 */

// Forward declarations:
class Factor;

class Expression {
public:
	Expression();

	/*
	 * Typical-use constructor.
	 * @param numTerms Number of terms in this product.
	 * @param terms Array of factors in this product.
	 */
	Expression( int numTerms, Factor* terms );

private:
	/*
	 * Number of terms in this expression (length of array terms).
	 */
	int numTerms;

	/*
	 * Array of terms in this expression.
	 */
	Factor* terms;

};

/*
 * Representation of a product of factors. Expressions are already in a fully
 * distributed form ("standard representation"), so therefore this
 * implementation is very basic and used simply as a container.
 */
class Product {
public:
	Product();

	void addFactor( Factor f );

	std::vector<Factor>::const_iterator getIterator() const;

	friend std::ostream& operator<<( std::ostream& os, const Product& prod );

	std::string to_string() const;

protected:
	std::vector<Factor> factors;

};

class Sum {
public:
	Sum();

	void addProduct( Product p );

	std::vector<Product>::const_iterator getIterator() const;

	friend std::ostream& operator<<( std::ostream& os, const Sum& s );

	std::string to_string() const;

protected:
	std::vector<Product> products;
};

/*
 * ****************************************************************************
 * BASE AND DERIVED CLASS DEFINITIONS FOR SCALAR AND MATRIX OBJECTS
 * ****************************************************************************
 */

/*
 * Base class for the representation of a factor in a product. In the context
 * of the current perturbation theory, derived classes can be representations
 * for A, D, or a generalized Fourier sum.
 */
class Factor {

public:
	Factor();

	Factor( std::vector<int> indices );

	/*
	 * Numerically evaluate the factor for a given set of indices. The sum
	 * over free indices is implemented externally; numerical values of
	 * indices are taken to be explicit. The base class contains no
	 * mathematically valid implementation. All derived classes must
	 * override with their own appropriate definition.
	 * @param numIndices The number of indices passed to the method.
	 * @param indices Array of integers that represent the explicit indices.
	 */
	double eval( std::vector<int> indices, pt_util::PTSystemParameters params );

	std::string to_string() const;

	friend std::ostream& operator<<( std::ostream& os, const Factor& f );

	void setFlavorLabel( std::string label );

	std::string getFlavorLabel();

protected:
	/*
	 * Number of indices this term carries (length of array indices).
	 */
	int numIndices;

	/*
	 * Array of indices that this term carries. Note that depending of the
	 * stage of the calculation that this object is used in, these indices
	 * may be free or explicit. This indices are not used when eval() is
	 * called, since at that point all indices must be explicit.
	 */
	std::vector<int> indices;

	std::string stringRepresentation;

	std::string flavorLabel;
};

/*
 * Representation of the scalar A in terms of weak-coupling perturbation
 * theory. Unlike the analytic implementation, an order n is assigned
 * to a single instance to represent A^n.
 */
class TermA : public Factor {
public:
	/*
	 * Constructor for TermA.
	 * @param order The order of the scalar A.
	 */
	TermA( int order );

	std::string to_string() const;

	double eval( std::vector<int> indices, pt_util::PTSystemParameters params );

private:
	/*
	 * Mathematical order (power) of A, to represent A^{order}.
	 */
	int order;
};

/*
 * Representation of the diagonal matrix D in terms of weak-coupling
 * perturbation theory.
 */
class TermD : public Factor {
public:
	/*
	 * Constructor for TermD.
	 * @param indices Array of indices that this term carries. Depending on the
	 * 				  stage of the calculation, these indices may be free or
	 * 				  explicit.
	 */
	TermD( std::vector<int> indices );  // Recall that D always carries a single index in
							// momentum space, so numIndices = 1.

	std::string to_string() const;

	double eval( std::vector<int> indices, pt_util::PTSystemParameters params );
};

class FourierSum : public Factor {
public:
	FourierSum( std::vector<int> indices );

	std::string to_string() const;

	double eval( std::vector<int> indices, pt_util::PTSystemParameters params );
};

class CoefficientFloat : public Factor {
public:
	CoefficientFloat( double value );

	std::string to_string() const;

	double eval( std::vector<int> indices, pt_util::PTSystemParameters params );

private:
	double value;
};

class ExpressionInterpreter {
public:
	ExpressionInterpreter();

	Sum parseExpression( std::string ss );

};

#endif /* PTSYMBOLICOBJECTS_H_ */
