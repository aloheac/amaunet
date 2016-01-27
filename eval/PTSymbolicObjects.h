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

/*
 * Forward declarations.
 */

class Factor;

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

	std::vector<Factor>::const_iterator getEnd() const;

	friend std::ostream& operator<<( std::ostream& os, const Product& prod );

	std::string to_string() const;

	int getNumberOfFactors() const;

	int getOrderInA() const;

	int getNumberOfUniqueIndices() const;

	void finalize();

//protected:
	std::vector<Factor> factors;

	int orderInA;

	int numUniqueIndices;

	bool isFinalized;

	void calcOrderInA();

	void calcNumberOfUniqueIndices();

};

class Sum {
public:
	Sum();

	void addProduct( Product p );

	std::vector<Product>::const_iterator getIterator() const;

	std::vector<Product>::const_iterator getEnd() const;

	friend std::ostream& operator<<( std::ostream& os, const Sum& s );

	std::string to_string() const;

	int getNumberOfTerms() const;

//protected:
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

	virtual ~Factor() { }

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
	virtual double eval( std::vector<int> indices, pt_util::PTSystemParameters params );

	virtual std::string to_string() const;

	friend std::ostream& operator<<( std::ostream& os, const Factor& f );

	void setFlavorLabel( std::string label );

	std::string getFlavorLabel();

	char getFactorType() const;

	std::vector<int>::const_iterator getIndexIterator() const;

	std::vector<int>::const_iterator getIndexEnd() const;

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

	char factorType;
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

	//virtual ~TermA() { }

	std::string to_string() const;

	double eval( std::vector<int> indices, pt_util::PTSystemParameters params );

	int getOrder() const;

//private:
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

	void setFlavorLabel( std::string label );  // Override implementation in Factor.

	double eval( std::vector<int> indices, pt_util::PTSystemParameters params );
};

class DeltaBar : public Factor {
public:
	DeltaBar( std::vector<int> indices );

	std::string to_string() const;

	friend std::ostream& operator<<( std::ostream& os, const DeltaBar& db );

};

class FourierSum : public Factor {
public:
	FourierSum( std::vector<int> indices );

	std::string to_string() const;

	friend std::ostream& operator<<( std::ostream& os, const FourierSum& fs );

	double eval( std::vector<int> indices, pt_util::PTSystemParameters params );
};

class CoefficientFloat : public Factor {
public:
	CoefficientFloat( double value );

	std::string to_string() const;

	friend std::ostream& operator<<( std::ostream& os, const CoefficientFloat& cf );

	double eval( std::vector<int> indices, pt_util::PTSystemParameters params );

	void setValue( double v );

	double getValue();

private:
	double value;
};

class ExpressionInterpreter {
public:
	ExpressionInterpreter();

	Sum parseExpression( std::string ss );

};

#endif /* PTSYMBOLICOBJECTS_H_ */
