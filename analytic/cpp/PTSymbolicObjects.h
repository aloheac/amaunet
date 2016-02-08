/* ***********************************************************************
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * High-order Perturbation Theory Analytics
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Symbolic Peturbation Theory Expression Objects Header
 *
 * v. 0.1		02 Feb 2016
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at CHapel Hill
 * ***********************************************************************
 */

#ifndef PTSYMBOLICOBJECTS_H_
#define PTSYMBOLICOBJECTS_H_

#include <string>
#include <vector>

/*
 * ***********************************************************************
 * FORWARD DECLARATIONS
 * ***********************************************************************
 */

class Sum;
class Product;
class SymbolicTerm;
void unpackTrivialExpression( SymbolicTerm* & );

/*
 * ***********************************************************************
 * EXPRESSION AND TERM BASE CLASSES
 * ***********************************************************************
 */

class SymbolicTerm {
public:

	SymbolicTerm();

	SymbolicTerm( const SymbolicTerm* st );

	virtual ~SymbolicTerm();

	virtual bool operator==( const SymbolicTerm &other ) const;

	virtual Sum getDerivative();

	virtual void simplify();

	virtual void reduceTree();

	virtual const std::string to_string() const;

	void setAsNonInteracting();

	bool isTermInteracting();

	char* getFlavorLabel();

	int* getIndices();

	void setIndices( int* newIndices );

	char getTermID();

	virtual SymbolicTerm* copy();

protected:

	int derivativeOrder;

	bool isInteracting;

	char* flavorLabel;

	int* indices;

	char termID;
};

/*
 * ***********************************************************************
 * TERM DERIVED CLASSES - SCALAR AND MATRIX OBJECTS
 * ***********************************************************************
 */

class GenericTestTerm : public SymbolicTerm {
public:

	GenericTestTerm( int thisId, int thisDerivativeOrder );

	const std::string to_string() const;

	Sum getDerivative();

	GenericTestTerm* copy();

	int id;
};

class MatrixM : public SymbolicTerm {

	friend class MatrixB;

public:

	MatrixM();

	MatrixM( char* flavorLabel );

	MatrixM( const MatrixM* m );

	bool operator==( const MatrixM &other ) const;

	Sum getDerivative();

	const std::string to_string() const;

	MatrixM* copy();
};

class MatrixB : public SymbolicTerm {
public:

	MatrixB();

	MatrixB( char* flavorLabel );

	MatrixB( MatrixB* b );

	bool operator==( const MatrixB &other ) const;

	Sum getDerivative();

	const std::string to_string() const;

	MatrixB* copy();
};

class MatrixK : public SymbolicTerm {
public:

	MatrixK();

	MatrixK( char* flavorLabel );

	MatrixK( const MatrixK* k );

	bool operator==( const MatrixK &other ) const;

	const std::string to_string() const;

	MatrixK* copy();

	void fourierTransform();

private:

	bool isFourierTransformed;
};

class MatrixS : public SymbolicTerm {
public:

	MatrixS();

	MatrixS( MatrixS* s );

	const std::string to_string() const;

	MatrixS* copy();

};

class DetM : public SymbolicTerm {
public:

	DetM( MatrixM matrix );

	DetM( MatrixM matrix, bool inverted );

	DetM( const DetM* m );

	bool isDetInverted();

	const std::string to_string();

	bool operator==( const DetM &other ) const;

	Sum getDerivative();

	DetM* copy();

private:

	bool isInverted;

};

class TermA : public SymbolicTerm {
public:

	TermA();

	TermA( const TermA* A );

	const std::string to_string() const;

	bool operator==( const TermA &other ) const;

	TermA* copy();
};

class CoefficientFloat : public SymbolicTerm {
public:

	CoefficientFloat( double val );

	CoefficientFloat( const CoefficientFloat* f );

	const std::string to_string() const;

	CoefficientFloat* copy();

	Sum getDerivative();

	double eval();

private:

	double value;
};

class CoefficientFraction : public SymbolicTerm {
public:

	CoefficientFraction( int n, int d );

	CoefficientFraction( CoefficientFraction* f );

	const std::string to_string() const;

	CoefficientFraction* copy();

	Sum getDerivative();

	double eval();

private:

	int num, den;
};

/*
 * ***********************************************************************
 * TERM DERIVED CLASSES - EXPRESSION CLASSES
 * ***********************************************************************
 */

class Sum : public SymbolicTerm {

	friend void unpackTrivialExpression( SymbolicTerm* & );

public:

	Sum();

	Sum( std::vector<SymbolicTerm*> thisTerms );

	Sum( SymbolicTerm* term );

	Sum( const Sum* s );

	~Sum();

	const std::string to_string() const;

	Sum* copy();

	void simplify();

	void reduceTree();

	Sum getDerivative();

	Sum* getExpandedExpr();

	void addTerm( SymbolicTerm* thisTerm );

	int getNumberOfTerms();

	void setAsNonInteracting();

	bool operator==( const Sum &other ) const;

	std::vector<SymbolicTerm*>::const_iterator getIteratorBegin() const;

	std::vector<SymbolicTerm*>::const_iterator getIteratorEnd() const;

	std::vector<SymbolicTerm*>::iterator getIteratorBegin();

	std::vector<SymbolicTerm*>::iterator getIteratorEnd();

private:

	std::vector<SymbolicTerm*> terms;
};

class Product : public SymbolicTerm {

	friend void unpackTrivialExpression( SymbolicTerm* & );

public:

	Product();

	Product( std::vector<SymbolicTerm*> t );

	Product( SymbolicTerm* term );

	Product( const Product* p );

	~Product();

	const std::string to_string() const;

	Product* copy();

	void simplify();

	void reduceTree();

	Sum getDerivative();

	Sum* getExpandedExpr();

	void addTerm( SymbolicTerm* t );

	int getNumberOfTerms();

	void setAsNonInteracting();

	bool operator==( const Sum &other ) const;

	bool containsSum();

	std::vector<SymbolicTerm*>::const_iterator getIteratorBegin() const;

	std::vector<SymbolicTerm*>::const_iterator getIteratorEnd() const;

	std::vector<SymbolicTerm*>::iterator getIteratorBegin();

	std::vector<SymbolicTerm*>::iterator getIteratorEnd();

//private:

	std::vector<SymbolicTerm*> terms;

};

class Trace : public SymbolicTerm {
public:

	Trace( SymbolicTerm* expr );

	Trace( const Trace* );

	~Trace();

	const std::string to_string() const;

	void simplify();

	void reduceTree();

	Sum getDerivative();

	void setAsNonInteracting();

	bool operator==( const Trace &other ) const;

	void rewriteInKSFormalism();

private:

	SymbolicTerm* expr;
};

class Delta : public SymbolicTerm {
public:

	Delta( int a, int b );

	Delta( int a, int b, bool type );

	Delta( const Delta* d );

	const std::string to_string() const;

	bool operator==( const Delta &other ) const;

};

class FourierSum : public SymbolicTerm {
public:

	FourierSum( int* i, int Korder );

	FourierSum( const FourierSum * fs );

	~FourierSum();

	const std::string to_string();

	bool operator==( const FourierSum &other ) const;

private:

	int* indices;

	int order;

};

/*
 * ***********************************************************************
 * GENERIC HELPER FUNCTIONS
 * ***********************************************************************
 */

void unpackTrivialExpression( SymbolicTerm*& expr );

bool isZeroTrace( SymbolicTerm* tr );

int getProductAOrder( Product &prod );

int getDualProductAOrder( Product &prodA, Product &prodB );

/*
 * ***********************************************************************
 * EXPRESSION MANIPULATION FUNCTIONS
 * ***********************************************************************
 */

Sum distributeTrace( Trace tr );

Sum distributeAllTraces( Sum expr );

Sum truncateAOrder( Sum expr, int highestOrder );

Sum truncateOddOrders( Sum expr );

Sum rewriteSumInKSFormalism( Sum expr );

Sum indexExpression( Sum expr );

/*
 * ***********************************************************************
 * INPUT REDIRECTION OPERATOR OVERLOADS
 * ***********************************************************************
 */

std::ostream& operator<<( std::ostream& os, const SymbolicTerm &obj );

std::ostream& operator<<( std::ostream& os, const TermA &obj );

std::ostream& operator<<( std::ostream& os, const MatrixM &obj );

#endif /* PTSYMBOLICOBJECTS_H_ */