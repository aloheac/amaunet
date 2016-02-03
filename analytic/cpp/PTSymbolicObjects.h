/* ***********************************************************************
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * Numerical CUDA Perturbation Theory Fourier Transform Evaluation
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Symbolic Peturbation Theory Expression Objects Header
 *
 * v. 0.1		27 Jan 2015
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

	virtual Sum* getDerivative();

	virtual void simplify();

	virtual void reduceTree();

	virtual const std::string to_string();

	void setAsNonInteracting();

	bool isTermInteracting();

	char* getFlavorLabel();

	int* getIndices();

	void setIndices( int* newIndices );

	char getTermID();

	SymbolicTerm* copy();

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

class MatrixM : public SymbolicTerm {

	friend class MatrixB;

public:

	MatrixM();

	MatrixM( char* flavorLabel );

	MatrixM( const MatrixM* m );

	bool operator==( const MatrixM &other ) const;

	Sum* getDerivative();

	const std::string to_string();
};

class MatrixB : public SymbolicTerm {
public:

	MatrixB();

	MatrixB( char* flavorLabel );

	MatrixB( MatrixB* b );

	bool operator==( const MatrixB &other ) const;

	Sum* getDerivative();

	const std::string to_string();
};

class MatrixK : public SymbolicTerm {
public:

	MatrixK();

	MatrixK( char* flavorLabel );

	MatrixK( const MatrixK* k );

	bool operator==( const MatrixK &other ) const;

	const std::string to_string();

	void fourierTransform();

private:

	bool isFourierTransformed;
};

class MatrixS : public SymbolicTerm {
public:

	MatrixS();

	MatrixS( MatrixS* s );

	const std::string to_string();

};

class DetM : public SymbolicTerm {
public:

	DetM( MatrixM matrix );

	DetM( MatrixM matrix, bool inverted );

	DetM( const DetM* m );

	bool isDetInverted();

	const std::string to_string();

	bool operator==( const DetM &other ) const;

	Sum* getDerivative();

private:

	bool isInverted;

};

class TermA : public SymbolicTerm {
public:

	TermA();

	TermA( const TermA* A );

	const std::string to_string();

	bool operator==( const TermA &other ) const;
};

class CoefficientFloat : public SymbolicTerm {
public:

	CoefficientFloat( double val );

	CoefficientFloat( const CoefficientFloat* f );

	const std::string to_string();

	Sum* getDerivative();

	double eval();

private:

	double value;
};

class CoefficientFraction : public SymbolicTerm {
public:

	CoefficientFraction( int n, int d );

	CoefficientFraction( CoefficientFraction* f );

	const std::string to_string();

	Sum* getDerivative();

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
public:

	Sum();

	Sum( std::vector<SymbolicTerm*> thisTerms );

	Sum( const Sum* s );

	~Sum();

	const std::string to_string();

	void simplify();

	void reduceTree();

	Sum* getDerivative();

	Sum* getExpandedExpr();

	void addTerm( SymbolicTerm* thisTerm );

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
public:

	Product();

	Product( std::vector<SymbolicTerm*> t );

	Product( const Product* p );

	~Product();

	const std::string to_string();

	void simplify();

	void reduceTree();

	Sum* getDerivative();

	Sum* getExpandedExpr();

	void addTerm( SymbolicTerm* t );

	void setAsNonInteracting();

	bool operator==( const Sum &other ) const;

	bool containsSum();

	std::vector<SymbolicTerm*>::const_iterator getIteratorBegin() const;

	std::vector<SymbolicTerm*>::const_iterator getIteratorEnd() const;

	std::vector<SymbolicTerm*>::iterator getIteratorBegin();

	std::vector<SymbolicTerm*>::iterator getIteratorEnd();

private:

	std::vector<SymbolicTerm*> terms;

};

class Trace : public SymbolicTerm {
public:

	Trace( SymbolicTerm* expr );

	Trace( const Trace* );

	~Trace();

	const std::string to_string();

	void simplify();

	void reduceTree();

	Sum* getDerivative();

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

	const std::string to_string();

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

void unpackTrivialExpression( SymbolicTerm* expr );

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

#endif /* PTSYMBOLICOBJECTS_H_ */
