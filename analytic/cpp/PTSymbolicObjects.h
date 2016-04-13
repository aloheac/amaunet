/* ***********************************************************************
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * High-order Perturbation Theory Analytics
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Symbolic Peturbation Theory Expression Objects Header
 *
 * v. 0.1		09 Feb 2016
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at Chapel Hill
 * ***********************************************************************
 */

#ifndef PTSYMBOLICOBJECTS_H_
#define PTSYMBOLICOBJECTS_H_

#include <string>
#include <vector>
#include <memory>
#include <map>


/*
 * ***********************************************************************
 * FORWARD DECLARATIONS
 * ***********************************************************************
 */

class Sum;
class Product;
class Trace;
class SymbolicTerm;
class GenericTestTerm;
class MatrixM;
class MatrixB;
class MatrixK;
class MatrixS;
class DetM;
class TermA;
class CoefficientFloat;
class CoefficientFraction;
class Delta;
class FourierSum;
class IndexContraction;
class DeltaContractionSet;
bool unpackTrivialExpression( std::shared_ptr<SymbolicTerm> & );

/*
 * ***********************************************************************
 * TYPEDEF DECLARATIONS
 * ***********************************************************************
 */

typedef std::shared_ptr<Sum> SumPtr;
typedef std::shared_ptr<Product> ProductPtr;
typedef std::shared_ptr<Trace> TracePtr;
typedef std::shared_ptr<SymbolicTerm> SymbolicTermPtr;
typedef std::shared_ptr<GenericTestTerm> GenericTestTermPtr;
typedef std::shared_ptr<MatrixM> MatrixMPtr;
typedef std::shared_ptr<MatrixB> MatrixBPtr;
typedef std::shared_ptr<MatrixK> MatrixKPtr;
typedef std::shared_ptr<MatrixS> MatrixSPtr;
typedef std::shared_ptr<DetM> DetMPtr;
typedef std::shared_ptr<TermA> TermAPtr;
typedef std::shared_ptr<Delta> DeltaPtr;
typedef std::shared_ptr<CoefficientFloat> CoefficientFloatPtr;
typedef std::shared_ptr<CoefficientFraction> CoefficientFractionPtr;
typedef std::shared_ptr<FourierSum> FourierSumPtr;

/*
 * ***********************************************************************
 * ENUMERATED TYPES
 * ***********************************************************************
 */

enum class TermTypes : char {
	INVALID_TERM = '0',
	GENERIC_TEST_TERM = 'g',
	MATRIX_M = 'M',
	MATRIX_B = 'B',
	MATRIX_K = 'K',
	MATRIX_S = 's',
	DET_M = 'D',
	TERM_A = 'A',
	COEFFICIENT_FLOAT = 'L',
	COEFFICIENT_FRACTION = 'R',
	SUM = 'S',
	PRODUCT = 'P',
	TRACE = 'T',
	DELTA = 'd',
	FOURIER_SUM = 'F'
};

/*
 * ***********************************************************************
 * EXPRESSION AND TERM BASE CLASSES
 * ***********************************************************************
 */

class SymbolicTerm {

	friend bool unpackTrivialExpression( SymbolicTermPtr & );

public:

	SymbolicTerm();

	virtual ~SymbolicTerm();

	virtual bool operator==( const SymbolicTerm &other ) const;

	virtual Sum getDerivative();

	virtual void simplify();

	virtual void reduceTree();

	virtual const std::string to_string() const;

	virtual void setAsNonInteracting();

	bool isTermInteracting();

	const char* getFlavorLabel();

	int* getIndices();

	void setIndices( int i, int j );

	TermTypes getTermID();

	virtual SymbolicTermPtr copy();

protected:

	int derivativeOrder;

	bool isInteracting;

	const char* flavorLabel;

	int indices[ 2 ];

	TermTypes termID;
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

	SymbolicTermPtr copy();

	int id;
};

class MatrixM : public SymbolicTerm {

	friend class MatrixB;

public:

	MatrixM();

	MatrixM( const char* flavorLabel );

	bool operator==( const MatrixM &other ) const;

	Sum getDerivative();

	const std::string to_string() const;

	SymbolicTermPtr copy();
};

class MatrixB : public SymbolicTerm {
public:

	MatrixB();

	MatrixB( const char* flavorLabel );

	bool operator==( const MatrixB &other ) const;

	Sum getDerivative();

	const std::string to_string() const;

	SymbolicTermPtr copy();
};

class MatrixK : public SymbolicTerm {
public:

	MatrixK();

	MatrixK( const char* flavorLabel );

	bool operator==( const MatrixK &other ) const;

	const std::string to_string() const;

	SymbolicTermPtr copy();

	void fourierTransform();

private:

	bool isFourierTransformed;
};

class MatrixS : public SymbolicTerm {
public:

	MatrixS();

	const std::string to_string() const;

	SymbolicTermPtr copy();

};

class DetM : public SymbolicTerm {
public:

	DetM( MatrixM matrix );

	DetM( MatrixM matrix, bool inverted );

	bool isDetInverted();

	void setAsNonInteracting();

	const std::string to_string() const;

	bool operator==( const DetM &other ) const;

	Sum getDerivative();

	SymbolicTermPtr copy();

private:

	bool isInverted;

	MatrixM matrix;

};

class TermA : public SymbolicTerm {
public:

	TermA();

	TermA( const TermA* A );

	const std::string to_string() const;

	bool operator==( const TermA &other ) const;

	SymbolicTermPtr copy();
};

class CoefficientFraction : public SymbolicTerm {
public:

	CoefficientFraction();

	CoefficientFraction( double n, double d );

	CoefficientFraction( CoefficientFraction* f );

	const std::string to_string() const;

	CoefficientFraction operator*( const CoefficientFraction& obj ) const;

	CoefficientFraction operator+( const CoefficientFraction& obj ) const;

	CoefficientFraction operator*( const CoefficientFloat& obj ) const;

	CoefficientFraction operator+( const CoefficientFloat& obj ) const;

	CoefficientFraction& operator*=( const CoefficientFraction& obj );

	CoefficientFraction& operator+=( const CoefficientFraction& obj );

	SymbolicTermPtr copy();

	Sum getDerivative();

	double eval();

	void reduce();

private:

	double num, den;
};

class CoefficientFloat : public SymbolicTerm {

	friend class CoefficientFraction; // TODO: Replace friend declaration with making eval() const.

public:

	CoefficientFloat( double val );

	CoefficientFloat( const CoefficientFloat* f );

	const std::string to_string() const;

	CoefficientFloat operator*( const CoefficientFloat& obj ) const;

	CoefficientFloat operator+( const CoefficientFloat& obj ) const;

	CoefficientFraction operator*( const CoefficientFraction& obj ) const;

	CoefficientFraction operator+( const CoefficientFraction& obj ) const;

	CoefficientFloat& operator*=( const CoefficientFloat& obj );

	CoefficientFloat& operator+=( const CoefficientFloat& obj );

	SymbolicTermPtr copy();

	Sum getDerivative();

	double eval();

private:

	double value;
};

/*
 * ***********************************************************************
 * TERM DERIVED CLASSES - EXPRESSION CLASSES
 * ***********************************************************************
 */

class Sum : public SymbolicTerm {

	friend bool unpackTrivialExpression( SymbolicTermPtr & );

	friend Sum truncateOddOrders( SymbolicTermPtr expr );

	friend class Product;

public:

	Sum();

	Sum( std::vector<SymbolicTermPtr> thisTerms );

	Sum( SymbolicTermPtr term );

	Sum( const Sum &s );

	~Sum();

	Sum& operator=( const Sum &rhs );

	const std::string to_string() const;

	SymbolicTermPtr copy();

	void simplify();

	void reduceTree();

	Sum getDerivative();

	Sum getExpandedExpr();

	void addTerm( SymbolicTermPtr thisTerm );

	int getNumberOfTerms();

	void setAsNonInteracting();

	bool operator==( const Sum &other ) const;

	void reduceFourierSumIndices();

	void clear();

	std::vector<SymbolicTermPtr>::const_iterator getIteratorBegin() const;

	std::vector<SymbolicTermPtr>::const_iterator getIteratorEnd() const;

	std::vector<SymbolicTermPtr>::iterator getIteratorBegin();

	std::vector<SymbolicTermPtr>::iterator getIteratorEnd();

private:

	std::vector<SymbolicTermPtr> terms;
};

class Product : public SymbolicTerm {

	friend bool unpackTrivialExpression( SymbolicTermPtr & );

	friend class Sum;

	friend void indexExpression( SymbolicTermPtr expr );

public:

	Product();

	Product( std::vector<SymbolicTermPtr> t );

	Product( SymbolicTermPtr term );

	~Product();

	Product& operator=( const Product &rhs );

	const std::string to_string() const;

	SymbolicTermPtr copy();

	void simplify();

	void reduceTree();

	Sum getDerivative();

	Sum getExpandedExpr();

	void addTerm( SymbolicTermPtr t );

	int getNumberOfTerms();

	void setAsNonInteracting();

	bool operator==( const Sum &other ) const;

	bool containsSum();

	void zero();

	void clear();

	void reduceFourierSumIndices();

	std::vector<SymbolicTermPtr>::const_iterator getIteratorBegin() const;

	std::vector<SymbolicTermPtr>::const_iterator getIteratorEnd() const;

	std::vector<SymbolicTermPtr>::iterator getIteratorBegin();

	std::vector<SymbolicTermPtr>::iterator getIteratorEnd();

private:

	std::vector<SymbolicTermPtr> terms;

};

class Trace : public SymbolicTerm {

	friend bool isZeroTrace( SymbolicTermPtr );

	friend Sum distributeTrace( SymbolicTermPtr tr );

	friend void indexExpression( SymbolicTermPtr expr );

public:

	Trace( SymbolicTermPtr expr );

	~Trace();

	const std::string to_string() const;

	SymbolicTermPtr copy();

	void simplify();

	void reduceTree();

	Sum getDerivative();

	void setAsNonInteracting();

	bool operator==( const Trace &other ) const;

	void rewriteInKSFormalism();

private:

	SymbolicTermPtr expr;
};

class Delta : public SymbolicTerm {
public:

	Delta( int a, int b );

	Delta( int a, int b, bool type );

	const std::string to_string() const;

	SymbolicTermPtr copy();

	bool operator==( const Delta &other ) const;

	bool isDeltaBar();

private:

	bool isBar;

};

class FourierSum : public SymbolicTerm {
public:

	FourierSum( std::vector<IndexContraction> i, int orderInK );

	~FourierSum();

	const std::string to_string() const;

	SymbolicTermPtr copy();

	bool operator==( const FourierSum &other ) const;

	void reduceDummyIndices();

private:

	std::vector<IndexContraction> indices;

	int order;

};

/*
 * ***********************************************************************
 * GENERIC HELPER FUNCTIONS
 * ***********************************************************************
 */

bool unpackTrivialExpression( SymbolicTermPtr& expr );

bool isZeroTrace( SymbolicTermPtr tr );

int getProductAOrder( SymbolicTermPtr prod );

int getDualProductAOrder( Product &prodA, Product &prodB );

int getTerminatedContraction( std::map<int, int> contractedIndexMapping, int index );

std::map<int,int> constructContractionDictionary( DeltaContractionSet contractions );

bool areTermsCommon( SymbolicTermPtr termA, SymbolicTermPtr termB );

int gcd( int a, int b );

/*
 * ***********************************************************************
 * EXPRESSION MANIPULATION FUNCTIONS
 * ***********************************************************************
 */

Sum distributeTrace( SymbolicTermPtr tr );

Sum distributeAllTraces( SymbolicTermPtr expr );

Sum truncateAOrder( SymbolicTermPtr expr, int highestOrder );

Sum truncateOddOrders( SymbolicTermPtr expr );

void rewriteSumInKSFormalism( SymbolicTermPtr expr );

void indexExpression( SymbolicTermPtr expr );

Sum fourierTransformExpression( SymbolicTermPtr expr );

Sum combineLikeTerms( Sum &expr );

Sum combineLikeTerms( Sum &expr, int groupSize );

/*
 * ***********************************************************************
 * INPUT REDIRECTION OPERATOR OVERLOADS
 * ***********************************************************************
 */

std::ostream& operator<<( std::ostream& os, const SymbolicTerm &obj );

std::ostream& operator<<( std::ostream& os, const TermA &obj );

std::ostream& operator<<( std::ostream& os, const MatrixM &obj );

std::ostream& operator<<( std::ostream& os, const TermTypes &obj );

std::ostream& operator<<( std::ostream& os, const std::map<int, int> &obj );

#endif /* PTSYMBOLICOBJECTS_H_ */
