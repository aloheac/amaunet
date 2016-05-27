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
#include <set>

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
class MatrixK;
class MatrixS;
class TermA;
class TermE;
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
typedef std::shared_ptr<MatrixK> MatrixKPtr;
typedef std::shared_ptr<MatrixS> MatrixSPtr;
typedef std::shared_ptr<TermA> TermAPtr;
typedef std::shared_ptr<TermE> TermEPtr;
typedef std::shared_ptr<Delta> DeltaPtr;
typedef std::shared_ptr<CoefficientFloat> CoefficientFloatPtr;
typedef std::shared_ptr<CoefficientFraction> CoefficientFractionPtr;
typedef std::shared_ptr<FourierSum> FourierSumPtr;

/*
 * ***********************************************************************
 * ENUMERATED TYPES
 * ***********************************************************************
 */

/**
 * Enumerations used to uniquely define each type of symbolic term that may appear in an expression. An instance of
 * TermTypes is most typically used in the termID field of the SymbolicTerm class. Classes that are derived from
 * SymbolicTerm indicates which subclass the particular instance is so that static casts to the derived class type
 * can be made after checking which TermTypes enum is held in termID.
 *
 * Note that if INVALID_TERM ever appears during normal execution, an error has very likely occured.
 */
enum class TermTypes : char {
	INVALID_TERM = '0',
	GENERIC_TEST_TERM = 'g',
	MATRIX_K = 'K',
	MATRIX_S = 's',
	TERM_A = 'A',
	TERM_E = 'E',
	COEFFICIENT_FLOAT = 'L',
	COEFFICIENT_FRACTION = 'R',
	SUM = 'S',
	PRODUCT = 'P',
	TRACE = 'T',
	DELTA = 'd',
	FOURIER_SUM = 'F',
	DEBUG_TRACE = 'Z'
};

/*
 * ***********************************************************************
 * EXPRESSION AND TERM BASE CLASSES
 * ***********************************************************************
 */

/**
 * Base class for all expressions and symbolic objects that may appear in an expression, given the current formalism of
 * perturbation theory.
 */
class SymbolicTerm {

	friend bool unpackTrivialExpression( SymbolicTermPtr & );

public:

	/**
	 * SymbolicTerm default constructor. Accepts no arguments.
	 */
	SymbolicTerm();

	/**
	 * SymbolicTerm destructor.
	 */
	virtual ~SymbolicTerm();

	/**
	 * Equality operator overload. This method should always be overridden by derived classes.
	 * @param other SymbolicTerm instance which is being compared.
	 * @return true if other is mathematically and symbolically equivalent to this instance, false otherwise.
	 */
	virtual bool operator==( const SymbolicTerm &other ) const;

	/**
	 * Mathematically simplifies the expression in a trivial manner. By the present implementation, all Product
	 * expressions that contain an identifiable expression of zero itself is reduced to zero, and factors of one
	 * are removed from the Product. Any terms in a Sum which are zero are removed. All implementations are provided
	 * by derived classes.
	 */
	virtual void simplify();

	/**
	 * Reduces the expression tree of the data structure to canonical form. All expression manipulation functions
	 * typically expect passed expressions to be in canonical form.
	 */
	virtual void reduceTree();

	/**
	 * Returns a pretty-printed representation of the expression. All derived classes should override this method.
	 * @return A pretty-printed representation of the expression.
	 */
	virtual const std::string to_string() const;

	/**
	 * Gets the particle flavor label assigned to this expression. Flavor labels are not appropriate for all types of
	 * expressions, and is only used in cases where matrices are dependent on the particle species it refers to.
	 * @return Flavor label assigned to this expression.
	 */
	const char* getFlavorLabel();

	/**
	 * Gets the indices assigned to this expression. Two and only two indices may be assigned to an expression under
	 * the present formalism. Exactly two integers are allocated for this array for each instance. The indices are
	 * commonly referred to as i and j.
	 * @return Two element array of indices assigned to this expression.
	 */
	int* getIndices();

	/**
	 * Sets the two indices assigned to this expression. Indices are held in the two-element int array field indices.
	 * @param i The first index.
	 * @param j The second index.
	 */
	void setIndices( int i, int j );

	/**
	 * Gets the term type identifier, which is a valid instance of the TermTypes enumeration. All derived classes set
	 * the term identifier associated with the derived class in constructors.
	 * @return The term type identifier assigned to this symbolic object.
	 */
	TermTypes getTermID();

	/**
	 * Generates a deep copy of this SymbolicTerm instance on the heap and returns a shared smart pointer to the copy.
	 * All derived classes should implement a method that overrides the base class definition. If the base class method
	 * is called directly, the symbolic object will have the term identifier of an invalid term.
	 * @return A shared smart pointer to the copy on the heap.
	 */
	virtual SymbolicTermPtr copy();

protected:

	/**
	 * The particle flavor label (e.g. "up", "dn") associated with this symbolic object, if applicable.
	 */
	const char* flavorLabel;

	/**
	 * The coordinate- or momentum-space indices associated with this symbolic object, if applicable. Under the current
	 * formalism any object which carries indices, including matrices and delta functions, should carry exactly two
	 * indices. Abstract indices are generally referred to with integers in this formalism.
	 */
	int indices[ 2 ];

	/**
	 * The term type identifier associated with this symbolic object.
	 */
	TermTypes termID;
};

/*
 * ***********************************************************************
 * TERM DERIVED CLASSES - SCALAR AND MATRIX OBJECTS
 * ***********************************************************************
 */

/**
 * Symbolic representation of the matrix K in terms of the current perturbation theory formalism. Class is derived from
 * the base class SymbolicTerm. This matrix is defined as
 *
 * 	K \equiv M_0^{-1} \mathcal{T}
 *
 * 	where M_0 indicates the matrix M in the non-interacting limit (A = 0) and \mathcal{T} is the matrix of kinetic
 * 	energy operators. This matrix is later Fourier transformed to a diagonal representation in momentum space. The
 * 	pretty-printed representation of this symbolic object is
 *
 * 	K_##_( a, b )
 *
 * 	where ## is the flavor label assigned to this matrix, and ( a, b ) are the corresponding coordinate-space indices
 * 	of this element, after trace operators have been removed after indexing. Note that the index set ( 0, 0 )
 * 	corresponds to an invalid set of indices where no indices have yet been assigned.
 */
class MatrixK : public SymbolicTerm {

public:

	/**
	 * The default constructor for MatrixK. Assigns the empty string for the flavor label.
	 */
	MatrixK();

	/**
	 * Constructor which assigns a given flavor label to the matrix.
	 * @param flavorLabel The flavor label to assign.
	 */
	MatrixK( const char* flavorLabel );

	/**
	 * Equality operator overload. Two MatrixK objects are considered to be equivalent if the flavor labels of both
	 * objects are equal and if both objects are or are not symbolically Fourier transformed. The indices of both
	 * objects do not need to be equal since the objects commute after the expression in indexed. This operator is
	 * generally used only when combining like terms, and becomes important when considering polarized systems.
	 * @param other The other instance of MatrixK to check for mathematical equivalance.
	 * @return true if this instance and other are mathematically equivalent, false otherwise.
	 */
	bool operator==( const MatrixK &other ) const;

	/**
	 * Gets the pretty-printed representation of this instance of MatrixK. See the class description for clarification
	 * on this representation.
	 * @return The pretty-printed representation of this instance of MatrixK.
	 */
	const std::string to_string() const;

	/**
	 * Generates a deep copy of this instance on the heap and returns a smart shared pointer to the copy.
	 * @return A smart shared pointer to the generated copy.
	 */
	SymbolicTermPtr copy();

	/**
	 * Marks this instance of MatrixK as Fourier transformed to momentum space. The result is refered to as the
	 * propagator D in the formalism.
	 */
	void fourierTransform();

private:

	/**
	 * Boolean flag indicating whether this instance is symbolically Fourier transformed to momentum space. This field
	 * is always initialized to false.
	 */
	bool isFourierTransformed;
};

class MatrixS : public SymbolicTerm {
public:

	MatrixS();

	const std::string to_string() const;

	SymbolicTermPtr copy();

};

class TermA : public SymbolicTerm {
public:

	TermA();

	TermA( const TermA* A );

	const std::string to_string() const;

	bool operator==( const TermA &other ) const;

	SymbolicTermPtr copy();
};

class TermE : public SymbolicTerm {

public:

	TermE( int thisOrder );

	TermE( int thisOrder, const char* thisFlavorLabel );

	const std::string to_string() const;

	bool operator==( const TermE &other ) const;

	SymbolicTermPtr copy();

	int getOrder();

	SymbolicTermPtr getFullExpression();

private:

	int order;

};

class CoefficientFraction : public SymbolicTerm {
public:

	CoefficientFraction();

	CoefficientFraction( double n, double d );

	const std::string to_string() const;

	CoefficientFraction operator*( const CoefficientFraction& obj ) const;

	CoefficientFraction operator+( const CoefficientFraction& obj ) const;

	CoefficientFraction operator*( const CoefficientFloat& obj ) const;

	CoefficientFraction operator+( const CoefficientFloat& obj ) const;

	CoefficientFraction& operator*=( const CoefficientFraction& obj );

	CoefficientFraction& operator+=( const CoefficientFraction& obj );

	SymbolicTermPtr copy();

	double eval();

	void reduce();

private:

	double num, den;
};

class CoefficientFloat : public SymbolicTerm {

	friend class CoefficientFraction; // TODO: Replace friend declaration with making eval() const.

public:

	CoefficientFloat( double val );

	const std::string to_string() const;

	CoefficientFloat operator*( const CoefficientFloat& obj ) const;

	CoefficientFloat operator+( const CoefficientFloat& obj ) const;

	CoefficientFraction operator*( const CoefficientFraction& obj ) const;

	CoefficientFraction operator+( const CoefficientFraction& obj ) const;

	CoefficientFloat& operator*=( const CoefficientFloat& obj );

	CoefficientFloat& operator+=( const CoefficientFloat& obj );

	SymbolicTermPtr copy();

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

	Sum getExpandedExpr();

	void addTerm( SymbolicTermPtr thisTerm );

	int getNumberOfTerms();

	bool operator==( const Sum &other ) const;

	void reduceFourierSumIndices();

	void clear();
        
        void combineCoefficients();

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

	Sum getExpandedExpr();

	void addTerm( SymbolicTermPtr t );

	int getNumberOfTerms();

	bool operator==( const Sum &other ) const;

	bool containsSum();

	void zero();

	void clear();

	void reduceFourierSumIndices();
        
        void combineCoefficients();

	std::vector<SymbolicTermPtr>::const_iterator getIteratorBegin() const;

	std::vector<SymbolicTermPtr>::const_iterator getIteratorEnd() const;

	std::vector<SymbolicTermPtr>::iterator getIteratorBegin();

	std::vector<SymbolicTermPtr>::iterator getIteratorEnd();

private:

	std::vector<SymbolicTermPtr> terms;

};

class Trace : public SymbolicTerm {

	friend bool isZeroTrace( SymbolicTermPtr );

	friend void indexExpression( SymbolicTermPtr expr );

public:

	Trace( SymbolicTermPtr expr );

	~Trace();

	const std::string to_string() const;

	SymbolicTermPtr copy();

	void simplify();

	void reduceTree();

	bool operator==( const Trace &other ) const;
        
        bool operator<( const Trace &other ) const;

	bool operator>( const Trace &other ) const;

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

std::vector< std::set<int> > groupContractions( std::vector< std::set<int> > groups );

std::map<int,int> constructContractionDictionary( DeltaContractionSet contractions );

bool areTermsCommon( SymbolicTermPtr termA, SymbolicTermPtr termB );

unsigned long gcd( unsigned long a, unsigned long b );

int factorial( int n );

/*
 * ***********************************************************************
 * EXPRESSION MANIPULATION FUNCTIONS
 * ***********************************************************************
 */

Sum truncateAOrder( SymbolicTermPtr expr, int highestOrder );

Sum truncateOddOrders( SymbolicTermPtr expr );

void indexExpression( SymbolicTermPtr expr );

Sum fourierTransformExpression( SymbolicTermPtr expr );

Sum combineLikeTerms( Sum &expr );

Sum combineLikeTerms( Sum &expr, int groupSize );

Sum generateExponentialSeries( int order, Product x );

Sum generateDeterminantExpansion( int order, const char* flavorLabel, bool insertFullE );

Sum sortTracesByOrder( Sum &expr );

/*
 * ***********************************************************************
 * INPUT REDIRECTION OPERATOR OVERLOADS
 * ***********************************************************************
 */

std::ostream& operator<<( std::ostream& os, const SymbolicTerm &obj );

std::ostream& operator<<( std::ostream& os, const TermA &obj );

std::ostream& operator<<( std::ostream& os, const TermTypes &obj );

std::ostream& operator<<( std::ostream& os, const std::map<int, int> &obj );

#endif /* PTSYMBOLICOBJECTS_H_ */
