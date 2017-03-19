/* ***********************************************************************
 * Amaunet: High-order Lattice Perturbation Theory
 *          for Non-Relativistic Quantum Matter
 *
 * High-order Perturbation Theory Analytics
 * Weak-coupling Expansion for Fermionic Contact Interactions
 *
 * Symbolic Perturbation Theory Expression Objects Header
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
	 * Copy constructor for SymbolicTerm. Simply performs a member-wise copy of the passed object, but is provided for
	 * explicit clarity. Derived classes should implement their own copy constructor where appropriate.
	 * @param other Instance of SymbolicTerm to copy.
	 */
	SymbolicTerm( const SymbolicTerm &other );

	/**
	 * SymbolicTerm destructor.
	 */
	virtual ~SymbolicTerm();

	/**
	 * Assignment operator for SymbolicTerm. Derived classes should override and implement their own assignment
	 * operators.
	 * @param rhs Instance of SymbolicTerm to copy.
	 */
	virtual SymbolicTerm& operator=( const SymbolicTerm& rhs );

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
	 * Copy constructor for MatrixK.
	 * @param other Instance of MatrixK to copy.
	 */
	MatrixK( const MatrixK &other );

	/**
	 * Assignment operator overload for MatrixK.
	 * @param rhs Instance of MatrixK to copy.
	 */
	MatrixK& operator=( const MatrixK &rhs );

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

/**
 * Symbolic representation of the matrix S in terms of the current perturbation theory formalism. Class is derived from
 * SymbolicTerm. The matrix S is diagonal with elements \sin\sigma_j, where the index j refers to a coordinate space
 * lattice site. This matrix arises from the Hubbard-Stratonovich form of the contact interaction, and its form will
 * need to change in a non-trivial way if we would like to introduce a new interaction.
 *
 * The pretty-printed form of MatrixS is "S_( a, b )", where ( a, b ) are the coordinate-space indices of the matrix.
 * Since the matrix is diagonal, a \delta_{a, b} replaces the matrix under path integration.
 *
 * Note that since no additional data members are added to this derived class, the compiler-provided copy constructor
 * and assignment operator overloads are sufficient.
 */
class MatrixS : public SymbolicTerm {

public:

	/**
	 * Default constructor for MatrixS.
	 */
	MatrixS();

	/**
	 * Gets the pretty-printed representation of this instance of MatrixS. The the class description for further details
	 * on the representation.
	 * @return The pretty-printed representation.
	 */
	const std::string to_string() const;

	/**
	 * Generates a deep copy of this instance on the heap and returns a smart shared pointer to the copy.
	 * @return A smart shared pointer to the generated copy.
	 */
	SymbolicTermPtr copy();
};

/**
 * Symbolic representation of the scalar A in terms of the current perturbation theory formalism. The class is
 * derived from SymbolicTerm. This scalar parameter, which is the expansion variable, is defined such that
 *
 * A \equiv \sqrt{2( e^{\tau g} - 1 )}
 *
 * where \tau is the temporal lattice spacing, and g is the bare coupling. The pretty-printed representation is
 * always simply "A". Instances of TermA are simply used to track the order of terms in the perturbation expansion.
 */
class TermA : public SymbolicTerm {
public:

	/**
	 * Default constructor for TermA.
	 */
	TermA();

	/**
	 * Gets the pretty-printed representation of TermA. Always returns "A".
	 * @return String of the pretty-printed representation of TermA.
	 */
	const std::string to_string() const;

	/**
	 * Equality operator overload for TermA. Always returns true since all TermA instances are mathematically
	 * equivalent.
	 * @return Boolean value of true.
	 */
	bool operator==( const TermA &other ) const;

	/**
	 * Generates a deep copy of this instance on the heap and returns a smart shared pointer to the copy.
	 * @return A smart shared pointer to the generated copy.
	 */
	SymbolicTermPtr copy();
};

/**
 * Symbolic representation of the scalar E_k in terms of the current perturbation formalism, where E_k is the trace
 * over the k-th order matrix product KS, such that
 *
 * E_k \equiv \frac{(-1)^{k + 1}{k} \tr\left[(K\mathcal{S}[\sigma])^k\right]
 */
class TermE : public SymbolicTerm {

public:

	/**
	 * Constructor for TermE. The order k of the term must always be specified where k is an integer and k >= 1.
	 * @param thisOrder Order k of the term E_k.
	 */
	TermE( int thisOrder );

	/**
	 * Constructor for TermE. Additionally specifies the particle flavor label that is assigned to all matricies K
	 * that appear in the product (KS)^k.
	 * @param thisOrder Order k of the term E_k.
	 * @param thisFlavorLabel Flavor label to be assigned to all MatrixK instances that appear in this term.
	 */
	TermE( int thisOrder, const char* thisFlavorLabel );

	/**
	 * Gets the pretty-printed representation TermE. Returns "E_k" where k is replaced by the value of thisOrder.
	 * @return String of the pretty-printed representation of TermE.
	 */
	const std::string to_string() const;

	/**
	 * Equality operator overload for TermE. Two instances of TermE are considered equal if the orders and flavor labels
	 * of this instance and other are both equal.
	 * @param other Other instance of TermE to compare for equality.
	 * @return true if order and flavorLabel of both instances are equal, false otherwise.
	 */
	bool operator==( const TermE &other ) const;

	/**
	 * Generates a deep copy of this instance on the heap and returns a smart shared pointer to the copy.
	 * @return A smart shared pointer to the generated copy.
	 */
	SymbolicTermPtr copy();

	/**
	 * Gets the order of this instance of TermE.
	 * @return The order of this term.
	 */
	unsigned int getOrder();

	/**
	 * Generates the full representation of the term E_k at the previously assigned order. The full representation is
	 * in terms of a scalar coefficient, a trace operator, and instances of MatrixK and MatrixS with the given flavor
	 * label assignments. The expression is placed on the heap with a smart pointer to the object returned, which is
	 * an instance of ProductPtr.
	 * @return A shared smart pointer to the full expression representation.
	 */
	SymbolicTermPtr getFullExpression();

private:

	/**
	 * Order of the E term, or the order of the product KS under the trace. Always is an integer greater then zero.
	 */
	unsigned int order;

};

/**
 * Symbolic representation of a fraction or ratio between two scalar values. The fraction is specified in terms of a
 * numerator and denominator.
 */
class CoefficientFraction : public SymbolicTerm {

public:

	/**
	 * Default constructor for CoefficientFraction. Sets the fraction equal to zero, where the numerator is 0 and the
	 * denominator is 1.
	 */
	CoefficientFraction();

	/**
	 * Constructs an instance of CoefficientFraction with a specified numerator and denominator.
	 * @param n Numerator of the fraction.
	 * @param d Denominator of the fraction.
	 */
	CoefficientFraction( double n, double d );

	/**
	 * Gets the pretty-printed representation of this instance of CoefficientFraction. Returns "n / d", where n is
	 * the string representation of the numerator, and d is the string representation of the denominator (as defined by
	 * operator<<).
	 * @return The pretty-printed string representation of this instance.
	 */
	const std::string to_string() const;

	/**
	 * Multiplication operator overload between two instances of CoefficientFraction. Returns the resulting product
	 * as a new CoefficientFraction which has been simplified by a call to reduce(). Note that any two
	 * CoefficientFraction objects commute.
	 * @param obj The CoefficientFraction term to multiply against this instance.
	 * @return The product of this instance and obj.
	 */
	CoefficientFraction operator*( const CoefficientFraction& obj ) const;

	/**
	 * Addition operator overload between two instances of CoefficientFraction. Note that any two CoefficientFraction
	 * objects commute. Returns the resulting sum as a new CoefficientFraction which as been simplified by a call to
	 * reduce().
	 * @param obj The CoefficientFraction term to add against this instance.
	 * @return The sum of this instance and obj.
	 */
	CoefficientFraction operator+( const CoefficientFraction& obj ) const;

	/**
	 * Multiplication operator overload between this instance of CoefficientFraction and another CoefficientFloat
	 * object. Note that CoefficientFraction and CoefficientFloat objects commute. Returns the resulting product as a
	 * new CoefficientFraction which has been simplified by a call to reduce().
	 * @param obj The CoefficientFloat instance to multiply against this instance.
	 * @return The product of this instance and obj.
	 */
	CoefficientFraction operator*( const CoefficientFloat& obj ) const;

	/**
	 * Addition operator overload between this instance of CoefficientFraction and another CoefficientFloat object.
	 * Note that CoefficientFraction and CoefficientFloat objects commute. Returns the resulting sum as a new
	 * CoefficientFraction which has been simplified by a call to reduce().
	 * @param obj The CoefficientFloat instance to add against this instance.
	 * @return The sum of this instance and obj.
	 */
	CoefficientFraction operator+( const CoefficientFloat& obj ) const;

	/**
	 * Operator *= overload between two instances of CoefficientFraction. Implemented in terms of operator*.
	 * @param obj The CoefficientFraction instance to multiply against this instance.
	 * @return This instance of CoefficientFraction, whose value is the previous value of this instance times the value
	 *         of obj.
	 */
	CoefficientFraction& operator*=( const CoefficientFraction& obj );

	/**
	 * Operator += overload between two instances of CoefficientFraction. Implemented in terms of operator+.
	 * @param obj The CoefficientFraction instance to add against this instance.
	 * @return This instance of CoefficientFraction, whose value is the previous value of this instance plus the value
	 *         of obj.
	 */
	CoefficientFraction& operator+=( const CoefficientFraction& obj );

	/**
	 * Generates a deep copy of this instance on the heap and returns a smart shared pointer to the copy.
	 * @return A smart shared pointer to the generated copy.
	 */
	SymbolicTermPtr copy();

	/**
	 * Gets the floating-point numerical result of dividing the numerator by the denominator.
	 * @return The numerical evaluation of num / den.
	 */
	double eval() const;

	/**
	 * Reduces this fraction to lowest terms by the greatest common divisor.
	 */
	void reduce();

private:

	/**
	 * Numerator (num) and denominator (den) of this fraction.
	 */
	double num, den;
};

/**
 * Symbolic representation of a signed floating-point coefficient in terms of the current perturbation theory
 * formalism. An instance of this class is most typically used to represent a coefficient of -1.
 */
class CoefficientFloat : public SymbolicTerm {

	friend class CoefficientFraction; // TODO: Replace friend declaration with making eval() const.

public:

	/**
	 * Constructor for CoefficientFloat. Accepts the value of the floating-point coefficient.
	 * @param val Numerical value of this floating-point coefficient.
	 */
	CoefficientFloat( double val );

	/**
	 * Gets the pretty-printed representation of this term. Always returns the string reprentation of the value. A copy
	 * of the value is truncated below a particular threshold to avoid machine epsillon precision errors. This
	 * threshold value is set to 1e-10.
	 * @return The string representation of this CoefficientFloat instance.
	 */
	const std::string to_string() const;

	/**
	 * Multiplication operator overload for products between two instances of CoefficientFloat.
	 * @param obj Other CoefficientFloat instance to multiply with this instance.
	 * @return A new CoefficientFloat instance whose value is the product of this instance's value and the value of obj.
	 */
	CoefficientFloat operator*( const CoefficientFloat& obj ) const;

	/**
	 * Addition operator overload for sums between two instances of CoefficientFloat.
	 * @param obj Other CoefficientFloat instance to add with this instance.
	 * @return A new CoefficientFloat instance whose value is the sum of this instance's value and the value of obj.
	 */
	CoefficientFloat operator+( const CoefficientFloat& obj ) const;

	/**
	 * Multiplication operator overload for products between an instance of CoefficientFloat and an instance of
	 * CoefficientFraction.
	 * @param obj CoefficientFraction instance to multiply with this instance of CoefficientFloat.
	 * @return A new CoefficientFraction instance whose value is the product of this instance's value and the value
	 * of obj.
	 */
	CoefficientFraction operator*( const CoefficientFraction& obj ) const;

	/**
	 * Addition operator overload for sums between an instance of CoefficientFloat and an instance of
	 * CoefficientFraction.
	 * @param obj CoefficientFraction instance to add with this instance of CoefficientFloat.
	 * @return A new CoefficientFraction instance whose value is the sum of this instance's value and the value of obj.
	 */
	CoefficientFraction operator+( const CoefficientFraction& obj ) const;

	/**
	 * Multiplication-equals operator overload for products between this instance of CoefficientFloat and another
	 * instance of CoefficientFloat.
	 * @param obj CoefficientFloat instance to multiply with this instance of CoefficientFloat.
	 * @return A reference to this instance of CoefficientFloat after the product is performed.
	 */
	CoefficientFloat& operator*=( const CoefficientFloat& obj );

	/**
	 * Addition-equals operator overload for sums between this instance of CoefficientFloat and another instance of
	 * CoefficientFloat.
	 * @param obj CoefficientFloat instance to add with this instance of CoefficientFloat.
	 * @return A reference to this instance of CoefficientFloat after the sum is performed.
	 */
	CoefficientFloat& operator+=( const CoefficientFloat& obj );

	/**
	 * Generates a deep copy of this instance on the heap and returns a smart shared pointer to the copy.
	 * @return A smart shared pointer to the generated copy.
	 */
	SymbolicTermPtr copy();

	/**
	 * Gets the floating-point numerical value of this instance. Implemented here for polymorphism purposes.
	 * @return The numerical value of this instance.
	 */
	double eval() const;

private:

	/**
	 * The value this floating-point coefficient holds.
	 */
	double value;
};

/*
 * ***********************************************************************
 * TERM DERIVED CLASSES - EXPRESSION CLASSES
 * ***********************************************************************
 */

/**
 * Symbolic representation of a sum of terms. Sum instances may be nested within other Sums, such that they are
 * "sums of one or more sums". Sums may be simplified and expanded through a variety of provided methods.
 */
class Sum : public SymbolicTerm {

	friend bool unpackTrivialExpression( SymbolicTermPtr & );

	friend Sum truncateOddOrders( SymbolicTermPtr expr );

	friend class Product;

public:

	/**
	 * Default constructor. Creates an empty sum.
	 */
	Sum();

	/**
	 * Constructs a Sum whose terms are the elements of the passed vector.
	 * @param thisTerms A vector of SymbolicTermPtr objects which each reference a term in this sum.
	 */
	Sum( std::vector<SymbolicTermPtr> thisTerms );

	/**
	 * Constructs a Sum and pushes back the single passed term into the sum.
	 * @param term The term to add to the empty sum.
	 */
	Sum( SymbolicTermPtr term );

	/**
	 * Copy constructor for a Sum.
	 * @param s The Sum to copy.
	 */
	Sum( const Sum &s );

	/**
	 * Destructor.
	 */
	~Sum();

	/**
	 * Assignment operator overload for a Sum.
	 * @param rhs The Sum on the right-hand side of the assignment.
	 * @return A reference to this instance.
	 */
	Sum& operator=( const Sum &rhs );

	/**
	 * Gets a pretty-printed representation of this Sum.
	 * @return A string that holds the pretty-printed representation of this Sum.
	 */
	const std::string to_string() const;

	/**
	 * Generates a deep copy of this instance on the heap and returns a smart shared pointer to the copy.
	 * @return A smart shared pointer to the generated copy.
	 */
	SymbolicTermPtr copy();

	/**
	 * Mathematically simplifies the sum by removing any CoefficientFloat or CoefficientFraction terms which evaluate
	 * to zero. Trivial terms are unpacked in the process.
	 */
	void simplify();

	/**
	 * Simplifies the structure representation of this sum by recursively unpacking trivial expressions.
	 */
	void reduceTree();

	/**
	 * Computes and returns the fully expanded expression.
	 * @return The expanded expression.
	 */
	Sum getExpandedExpr();

	/**
	 * Adds a term to teh end of the sum.
	 * @param thisTerm A SymbolicTermPtr which references the term to add to this sum.
	 */
	void addTerm( SymbolicTermPtr thisTerm );

	/**
	 * Gets the length of the internal vector holding SymbolicTermPtr terms. This may not be the mathematical
	 * number of terms in the sum, unless it is fully expanded and has a fully reduced tree.
	 * @return The number of terms in the sum as currently represented.
	 */
	int getNumberOfTerms();

    SymbolicTermPtr getTerm( int i );

	/**
	 * Calls the reduceFourierSumIndices() method on each term of the sum which is a Product. The expression
	 * should be expanded and reduced before using this method.
	 */
	void reduceFourierSumIndices();

	/**
	 * Clears all terms in this sum.
	 */
	void clear();

	/**
	 * Calls the combineCoefficients() method on each term of the sum which is a Product. The expression
	 * should be expanded and reduced before using this method.
	 */
	void combineCoefficients();

	/**
	 * Gets the beginning iterator of the internal vector holding the terms of this sum.
	 * @return The beginning iterator of this sum.
	 */
	std::vector<SymbolicTermPtr>::const_iterator getIteratorBegin() const;

	/**
	 * Gets the ending iterator of the internal vector holding the terms of this sum.
	 * @return The ending iterator of this sum.
	 */
	std::vector<SymbolicTermPtr>::const_iterator getIteratorEnd() const;

	/**
	 * Gets the beginning iterator of the internal vector holding the terms of this sum.
	 * @return The beginning iterator of this sum.
	 */
	std::vector<SymbolicTermPtr>::iterator getIteratorBegin();

	/**
	 * Gets the beginning iterator of the internal vector holding the terms of this sum.
	 * @return The beginning iterator of this sum.
	 */
	std::vector<SymbolicTermPtr>::iterator getIteratorEnd();

private:

	/**
	 * The vector of SymbolicTermPtr objects which reference the terms of this sum.
	 */
	std::vector<SymbolicTermPtr> terms;
};

/**
 * Symbolic representation of a product of terms (or factors). Product instances may be nested within other Products,
 * and may represent products of sums. Products may be expanded and simplified through a variety of provided methods.
 */
class Product : public SymbolicTerm {

	friend bool unpackTrivialExpression( SymbolicTermPtr & );

	friend class Sum;

	friend void indexExpression( SymbolicTermPtr expr );

public:

	/**
	 * Default constructor for a Product. Creates a product with zero factors.
	 */
	Product();

	/**
	 * Constructs a Product whose factors are the elements of the passed vector.
	 * @param thisTerms A vector of SymbolicTermPtr objects which each reference a factor in this product.
	 */
	Product( std::vector<SymbolicTermPtr> t );

	/**
	 * Constructs a Product whose only factor is the passed term.
	 * @param term A SymbolicTermPtr object which references the factor to insert into this product.
	 */
	Product( SymbolicTermPtr term );

	/**
	 * Destructor for Product instances.
	 */
	~Product();

	/**
	 * Assignment operator overload for Product objects.
	 * @param rhs The right-hand side of the assignment expression.
	 * @return A reference to this instance.
	 */
	Product& operator=( const Product &rhs );

	/**
	 * Gets the pretty-printed representation of this Product instance as a string.
	 * @return A string which hold the pretty-printed representation of this Product.
	 */
	const std::string to_string() const;

	/**
	 * Generates a deep copy of this instance on the heap and returns a smart shared pointer to the copy.
	 * @return A smart shared pointer to the generated copy.
	 */
	SymbolicTermPtr copy();

	/**
	 * Mathematically simplifies this product by evaluating the entire product to zero if any CoefficientFloat or
	 * CoefficientFraction term evaluates to zero, and by removing any CoefficientFloat or CoefficientFraction term
	 * which evaluates to one.
	 */
	void simplify();

	/**
	 * Reduces (or flattens) the tree representation of this Product by recursively unpacking trivial expressions.
	 */
	void reduceTree();

    /**
     * Computes the fully expanded expression.
     * @return The fully expanded expression.
     */
	Sum getExpandedExpr();

    /**
     * Adds a factor to the end of this product.
     * @param t The SymbolicTermPtr object which references the term to add to this product.
     */
	void addTerm( SymbolicTermPtr t );

    /**
	 * Gets the length of the internal vector holding SymbolicTermPtr terms. This may not be the mathematical
	 * number of terms in the product, unless it is fully expanded and has a fully reduced tree.
	 * @return The number of terms in the product as currently represented.
	 */
	int getNumberOfTerms();

	/**
	 * Assignment operator overload for when the right-hand side is an instance of a sum. This always returns false,
	 * which defines the fact that an instance of a Product can never be equivalent to an instance of a Sum.
	 * @param other An instance of a Sum.
	 * @return False.
	 */
	bool operator==( const Sum &other ) const;

    /**
     * Determines whether at least one of the factors in this Product is an instance of a Sum. Useful for determining
     * if it is necessary to expand the expression.
     * @return True if the Product contains a Sum, false otherwise.
     */
	bool containsSum();

    /**
     * Sets the internal vector of factors to a new vector whose only factor is a CoefficientFloat of value zero.
     */
	void zero();

    /**
     * Clears the internal vector of factors.
     */
	void clear();

	/**
	 * Iterates through each factor in the product and calls the method reduceDummyIndices() on all instances of
	 * FourierSum objects.
	 */
	void reduceFourierSumIndices();

	/**
	 * Multiplies all values of CoefficientFloat and CoefficientFraction factors together and inserts a new
	 * CoefficientFraction object with total value of all previous coefficient factors. The new object is inserted
	 * at the end of the Product.
	 */
	void combineCoefficients();

	/**
	 * Gets the beginning iterator of the internal vector holding the terms of this sum.
	 * @return The beginning iterator of this sum.
	 */
	std::vector<SymbolicTermPtr>::const_iterator getIteratorBegin() const;

	/**
	 * Gets the ending iterator of the internal vector holding the terms of this sum.
	 * @return The ending iterator of this sum.
	 */
	std::vector<SymbolicTermPtr>::const_iterator getIteratorEnd() const;

	/**
	 * Gets the beginning iterator of the internal vector holding the terms of this sum.
	 * @return The beginning iterator of this sum.
	 */
	std::vector<SymbolicTermPtr>::iterator getIteratorBegin();

	/**
	 * Gets the ending iterator of the internal vector holding the terms of this sum.
	 * @return The ending iterator of this sum.
	 */
	std::vector<SymbolicTermPtr>::iterator getIteratorEnd();

private:

	/**
	 * The vector of SymbolicTermPtr objects which reference the factors of this product.
	 */
	std::vector<SymbolicTermPtr> terms;

};

/**
 * Symbolic representation of a trace over an expression. This objects holds a pointer to a single expression object
 * whose trace is to be computed over.
 */
class Trace : public SymbolicTerm {

	friend bool isZeroTrace( SymbolicTermPtr );

	friend void indexExpression( SymbolicTermPtr expr );

public:

	/*
	 * Constructor.
	 * @param expr The pointer to an expresson whose trace is to be computed over.
	 */
	Trace( SymbolicTermPtr expr );

	/**
	 * Destructor.
	 */
	~Trace();

	/**
	 * Gets the pretty-printed representation of this Trace instance as a string.
	 * @return A string which hold the pretty-printed representation of this Trace.
	 */
	const std::string to_string() const;

	/**
	 * Generates a deep copy of this instance on the heap and returns a smart shared pointer to the copy.
	 * @return A smart shared pointer to the generated copy.
	 */
	SymbolicTermPtr copy();

	/**
	 * Calls simplify() on the object referenced by expr, whch mathematically simplifies this product by evaluating
	 * the entire product to zero if any CoefficientFloat or CoefficientFraction term evaluates to zero, and by
	 * removing any CoefficientFloat or CoefficientFraction term which evaluates to one.
	 */
	void simplify();

	/**
	 * Calls reduceTree() on the object referenced by expr, which reduces (or flattens) the tree representation of
	 * this Product by recursively unpacking trivial expressions.
	 */
	void reduceTree();

	/**
	 * Equivalance operator overload for when the right-hand side of the expression is an instance of a Trace.
	 * @param other The right-hand side of the equivalance comparison.
	 * @return True if *expr == other->expr, false otherwise.
	 */
	bool operator==( const Trace &other ) const;

	/**
	 * Less than operator overload for when the right-hand side of the expression is an instance of a Trace. Traces
	 * must be distributed over the fully expanded form of *expr. Natural ordering for Trace objects is determined by
	 * the number of terms in the fully-expanded form of *expr.
	 * @param other The right-hand side of the expression.
	 * @return True if the number of terms in this instance is less than the number of terms in other, false otherwise.
	 */
	bool operator<( const Trace &other ) const;

	/**
	 * Greater than operator for when the right-hand side of the expression is an instance of a Trace. Traces must
	 * be distributed over the fully expanded form of *expr. Natural ordering for Trace objects is determined by the
	 * number of terms in the fully-expanded form of *expr.
	 * @param other The right-hand side of the expression.
	 * @return True if the number of terms in this instance is greater than the number of terms in other, false
	 * otherwise.
	 */
	bool operator>( const Trace &other ) const;

private:

	/**
	 * Pointer to the expression object whose trace is to be taken over.
	 */
	SymbolicTermPtr expr;
};

/**
 * Symbolic representation of a Kronecker delta function of two indices.
 */
class Delta : public SymbolicTerm {
public:

	/**
	 * Constructor. Creates a standard delta function.
	 * @param a The first index.
	 * @param b The second index.
	 */
	Delta( int a, int b );

	/**
	 * Constructor. Allows for the specification
	 * @param a The first index.
	 * @param b The second index.
	 * @param type Boolean flag for whether this object is \bar{delta} (equivalent mathematically as 1 - \delta_{ij}),
	 * or not.
	 */
	Delta( int a, int b, bool type );

	/**
	 * Gets the pretty-printed representation of this Delta instance as a string.
	 * @return A string which hold the pretty-printed representation of this Trace.
	 */
	const std::string to_string() const;

	/**
	 * Generates a deep copy of this instance on the heap and returns a smart shared pointer to the copy.
	 * @return A smart shared pointer to the generated copy.
	 */
	SymbolicTermPtr copy();

	/**
	 * Equivalence operator overload for when the right-hand side of the expression is an instance of a Delta.
	 * Two Delta objects are considered equivalent if they are both barred or not barred, and if their ordered
	 * indices are equivalent.
	 * @param other The right hand side of the expression.
	 * @return True if both Delta objects are equivalent according to their indices and if it is \bar{\delta} or \delta,
	 * false otherwise.
	 */
	bool operator==( const Delta &other ) const;

	/**
	 * Returns whether this Delta instance is a \bar{\delta} or not.
	 */
	bool isDeltaBar();

private:

	/**
	 * Boolean flag that sets if this Delta instance is a \bar{\delta} or not.
	 */
	bool isBar;

};

/**
 * Symbolic representation of a Fourier sum, or a sum over coordinate space of a product of free propagators. See the
 * written documentation for the mathematical details behind this class.
 */
class FourierSum : public SymbolicTerm {
public:

	/**
	 * Constructor.
	 * @param i Vector of IndexContraction objects that are resulting momentum-space indices.
	 * @param orderInK Order in the matricies K in this expression. Should be equal to the length of the vector i.
	 */
	FourierSum( std::vector<IndexContraction> i, int orderInK );

	/**
	 * Destructor.
	 */
	~FourierSum();

	/**
	 * Gets the pretty-printed representation of this FourierSum instance as a string.
	 * @return A string which hold the pretty-printed representation of this Trace.
	 */
	const std::string to_string() const;

	/**
	 * Generates a deep copy of this instance on the heap and returns a smart shared pointer to the copy.
	 * @return A smart shared pointer to the generated copy.
	 */
	SymbolicTermPtr copy();

	/**
	 * Equivalance operator overload for when the right-hand side of the expression is an instance of a FourierSum.
	 * Two FourierSum instances are considered equivalent if their vectors of sorted indices are equivalent.
	 * @param other The right-hand side of the expression.
	 * @return True if the FourierSum instances are considered equivalent, false otherwise.
	 */
	bool operator==( const FourierSum &other ) const;

	/**
	 * Reduces dummy indices such that extraneous free indices are eliminated and an equivalent diagram is written in
	 * terms of the fewest number of indices possible.
	 */
	void reduceDummyIndices();

	/**
	 * Returns the vector of contracted indices.
	 * @return
	 */
    std::vector<IndexContraction> getContractionVector();

private:

    /**
     * The vector of contracted momentum-space indices in this Fourier sum over coordinate space.
     */
	std::vector<IndexContraction> indices;

    /**
     * Order in matrices K of this product.
     */
	int order;

};

/*
 * ***********************************************************************
 * GENERIC HELPER FUNCTIONS
 * ***********************************************************************
 */

/**
 * Determines if this expression is a trivial one; that is, this instance is a trivial sum with a single term.
 * If it is, the SymbolicTermPtr is redirected to the nested expression.
 * @param expr A SymbolicTermPtr that references a Product or Sum which should be unpacked.
 * @return Returns true if the expression was trivial and a simplification was performed, false otherwise.
 */
bool unpackTrivialExpression( SymbolicTermPtr& expr );

/**
 * Determines if the argument of the passed trace symbolically evaluates to zero. The passed argument must reference
 * an instance of a Trace. A trace is considered to evaluate to zero if the contained expression contains a Sum or
 * Product with zero terms, or if the string representation of the contained expression evaluates to "0" or "{0}".
 * One should usually call expr->simplify() before calling this function.
 * @param tr A pointer to a Trace object which is to be checked.
 * @return True if the trace symbolically evaluates to zero, false otherwise.
 */
bool isZeroTrace( SymbolicTermPtr tr );

/**
 * Determines the order of a Product in TermA objects (e.g., determine n for A^n in a product).
 * @param prod A pointer to a Product object to be evaluated.
 * @return The order n for A^n of the expression.
 */
int getProductAOrder( SymbolicTermPtr prod );

int getDualProductAOrder( Product &prodA, Product &prodB );

std::map<std::string, int> getFlavorLabelOrder( SymbolicTermPtr prod );

bool areFlavorLabelOrdersIdentical( std::map<std::string, int> orderA, std::map<std::string, int> orderB );

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

bool areDiagramsSimilar( std::vector<IndexContraction> diagramA, std::vector<IndexContraction> diagramB );

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

std::ostream& operator<<( std::ostream& os, const std::map<std::string, int> &obj );

#endif /* PTSYMBOLICOBJECTS_H_ */
