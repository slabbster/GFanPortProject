#ifndef POLYNOMIAL_H_INCLUDED
#define POLYNOMIAL_H_INCLUDED

class Polynomial;
class PolynomialSet;
//class PolynomialSetList;

#include <map>
#include "polynomialring.h"
//#include "linalg.h"
#include "term.h"
#include "termorder.h"


struct TermMapCompare
{
  inline bool operator()(const Monomial &s1, const Monomial &s2) const
  {
    return LexicographicTermOrder()(s1.exponent,s2.exponent);
  }
};

typedef std::map<Monomial,FieldElement,TermMapCompare> TermMap;

class Polynomial
{
  Term marked;
  bool isMarkedBool;
  int sugar;
  PolynomialRing theRing;
  /**
   * Used for splitting polynomials into two.
   */
  Polynomial half(bool secondHalf)const;
 public:
  PolynomialRing const &getRing()const{return theRing;}
  TermMap terms;
  Polynomial(const Term &t);
  Polynomial(PolynomialRing const &r);
  void computeInitialSugar();
  int getSugar()const;
  void madd(const Term &m, const Polynomial &p);
  void operator+=(const Polynomial &p);
  void operator-=(const Polynomial &p);
  friend Polynomial operator+(const Polynomial &p, const Polynomial &q);
  friend Polynomial operator-(const Polynomial &p, const Polynomial &q);
  void operator*=(const Term &p);
  void operator*=(const Monomial &m);
  void operator*=(FieldElement const &c);
  void operator*=(Polynomial const &p);
  friend Polynomial operator*(const Polynomial &p, Term const &t);
  friend Polynomial operator*(const Polynomial &p, Monomial const &m);
  friend Polynomial operator*(const Polynomial &p, FieldElement const &c);
  friend Polynomial operator*(const Polynomial &p, const Polynomial &q);
  int getNumberOfVariables()const;
  void changeNumberOfVariables(int n);
  void mark(class TermOrder const &termOrder);
  void mark(Monomial const &monomial);
  /**
   * Changes the marking of *this to that of p.
   */
  void copyMarking(Polynomial const &p);
  void scaleMarkedCoefficientToOne();
  bool isMonomial() const;
  const Term &getMarked()const{return marked;}
  bool isZero()const;
  int numberOfTerms()const;//linear time in worst case!
  IntegerVector exponentsSum()const;
  IntegerVector greatestCommonMonomialDivisor()const;
  int totalDegree()const;
  int64 degree(IntegerVector const &w)const;// evaluateTropically
  Polynomial homogenization(PolynomialRing const &newRing, IntegerVector const *w=0)const;
  /**
     Each variable x_i in the polynomial is substituted by w_ix_i.
   */
  Polynomial torusAct(class FieldVector const &w)const;
//  Polynomial derivative()const;//single variable
  Polynomial derivative(int j=0)const;//with respect to jth variable
  Polynomial deHomogenization()const;
  Polynomial deHomogenizationInSameRing()const;
  int numberOfVariablesInRing()const;//by looking at just one monomial
  bool checkMarking(TermOrder const &termOrder)const;
  /**
     Checks wether the polynomial is homogeneous with respect to the
     grading given by the vector.
   */
  bool isHomogeneous(IntegerVector const &v)const;
  bool isHomogeneous()const;
  void saturate(int variableNum=-1);//default is to saturate with respect to all variables
  bool isMarked()const;
  bool isValid(int numberOfVarialbes=-1)const; //debug testing
  FieldElement evaluate(const FieldElement &x)const;//single variable
/**
 * This routine returns a vector with one entry for each variable in the ring. The entry is 1
 * if the variable is used in the polynomial - 0 otherwise.
 */
  IntegerVector usedVariables()const;
  int maximalIndexOfVariableInSupport()const;
  Polynomial embeddedInto(PolynomialRing const &r2, list<int> const *chosenVariables=0)const;
  /* Creates a polynomial in r2 by extending or truncating the
     exponent vectors of the current polynomial appropriately. In
     particular this method can be used for dehomogenization */
  /**
   * Returns the index of the last variable used in this polynomial +1.
   */
  int numberOfVariablesInUseConsecutive()const;
  string toString(bool latex/*, bool mathMode*/)const;
  /**
   * For debugging
   */
  bool checkExponentVectors()const;
  double evaluateFloat(FloatVector const &x)const;
  ComplexNumber evaluateComplex(ComplexVector const &x)const;
};


class PolynomialCompare
{
 public:
  bool operator()(const Polynomial &a, const Polynomial &b)const;
};


class PolynomialCompareMarkedTerms
{
  TermOrder const &termOrder;
 public:
  PolynomialCompareMarkedTerms(TermOrder const &termOrder_);
  bool operator()(const Polynomial &a, const Polynomial &b)const;
};


class PolynomialCompareMarkedTermsReverse
{
  TermOrder const &termOrder;
 public:
  PolynomialCompareMarkedTermsReverse(TermOrder const &termOrder_);
  bool operator()(const Polynomial &a, const Polynomial &b)const;
};

/**
 * Used for sorting according to number of terms without removing duplicates.
 */
class PolynomialCompareNumberOfTermsStable
{
 public:
  bool operator()(const Polynomial &a, const Polynomial &b)const;
};

class PolynomialSet : public list<Polynomial>
{
  PolynomialRing theRing;
 public:
  PolynomialSet(PolynomialRing const &r):theRing(r){}
  PolynomialRing const &getRing()const
    {
      return theRing;
    }
  /**
   * Returns the index of the last variable used in an element of the set +1.
   */
  int numberOfVariablesInUseConsecutive()const;
  void changeNumberOfVariables(int n=-1);
  int numberOfVariablesInRing()const;//by looking at just one monomial
  //  Field &
  void scaleMarkedCoefficientsToOne();
  void mark(class TermOrder const &termOrder);
  void markAndScale(TermOrder const &termOrder);
  /**
   * Changes the markings of *this to that of g.
   */
  void copyMarkings(PolynomialSet const &g);

  bool checkMarkings(TermOrder const &termOrder)const;
  /**
     The set of vectors inducing the markings of the polynomials in
     the set is an open polyhedral cone (assuming that the markings
     are consistent). This function checks if the vector v is
     contained in the closure of that cone.
   */
  bool containsInClosedGroebnerCone(IntegerVector const &v)const;
  /**
     Checks whether the polynomials are homogeneous with respect to the
     grading given by the vector.
   */
  bool isHomogeneous(IntegerVector const &v)const;
  bool isHomogeneous()const;
  void unionPolynomial(const Polynomial &p);
  void unionSet(const PolynomialSet &s);
  int totalDegree()const;
  IntegerVector exponentsSum()const;
  void sort_();
  /**
   * This routine reorders the polynomials such that polynomials with few terms go first. Duplicates are not removed.
   */
  void simplestPolynomialsFirst();
  friend bool operator==(PolynomialSet const &a, PolynomialSet const &b); //remember to call sort_ before calling this procedure
  bool isEqualTo(PolynomialSet const &b)const; //remember to call sort_ before calling this procedure
  /**
   * Assuming that the set is a Groebner basis routine tests if its ideal equals <1>.
   */
  bool isUnitIdeal()const;
  PolynomialSet markedTermIdeal()const;
  void computeInitialSugar();
  bool isMarked()const;
  PolynomialSet homogenization(PolynomialRing const &newRing, IntegerVector const *w=0)const;
  PolynomialSet torusAct(FieldVector const &w)const;
/**
 * This routine computes a homogeneity space of the polynomial set (not the ideal) and based on this
 * substitutes 1 for certain variables to make polynomial set non-homogeneous in any grading. The returned
 * ideal sits inside new ring.
 */
  PolynomialSet multiDeHomogenization()const;
  PolynomialSet deHomogenization()const;
  PolynomialSet deHomogenizationInSameRing()const;
  /**
   * This function intersects the list of polynomials with a subpolynomial ring newRing.
   * The ring newRing is embedded into the current polynomial ring as the polynomial ring
   * on the first n variables if the chosenVariables list is not specified. Otherwise it is embedded
   * as the polynomial ring on the chosenVariables given in the list.
   */
  PolynomialSet polynomialRingIntersection(PolynomialRing const &newRing, list<int> const *chosenVariables=0)const;
  bool isValid()const; //debug testing
  /**
   * Factors out common monomial factors of the elements in the set. This does not saturate the ideal
   * generated by the polynomials!
   */
  void saturate(int variableNum=-1);//default is to saturate with respect to all variables
  PolynomialSet embeddedInto(PolynomialRing const &r2, list<int> const *chosenVariables=0)const;
  int totalNumberOfTerms()const;
  FloatVector evaluateFloat(FloatVector const &x)const;
  ComplexVector evaluateComplex(ComplexVector const &x)const;
/**
 * This function removes all zero polynomials from the set.
 */
  void removeZeros();
  /**
   * This function removes duplicates from the set by first sorting.
   */
  void removeDuplicates();
};

typedef list<PolynomialSet> PolynomialSetList;
/*class PolynomialSetList : public list<PolynomialSet>
{
};*/

#endif
