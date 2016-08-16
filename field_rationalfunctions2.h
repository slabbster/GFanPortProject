#ifndef FIELD_RATIONALFUNCTIONS2_INCLUDED
#define FIELD_RATIONALFUNCTIONS2_INCLUDED

#include "field.h"
#include <string>
#include "polynomial.h"

using namespace std;

class FieldRationalFunctions2Implementation : public FieldImplementation
{
  PolynomialRing thePolynomialRing;
  FieldElementImplementation *zHomomorphismImplementation(int n);/* Creates FieldElementImplementation object with refcount1 */
  FieldElement zHomomorphism(int n);
  const char *name();
  std::string toString()const;
public:
  virtual bool isRationals()const;
  PolynomialRing getPolynomialRing()const;
  FieldRationalFunctions2Implementation(PolynomialRing const &r);
//  PolynomialRing getPolynomialRing()const;
};



// Let's see how inheritance and slicing works together
class FieldRationalFunctions2 : public Field
{
public:
  FieldRationalFunctions2(PolynomialRing const &r);
  FieldElement polynomialToFraction(Polynomial const &p);
};

/**
 * Creates a polynomial ring where
 * After having
 */
PolynomialRing makeVariablesParameters(PolynomialRing const &r, int numberOfParameters);
Polynomial makeVariablesParameters(PolynomialRing const &genericRing, Polynomial const &p);
PolynomialSet makeVariablesParameters(PolynomialRing const &genericRing, PolynomialSet const &p);

#endif
