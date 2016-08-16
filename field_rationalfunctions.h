#ifndef FIELD_RATIONALFUNCTIONS_H_INCLUDED
#define FIELD_RATIONALFUNCTIONS_H_INCLUDED

#include "field.h"
#include <string>
#include "polynomial.h"

using namespace std;

class FieldRationalFunctionsImplementation : public FieldImplementation
{
  Field coefficientField;
  string parameterName;
  PolynomialRing thePolynomialRing;
  FieldElementImplementation *zHomomorphismImplementation(int n);/* Creates FieldElementImplementation object with refcount1 */
  FieldElement zHomomorphism(int n);
  const char *name();
  std::string toString()const;
public:
  virtual bool isRationals()const;
  FieldRationalFunctionsImplementation(Field const &f_, string const &parameterName_);
  PolynomialRing getPolynomialRing()const;
};



// Let's see how inheritance and slicing work together
class FieldRationalFunctions : public Field
{
public:
  FieldRationalFunctions(Field const &coefficientField, string const &parameterName);
  FieldElement exponent(int power);
  FieldElement fromCoefficientField(FieldElement const &c);
  /**
   This function gives a value in the coefficient field by substituting the perturbation variable by tvalue.
   In case the denominator becomes zero during this process the routine asserts.
   */
  FieldElement substitute(FieldElement const &e, FieldElement const &tvalue)const;
};

void testRationalFunctionField();//test routine
#endif
