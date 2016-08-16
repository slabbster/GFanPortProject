#ifndef FIELD_RATIONALS_H_INCLUDED
#define FIELD_RATIONALS_H_INCLUDED

#include <vector>
#include <cstdio> /* Always include cstdio before gmp.h.*/
#include <gmp.h> //remove?
#include "vektor.h"
#include "field.h"


class FieldRationalsImplementation : public FieldImplementation
{
  FieldElementImplementation *zHomomorphismImplementation(int n);/* Creates FieldElementImplementation object with refcount1 */
  FieldElement zHomomorphism(int n);
  const char *name();
  std::string toString()const;
public:
  virtual bool isRationals()const;
  FieldRationalsImplementation();
};

/*
class FieldRationals : public Field
{
public:
  FieldRationals();
};
*/
		       //extern FieldRationals Q;

extern Field Q;

IntegerVector primitiveVector(vector<FieldElement> const &v);
IntegerVector toIntegerVector(vector<FieldElement> const &v);// MUST HAVE INTEGER ENTRIES
int toInteger(FieldElement const &a);// MUST BE INTEGER
mpq_t *fieldElementToGmp(FieldElement const &c);
double fieldElementToFloatingPoint(FieldElement const&c);
/**
 * Will create a new FieldElement containing the value of c. The gmp representation of c is copied,
 * so c must be destructed by the caller.
 */
FieldElement fieldElementFromGmp(mpq_t *c);
/**
 * Will create a new FieldElement containing the value of c. The gmp representation of c is copied,
 * so c must be destructed by the caller.
 */
FieldElement fieldElementFromGmpZ(const mpz_t *c);
FieldElement gcd(FieldElement const &a, FieldElement const &b, FieldElement &s, FieldElement &t);
/**
 * Assumes that a and b are integers in the field of rational numbers.
 * Assumes that b is non-zero.
 * The remainder (if pointer is non-zero) is set to the remainder
 * of a divided by b in the interval [0,B[, where B=abs(b)=max(b,-b).
 * The return value equals the integer (a-remainder)/b.
 */
FieldElement integerDivision(FieldElement const &a, FieldElement const &b, FieldElement *remainder=0);
#endif
