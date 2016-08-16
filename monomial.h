#ifndef MONOMIAL_H_INCLUDED
#define MONOMIAL_H_INCLUDED

#include <string>
#include "vektor.h"
#include "polynomialring.h"

using namespace std;

class Monomial
{
  PolynomialRing theRing;
 public:
  int testtest;
  IntegerVector exponent;
  Monomial(PolynomialRing const &r):
    theRing(r),
    exponent(r.getNumberOfVariables())/* In Gfan 0.4 the exponent was not initialized. Does this slow the program now?*/
    {
    }
    Monomial(PolynomialRing const &r, IntegerVector const &v);
  PolynomialRing const &getRing()const{return theRing;}
  string toString(bool alwaysWriteSign, bool writeIfOne, bool latex/*, bool mathMode*/)const;
};

#endif
