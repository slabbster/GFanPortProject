#ifndef DIVISION_H_INCLUDED
#define DIVISION_H_INCLUDED

#include "polynomial.h"

Polynomial division(Polynomial p, PolynomialSet const &l, TermOrder const &termOrder, PolynomialSet *q=0);
Polynomial smartDivision(Polynomial p, PolynomialSet l, TermOrder const &termOrder);
Polynomial divisionLift(Polynomial p, PolynomialSet l, PolynomialSet lList, TermOrder const &termOrder, bool noMarking=false);

/* Note the input PolynomialSet l should be marked. The term order argument is
   only used for choosing a term in p and does not have to agree with the
   marking. */

bool isIdealContainedInIdeal(PolynomialSet const &generators, PolynomialSet const &groebnerBasis);
bool areIdealsEqual(PolynomialSet const &groebnerBasis1, PolynomialSet const &groebnerBasis2);

/* Note the second parameter for isIdealContainedInIdeal must be a marked Groebner Basis and so must both parameters for areIdealsEqual*/

/**
   Computes a vector in the relative interior of the cone of weights giving the specified marking of g. This is a useful vector to use for a termorder when doing polynomial division.
 */
IntegerVector termorderWeight(PolynomialSet const &g);

#endif
