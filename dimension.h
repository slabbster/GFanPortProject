#ifndef DIMENSION_H_INCLUDED
#define DIMENSION_H_INCLUDED

#include "polynomial.h"

PolynomialSet radicalOfMonomialIdeal(PolynomialSet const &monomialGenerators);
int krullDimensionOfMonomialIdeal(PolynomialSet const &monomialGenerators);
int krullDimension(PolynomialSet const &groebnerBasis);

#endif
