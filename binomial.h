#ifndef BINOMIAL_H_INCLUDED
#define BINOMIAL_H_INCLUDED

#include "vektor.h"
#include "polynomial.h"

IntegerVector binomialToIntegerVector(Polynomial const &p);
Polynomial integerVectorToBinomial(IntegerVector const &v, PolynomialRing const &r);

#endif
