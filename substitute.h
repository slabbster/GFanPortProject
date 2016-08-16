#ifndef SUBSTITUTE_H_INCLUDED
#define SUBSTITUTE_H_INCLUDED

#include "polynomial.h"
#include "matrix.h"

Polynomial multiplicativeChange(Polynomial const &p, IntegerMatrix const &mat);//marking not preserved
PolynomialSet multiplicativeChange(PolynomialSet const &g, IntegerMatrix const &mat);//markings not preserved

#endif
