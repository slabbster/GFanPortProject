#ifndef SINGULAR_H_INCLUDED
#define SINGULAR_H_INCLUDED

#include "polynomial.h"
#include "termorder.h"

void singularBuchberger(PolynomialSet *g, TermOrder const &termOrder);

#endif
