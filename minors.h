#ifndef MINORS_H_INCLUDED
#define MINORS_H_INCLUDED
#include "polynomial.h"

PolynomialSet minors(PolynomialRing const &R, int r, int d, int n, bool withNames, bool M2Convention);

#endif
