#ifndef SINGULARCONVERSION_H_INCLUDED
#define SINGULARCONVERSION_H_INCLUDED

#include <assert.h>

#define OM_NDEBUG
#define NDEBUG
#include "mod2.h"
#include "structs.h" // Singular structs
#include "ring.h"
#include "numbers.h"
#include "polys.h"
#include "longrat.h"
#include "ideals.h"
#include "kstd1.h"
#undef NDEBUG

#include "polynomialring.h"
#include "polynomial.h"
#include "field_rationals.h"
#include "printer.h"
#include "log.h"
#include <iostream>

ring singularRing(PolynomialRing const &r);
void freeSingularRing(ring R);
poly singularPolynomial(Polynomial const &p);
ideal singularPolynomialSet(PolynomialSet const &g);
FieldElement fromSingularCoefficient(PolynomialRing const &r, number c);
Polynomial fromSingularPolynomial(PolynomialRing const &r, poly &p);
PolynomialSet fromSingularIdeal(PolynomialRing const &r, ideal i);

#endif
