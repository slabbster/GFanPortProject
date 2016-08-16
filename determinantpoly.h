#ifndef DETERMINANTPOLY_H_DEFINED
#define DETERMINANTPOLY_H_DEFINED

#include "polynomial.h"

/**
 * Computes the initial form w.r.t. w of the determinant of the square matrix whose entries are the entries of g.
 * If takeDerivatives is true, then the matrix to be considered is the Jacobi matrix of g, which is then considered as vector.
 */
Polynomial initialFormOfDeterminant(PolynomialSet const &g, IntegerVector const &w, bool takeDerivatives);

/**
 * Computes all the codim by codim minors of the Jacobi matrix of g.
 */
PolynomialSet jacobiMinors(PolynomialSet const &g, int codim);

#endif


