#ifndef INTSINPOLYTOPE_H_INCLUDED
#define INTSINPOLYTOPE_H_INLCUDED

#include "vektor.h"
#include "matrix.h"

bool solveIntegerProgramIneq(IntegerMatrix const &M, IntegerVector const &rightHandSide, IntegerVector &solution);

/** Returns the integer points in a polytope.
    The polytope is of the form Mx<=rightHandSide.
    One point p in the polytope must be given as input.
 */
IntegerVectorList intsInPolytopeGivenIneqAndPt(IntegerMatrix const &M, IntegerVector const &rightHandSide, IntegerVector const &p);
IntegerVectorList intsInPolytopeGivenIneq(IntegerMatrix const &M, IntegerVector const &rightHandSide);

#endif
