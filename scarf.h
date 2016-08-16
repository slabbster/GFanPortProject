#ifndef SCARF_H_INCLUDED
#define SCARF_H_INCLUDED

#include "matrix.h"

IntegerVectorList neighbours(IntegerMatrix const &A);
bool satisfiesA1(IntegerMatrix const &A);
bool satisfiesA2(IntegerMatrix const &A);
bool satisfiesA3(IntegerMatrix const &A, IntegerVectorList const *N=0);
/** Checks Scarf's genericity condition for jth row of the matrix A*/
bool satisfiesA3i(IntegerMatrix const &A, int j, IntegerVectorList const *N=0);
IntegerVectorList orientedNeighbours(IntegerVectorList const& N, IntegerVector const &v);
IntegerVector kFlip(IntegerMatrix const &A, IntegerMatrix const &N, IntegerVector simplex, int vertex);
void traverseScarfComplex(IntegerMatrix const &A, IntegerMatrix const &N, IntegerVector simplex);
IntegerVector computeMaximalScarfSimplex(IntegerMatrix const &A, IntegerMatrix const &N);

#endif
