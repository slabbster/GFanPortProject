#ifndef TROPICALTRAVERSE_H_INCLUDED
#define TROPICALTRAVERSE_H_INCLUDED

#include "bergman.h"


/**
   Represents a tropical variety up to symmetry.

   New strategy.
 */

PolyhedralFan tropicalTraverse(PolynomialSet coneGroebnerBasis, PolynomialSet idealGroebnerBasis, SymmetryGroup const *symmetryGroup=0);
void coneChangeDebugger(PolynomialSet coneGroebnerBasis, PolynomialSet idealGroebnerBasis, IntegerVectorList const &ridges, IntegerVectorList const &rays);

#endif
