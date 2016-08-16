#ifndef MINKOWSKIDUAL_H_INCLUDED
#define MINKOWSKIDUAL_H_INCLUDED

#include "polyhedralfan.h"
#include "polynomial.h"
#include "symmetriccomplex.h"

/*
  The dual Minkowski routine is used for extracting the normal cones
  of the mixed faces of symmetric Minkowski sum. The point is that in
  order not to compute all faces one must do the face extraction in
  the oppisite order of what is currently done in PolyhedralFan.
  Equivalently, in the usual order for the dual polytope.
 */

/*
  cones contains only full dimensional cones.
 */
SymmetricComplex dualMinkowskiMixed(PolynomialSet const &g, SymmetryGroup const &sym, PolyhedralFan const &cones);

#endif
