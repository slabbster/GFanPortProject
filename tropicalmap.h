#ifndef TROPICALMAP_H_INCLUDED
#define TROPICALMAP_H_INCLUDED

#include "polyhedralfan.h"
#include "polynomial.h"

PolyhedralCone imageOfConeUnderTropicalMap(PolynomialSet const &polynomialMap, PolyhedralCone const &cone);
PolyhedralFan imageOfTropicalMapGivenLinearityCones(PolynomialSet const &polynomialMap, PolyhedralFan linearityCones);

/**
   This routine computes the image of the support of the fan domain
   under the tropicalization of the polynomial map.
 */
PolyhedralFan imageOfTropicalMap(PolynomialSet const &polynomialMap, PolyhedralFan domain);

#endif
