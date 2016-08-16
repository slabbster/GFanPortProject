#ifndef TROPICALBASIS_H_INCLUDED
#define TROPICALBASIS_H_INCLUDED

#include "polynomial.h"
#include "polyhedralfan.h"

PolynomialSet tropicalBasisOfCurve(int n, PolynomialSet g, PolyhedralFan *intersectionFan=0, int linealitySpaceDimension=-1);
 /*Assuming g is homogeneous and has dimension 1 modulo homogeneous
   space. intersectionFan specifies a fan to use for the temporary
   hypersurface intersections. This is useful since we then don't have
   to call tropicalPrincipalIntersection to get the variety after
   having computed a tropical basis. No PolyhedralFan needs to be
   specified. If specified, the fan does not have to be initialized in
   any way. Its data is simply overwritten.*/

#endif
