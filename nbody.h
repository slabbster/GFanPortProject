#ifndef NBODY_H_INCLUDED
#include "polynomial.h"

PolynomialSet AlbouyChencinerEquations(int N, bool withMasses, bool symmetric=true, bool withSVariables=false, bool saturate=true);
PolynomialSet DziobekEquations(PolynomialRing const &r, int N, bool withMasses, bool withSVariables=false, bool saturate=true);
PolynomialSet nbodyDeterminants(PolynomialRing const &r, int N, bool withMasses, int determinantSize);
Polynomial massEquation(PolynomialRing const &r, int N, bool withMasses, bool saturate=true);
PolynomialSet SEquations(PolynomialRing const &r, int N, bool withMasses, bool saturate=true);

#endif
