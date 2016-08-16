#ifndef SATURATION_H_INCLUDED
#define SATURATION_H_INCLUDED

#include "polynomial.h"


PolynomialSet nonHomogeneousSaturation(PolynomialSet const &s);//Computes (I:x_1..x_n^\infty), returns marked Groebner basis

// Not implemented:
//PolynomialSet colonIdeal(PolynomialSet const &ideal, polynomial f);// Computes a Gr\"obner basis for (I:f). The input set "ideal" does not have to be a Gr\"obner basis

 
PolynomialSet idealIntersection(PolynomialSet const &a, PolynomialSet const &b);//Computes the intersection of two ideals in the same polynomial ring


#endif
