#ifndef NEWTONPOLYTOPE_H_INCLUDED
#define NEWTONPOLYTOPE_H_INCLUDED

#include "vektor.h"
#include "polynomial.h"
#include "halfopencone.h"

IntegerVectorList newtonPolytope(Polynomial const &p);
HalfOpenConeList normalFan(int dimension, IntegerVectorList l, TermOrder const &t);
void removeNonExtremeVerticesOfPolytope(IntegerVectorList &polytope);


#endif
