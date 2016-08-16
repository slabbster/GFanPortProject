#ifndef GENERICWALK_H_INCLUDED
#define GENERICWALK_H_INCLUDED

#include "vektor.h"
#include "polynomial.h"
#include "termorder.h"

IntegerVectorList::const_iterator shootGenericRay(IntegerVectorList const &g, const TermOrder &source, const TermOrder &target);
PolynomialSet genericWalk(PolynomialSet const &start, const TermOrder &source, const TermOrder &target);
PolynomialSet genericWalkPerturbation(PolynomialSet const &start, const TermOrder &source, const TermOrder &target, int sourceDegree, int targetDegree);

#endif
