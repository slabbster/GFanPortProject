#ifndef TRAVERSER_STABLEINTERSECTION_H_INCLUDED
#define TRAVERSER_STABLEINTERSECTION_H_INCLUDED

#include "symmetrictraversal.h"
#include "polynomial.h"

class StableIntersectionTraverser: public ConeTraverser
{
	PolynomialSet coneGroebnerBasis;
	PolynomialSet idealGroebnerBasis;
	PolyhedralCone theCone;
	int n,d;
	void updatePolyhedralCone(int multiplicity);
public:
	StableIntersectionTraverser(PolynomialSet const &coneGroebnerBasis_, PolynomialSet const &idealGroebnerBasis_);
	virtual void changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector);
	virtual IntegerVectorList link(IntegerVector const &ridgeVector);
	virtual PolyhedralCone & refToPolyhedralCone();
};

#endif
