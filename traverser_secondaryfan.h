#ifndef TRAVERSER_SECONDARYFAN_H_INCLUDED
#define TRAVERSER_SECONDARYFAN_H_INCLUDED

#include "symmetrictraversal.h"
#include "triangulation2.h"

Triangulation2 triangulationWithFullDimensionalIntersection(Triangulation2 g, PolyhedralCone const &c);

class SecondaryFanTraverser: public ConeTraverser
{
	Triangulation2 theTriangulation;
	PolyhedralCone theCone;
	PolyhedralCone theRestrictingCone;
        bool isHomogeneous;
	int n,d;
	void updatePolyhedralCone();
public:
	SecondaryFanTraverser(Triangulation2 const &triangulation_);
	SecondaryFanTraverser(Triangulation2 const &triangulation_, PolyhedralCone const &restrictingCone);
	virtual void changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector);
	virtual IntegerVectorList link(IntegerVector const &ridgeVector);
	PolyhedralCone & refToPolyhedralCone();
};

#endif
