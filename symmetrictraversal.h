#ifndef SYMMETRICTRAVERSAL_H_INCLUDED
#define SYMMETRICTRAVERSAL_H_INCLUDED

#include "symmetriccomplex.h"
#include "polyhedralfan.h"

/*
 This file contains the generic algorithm for traversing a connected component of a pure fan up to symmetry.
 This will in time be the algorithm to use for all fan traversals which are not reverse search.
 */
class ConeTraverser
{
public:
	/**
	 * Go to the cone which is connected to the current facet through the ridge in direction ray.
	 * The "ridge" is a relative interior point of the ridge.
	 */
	virtual void changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)=0;
/**
 * Compute the link of the fan in the ridge given by the vector ridge IS THIS A FACET NORMAL OR AN INTERIOR POINT?
 * This gives a list of symmetry invariant points under the actions keeping the link fixed.
 */
	virtual IntegerVectorList link(IntegerVector const &ridgeVector)=0;
	virtual PolyhedralCone & refToPolyhedralCone()=0;
};

class SymmetricTarget
{
public:
//	virtual bool process(PolyhedralCone const &cone)=0;
	virtual bool process(ConeTraverser &traverser)=0;
};

class SymmetricTargetFanBuilder : public SymmetricTarget
{
	PolyhedralFan coneCollection;
public:
	PolyhedralFan const &getFanRef(){return coneCollection;}
//	SymmetricComplex toSymmetricComplex()const;
	SymmetricTargetFanBuilder(int n, SymmetryGroup const &sym);
	/* return false to exit */
//	bool process(PolyhedralCone const &cone);
	bool process(ConeTraverser &traverser);
};

void symmetricTraverse(ConeTraverser &traverser, SymmetricTarget &target, SymmetryGroup const *sym=0);

#endif
