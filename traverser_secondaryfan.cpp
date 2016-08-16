#include "traverser_secondaryfan.h"
#include <iostream>
#include "wallideal.h"
#include "tropical2.h"
#include "log.h"

Triangulation2 triangulationWithFullDimensionalIntersection(Triangulation2 g, PolyhedralCone const &c)
{
	debug<<"A\n";
	PolyhedralCone A=intersection(c,g.secondaryFanSupport());
	debug<<"B\n";
	IntegerVector v=A.getRelativeInteriorPoint();
	IntegerVectorList l=c.generatorsOfSpan();
	l.push_front(v);
	MatrixTermOrder T(l);
debug<<T;
	g.changeToTriangulationInducedBy(T);
	return g;
}

SecondaryFanTraverser::SecondaryFanTraverser(Triangulation2 const &triangulation_):
	theTriangulation(triangulation_),
	theCone(triangulation_.getN()),
	theRestrictingCone(triangulation_.getN())
	{
	  isHomogeneous=theTriangulation.isHomogeneous();
	  //	debug<<"SUPPORT"<<theTriangulation.secondaryFanSupport();


	  theRestrictingCone=intersection(theRestrictingCone,theTriangulation.secondaryFanSupport());
		//debug<<"SUPPORT"<<theTriangulation.secondaryFanSupport();
		//debug<<"RESTRICTING"<<theRestrictingCone;
		n=triangulation_.getN();
//		d=n;
		d=theRestrictingCone.dimension();
		updatePolyhedralCone();
	}

SecondaryFanTraverser::SecondaryFanTraverser(Triangulation2 const &triangulation_, PolyhedralCone const &restrictingCone):
	theTriangulation(triangulation_),
	theCone(triangulation_.getN()),
	theRestrictingCone(restrictingCone)
	{
          isHomogeneous=theTriangulation.isHomogeneous();
	theRestrictingCone=intersection(theRestrictingCone,theTriangulation.secondaryFanSupport());
		d=theRestrictingCone.dimension();
		n=triangulation_.getN();
		updatePolyhedralCone();
	}


void SecondaryFanTraverser::updatePolyhedralCone()
{
	IntegerVectorList empty;

	theCone=theTriangulation.secondaryCone().negated();
	//debug<<"THECONE"<<theCone;

	theCone=intersection(theCone,theRestrictingCone);
	theCone.canonicalize();

	//debug<<theCone;
if(theCone.dimension()!=d)
	{
		cerr<<"WARNING THIS CONE HAS THE WRONG DIMENSION!"<<endl;

		cerr<<"d"<<d<<"n"<<n;

		AsciiPrinter P(Stderr);
		theRestrictingCone.print(&P);
		assert(0);
	}
}

void SecondaryFanTraverser::changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)
{
//  debug<<"changeCone\n";
  assert(theCone.contains(ridgeVector));
  //      debug<<"ridgeVector"<<ridgeVector<<"\n";
  //      debug<<"rayVector"<<rayVector<<"\n";
  if(d==n)
    {
      /* If we are traversing a full-dimensional fan, then the rayVector is
       * the unique (up to scaling) normal vector of the lineality space.*/
      theTriangulation.flipNew(-rayVector);
    }
  else
    {
	IntegerVectorList m;
	m.push_back(ridgeVector);
	m.push_back(rayVector);
	m.push_back(IntegerVector::allOnes(n));//Tie break with degree order - this might not be necessary.
	MatrixTermOrder T(m);
	theTriangulation.changeToTriangulationInducedBy(T);
    }
  updatePolyhedralCone();
//	debug<<"END changeCone\n";
}

IntegerVectorList SecondaryFanTraverser::link(IntegerVector const &ridgeVector)
{
	IntegerVectorList ret;
	IntegerVector v=theCone.link(ridgeVector).getUniquePoint();

	ret.push_back(v);

	if(!theRestrictingCone.containsRelatively(ridgeVector))return ret;


if(isHomogeneous)
	  {
	        ret.push_back(-v);
	        return ret;
	  }
  // If the ideal is known to be homogeneous in a positive grading then we don't have to make the following test
	// Furthermore, for restricted traversals the following test test the ridge for having a positive vector - another definitions of flipability would be to check the lifting of the ridge to the full space.
	PolyhedralCone temp=theCone.faceContaining(ridgeVector);
	if(temp.containsPositiveVector())ret.push_back(-v);
//	cerr<<"THIS NEEDS TO BE FIXED!!!\n";
//	ret.push_back(-v);
	return ret;
}

PolyhedralCone & SecondaryFanTraverser::refToPolyhedralCone()
{
	return theCone;
}
