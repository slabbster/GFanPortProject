#include "polyhedralcone.h"
#include "termorder.h"
#include "buchberger.h"
#include "wallideal.h"
#include "tropical2.h"

class RestrictedGFanEnumeration{
  RestrictedGFanEnumeration(PolynomialSet const &I, PolyhedralCone const &c)
  {
    int n=c.ambientDimension();
    IntegerVectorList inequalities=c.getHalfSpaces();
    IntegerVectorList equations=c.getEquations();
    IntegerVector v=c.getRelativeInteriorPoint(); //positive orthant

    WeightReverseLexicographicTermOrder myOrder(v);

    PolynomialSet g=I;
    buchberger(&g,myOrder);
    //    minimize(&g);
    //    autoReduce(&g,myOrder);

    IntegerVectorList facets=wallInequalities(g);
    facets=algebraicTest(facets,g);

    inequalities.splice(inequalities.end(),facets);

    PolyhedralCone c2(inequalities,equations,n);

    c2.findFacets();

    IntegerVectorList facets2=c2.getHalfSpaces();
    for(IntegerVectorList::const_iterator i=facets2.begin();i!=facets2.end();i++)
      {
	IntegerVectorList equations2=equations;
	equations2.push_back(*i);
	PolyhedralCone c3(inequalities,equations2,n);
	IntegerVector v3=c3.getRelativeInteriorPoint();
	PolynomialSet g3=initialFormsAssumeMarked(g,v3);

	IntegerVectorList inequalities2=inequalities;
	inequalities.push_back(-(*i));
	IntegerVectorList empty;
	PolyhedralCone c4(inequalities2,empty,n);
	IntegerVector v4=c4.getRelativeInteriorPoint();
	WeightReverseLexicographicTermOrder myOrder4(v4);
	PolynomialSet g4=g3;
	buchberger(&g4);



	//Lift g to a basis compatible with termorder of g4
	PolynomialSet g5(theRing);
	for(PolynomialSet::const_iterator j=g4.begin();j!=g4.end();j++)
	  g5.push_back(divisionLift(*j, g3, g, LexicographicTermOrder()));
	autoreduce(g5);

	PolynomialSet g6=initialFormsAssumeMarked(g4,v4);
	//We now have a pair g5 and g6
	//We don't really need g6, do we?
      }
  }
};
