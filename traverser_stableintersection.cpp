#include "traverser_stableintersection.h"
#include <iostream>
#include "halfopencone.h"
#include "tropical.h"
#include "division.h"
#include "wallideal.h"
#include "groebnerengine.h"
#include "tropical2.h"
#include "multiplicity.h"
#include "mixedvolume.h"
#include "log.h"

StableIntersectionTraverser::StableIntersectionTraverser(PolynomialSet const &coneGroebnerBasis_, PolynomialSet const &idealGroebnerBasis_):
	coneGroebnerBasis(coneGroebnerBasis_),
	idealGroebnerBasis(idealGroebnerBasis_),
	theCone(coneGroebnerBasis_.getRing().getNumberOfVariables())
	{
		n=coneGroebnerBasis_.getRing().getNumberOfVariables();
		updatePolyhedralCone(mixedVolume(coneGroebnerBasis));
		PolyhedralCone homogeneitySpac=homogeneitySpace(coneGroebnerBasis);
		d=homogeneitySpac.dimensionOfLinealitySpace();
	}

void StableIntersectionTraverser::updatePolyhedralCone(int multiplicity)
{
	theCone=PolyhedralCone(fastNormals(wallInequalities(idealGroebnerBasis)),wallInequalities(coneGroebnerBasis),n);
	theCone.canonicalize();
	theCone.setMultiplicity(multiplicity);
}

static void checkSameLeadingTerms(PolynomialSet const &a, PolynomialSet const &b)
{
  assert(a.size()==b.size());
  PolynomialSet::const_iterator A=a.begin();
  for(PolynomialSet::const_iterator B=b.begin();B!=b.end();B++,A++)
    assert(A->getMarked().m.exponent==B->getMarked().m.exponent);
}

void StableIntersectionTraverser::changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)
{
	IntegerVectorList m;
	m.push_back(ridgeVector);
	m.push_back(rayVector);
	MatrixTermOrder T(m);

	idealGroebnerBasis.markAndScale(T);
	coneGroebnerBasis=initialFormsAssumeMarked(initialFormsAssumeMarked(idealGroebnerBasis,ridgeVector),rayVector);

	updatePolyhedralCone(mixedVolume(coneGroebnerBasis));
	checkSameLeadingTerms(idealGroebnerBasis,coneGroebnerBasis);
}

IntegerVectorList StableIntersectionTraverser::link(IntegerVector const &ridgeVector)
{
	assert(idealGroebnerBasis.containsInClosedGroebnerCone(ridgeVector));
	PolynomialSet tempIdeal=initialFormsAssumeMarked(idealGroebnerBasis,ridgeVector);

	tempIdeal.simplestPolynomialsFirst();

	PolyhedralFan theLink=tropicalHyperSurfaceIntersectionClosed(n,tempIdeal);

	IntegerVectorList rays1=theLink.getRaysInPrintingOrder(0);

	log2 {
	  cerr<<"Ray candidates:"<<endl;
	  AsciiPrinter(Stderr)<<rays1;
	}

	IntegerVectorList rays2;
	for(IntegerVectorList::const_iterator i=rays1.begin();i!=rays1.end();i++)
		if(mixedVolumePositive(initialForms(tempIdeal,*i)))
			{
			PolyhedralCone theRidge=theCone.faceContaining(ridgeVector);
                        IntegerVectorList linealityGenerators=theRidge.generatorsOfLinealitySpace();
//                        IntegerVectorList linealityGenerators=theRidge.linealitySpace().dualCone().getEquations();
			IntegerVectorList generators;generators.push_back(*i);
			PolyhedralCone ray=PolyhedralCone::givenByRays(generators,linealityGenerators,n);
			ray.canonicalize();
			rays2.push_back(ray.getUniquePoint());
			//			rays2.push_back(*i);
			}
	log2{
	cerr<<"True rays:"<<endl;
	AsciiPrinter(Stderr)<<rays2;
	}

	return rays2;
}

PolyhedralCone & StableIntersectionTraverser::refToPolyhedralCone()
{
	return theCone;
}

