#include "traverser_groebnerfan.h"
#include <iostream>
#include "division.h"
#include "wallideal.h"
#include "groebnerengine.h"
#include "tropical2.h"
#include "log.h"
#include "division.h"
#include "buchberger.h"

PolynomialSet groebnerBasisWithFullDimensionalIntersection(PolynomialSet g, PolyhedralCone const &c)
{
	IntegerVector v=c.getRelativeInteriorPoint();
	IntegerVectorList l=c.generatorsOfSpan();
	l.push_front(v);
	MatrixTermOrder T(l);
	buchberger(&g,T);
	return g;
}


GroebnerFanTraverser::GroebnerFanTraverser(PolynomialSet const &groebnerBasis_):
	groebnerBasis(groebnerBasis_),
	theCone(groebnerBasis_.getRing().getNumberOfVariables()),
	theRestrictingCone(groebnerBasis_.getRing().getNumberOfVariables()),
	isHomogeneous(false),
	isKnownToBeComplete(false)
	{
		n=groebnerBasis_.getRing().getNumberOfVariables();
		d=n;
		updatePolyhedralCone();
//cerr<<"DIMIMIMiM"<<d<<endl;
	}

GroebnerFanTraverser::GroebnerFanTraverser(PolynomialSet const &groebnerBasis_, PolyhedralCone const &restrictingCone):
	groebnerBasis(groebnerBasis_),
	theCone(groebnerBasis_.getRing().getNumberOfVariables()),
	theRestrictingCone(restrictingCone),
	isHomogeneous(false),
	isKnownToBeComplete(false)
	{
		d=theRestrictingCone.dimension();
		n=groebnerBasis_.getRing().getNumberOfVariables();
		updatePolyhedralCone();

		isHomogeneous=theCone.linealitySpace().containsPositiveVector();
		//d=n;
		//cerr<<"DIMIMIMiM"<<d<<endl;
	}

void GroebnerFanTraverser::updatePolyhedralCone()
{
//	IntegerVectorList empty;
//	theCone=PolyhedralCone::polyhedralConeWithKnownImpliedEquations(fastNormals(wallInequalities(groebnerBasis)),empty,n);
//	theCone=PolyhedralCone(fastNormals(wallInequalities(groebnerBasis)),empty,n);
//	theCone=PolyhedralCone(fastNormals(algebraicTest(wallInequalities(groebnerBasis),groebnerBasis)),empty,n);

	IntegerVectorList inequalities=fastNormals(wallInequalities(groebnerBasis));
	IntegerVectorList inequalities2=theRestrictingCone.getHalfSpaces();
	inequalities.splice(inequalities.begin(),inequalities2);
	theCone=PolyhedralCone::polyhedralConeWithKnownImpliedEquations(inequalities,theRestrictingCone.getEquations(),n);
//	theCone=intersection(theCone,theRestrictingCone);
	theCone.canonicalize();

/*	if(theCone.dimension()!=d)
	{
		cerr<<"WARNING THIS CONE HAS THE WRONG DIMENSION!"<<endl;

		cerr<<"d"<<d<<"n"<<n;

		AsciiPrinter P(Stderr);
		theRestrictingCone.print(&P);
		assert(0);
	}*/
}

void GroebnerFanTraverser::changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)
{
	assert(theCone.contains(ridgeVector));
	if(theRestrictingCone.isFullSpace())
	{
		//Notice that changeCone might be called even when no change is needed
		//This happens when we go to a ridge, but there is nothing to be done
		//In that case we return to the facet that we came from.
		//The call to changeCone should be interpreted as go to ridge, then cone back
		//to the facet in direction rayVector.
		{
			IntegerVectorList l=wallInequalities(groebnerBasis);
			bool testPerformed=false;
			for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
				if(dotLong(*i,ridgeVector)==0)
				{
					if(dotLong(rayVector,*i)>=0)return;//<<--------------------
					testPerformed=true;break;
				}
			assert(testPerformed);
		}
//		groebnerBasis=flip(groebnerBasis,-rayVector);
		{
			IntegerVectorList m;
			m.push_back(ridgeVector);
			m.push_back(rayVector);
			MatrixTermOrder T(m);
		groebnerBasis=flip(groebnerBasis,-rayVector,&T);
		}
	}
	else
	{
		PolynomialSet ridgeIdealOld=initialFormsAssumeMarked(groebnerBasis,ridgeVector);
		WeightReverseLexicographicTermOrder T(rayVector);
		PolynomialSet ridgeIdeal=GE_groebnerBasis(ridgeIdealOld,T,true);
		PolynomialSet g2(groebnerBasis.getRing());
		for(PolynomialSet::const_iterator j=ridgeIdeal.begin();j!=ridgeIdeal.end();j++)
	      g2.push_back(divisionLift(*j, ridgeIdealOld, groebnerBasis, T));
		groebnerBasis=GE_autoReduce(g2);
	}
	updatePolyhedralCone();
}

IntegerVectorList GroebnerFanTraverser::link(IntegerVector const &ridgeVector)
{
	IntegerVectorList ret;
	IntegerVector v=theCone.link(ridgeVector).getUniquePoint();

	ret.push_back(v);

	if(isKnownToBeComplete)
	  {
            ret.push_back(-v);
	    return ret;
	  }

	if(!theRestrictingCone.containsRelatively(ridgeVector))return ret;

	// If the ideal is known to be homogeneous in a positive grading then we don't have to make the following test
	// Furthermore, for restricted traversals the following test test the ridge for having a positive vector - another definitions of flipability would be to check the lifting of the ridge to the full space.
	if(isHomogeneous)
	{
//		debug<<"FDSFADFA-------------------------------\n";
		ret.push_back(-v);
		return ret;
	}
	{
	PolyhedralCone temp=theCone.faceContaining(ridgeVector);
	if(temp.containsPositiveVector())ret.push_back(-v);
	//cerr<<"THIS NEEDS TO BE FIXED!!!\n";
//	ret.push_back(-v);
	}
	return ret;
}

PolyhedralCone & GroebnerFanTraverser::refToPolyhedralCone()
{
	return theCone;
}

PolynomialSet &GroebnerFanTraverser::refToGroebnerBasisRepresentation()
{
	return groebnerBasis;
}

void GroebnerFanTraverser::setIsKnownToBeComplete(bool complete)
{
  isKnownToBeComplete=complete;
}
