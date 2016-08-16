#include "traverser_tropical.h"
#include <iostream>
#include "bergman.h"
#include "tropical.h"
#include "division.h"
#include "wallideal.h"
#include "groebnerengine.h"
#include "tropical2.h"
#include "multiplicity.h"
#include "buchberger.h"//in time use groebnerengine.h instead
#include "log.h"

static void checkSameLeadingTerms(PolynomialSet const &a, PolynomialSet const &b)
{
  assert(a.size()==b.size());
  PolynomialSet::const_iterator A=a.begin();
  for(PolynomialSet::const_iterator B=b.begin();B!=b.end();B++,A++)
    assert(A->getMarked().m.exponent==B->getMarked().m.exponent);
}

TropicalTraverser::TropicalTraverser(PolynomialSet const &coneGroebnerBasis_, PolynomialSet const &idealGroebnerBasis_):
	coneGroebnerBasis(coneGroebnerBasis_),
	idealGroebnerBasis(idealGroebnerBasis_),
	theCone(coneGroebnerBasis_.getRing().getNumberOfVariables())
	{
		n=coneGroebnerBasis_.getRing().getNumberOfVariables();
		updatePolyhedralCone();
		PolyhedralCone homogeneitySpac=homogeneitySpace(coneGroebnerBasis);
		d=homogeneitySpac.dimensionOfLinealitySpace();
	}

void TropicalTraverser::updatePolyhedralCone()
{
	//AsciiPrinter(Stderr)<<coneGroebnerBasis;
	//fprintf(Stderr,"%i",n);
	theCone=PolyhedralCone(fastNormals(wallInequalities(idealGroebnerBasis)),wallInequalities(coneGroebnerBasis),n);
	theCone.canonicalize();
	theCone.setMultiplicity(multiplicity(coneGroebnerBasis));
}

void TropicalTraverser::changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)
{
log2 {
	debug << "Interior point:"<<theCone.getUniquePoint()<<"\n";
	debug << "Ridge:"<<ridgeVector<<"Ray:"<<rayVector<<"\n";
}
	assert(idealGroebnerBasis.containsInClosedGroebnerCone(ridgeVector));
	  log2 cerr<<endl<<"Changing cone"<<endl;

	//  assert(!containsMonomial(coneGroebnerBasis));

	  AsciiPrinter P(Stderr);
	  //P<<idealGroebnerBasis;
	  //P<<ridgeVector<<rayVector;

	  PolynomialSet ridgeIdeal=initialFormsAssumeMarked(idealGroebnerBasis,ridgeVector);
	  PolynomialSet ridgeIdealOld=ridgeIdeal;
	  WeightReverseLexicographicTermOrder T(rayVector);

	  //  P<<ridgeIdeal;
	  log2 cerr<<"Computing initial Groebner basis"<<endl;
	  //  buchberger(&ridgeIdeal,T);

	  ridgeIdeal=GE_groebnerBasis(ridgeIdeal,T,true);
	  //printMarkedTermIdeal(ridgeIdeal,"ridgeIdeal");
	  coneGroebnerBasis=initialFormsAssumeMarked(ridgeIdeal,rayVector);

	  PolynomialSet g2(coneGroebnerBasis.getRing());
	  //  WeightTermOrder termOrder(termorderWeight(ridgeIdeal));
	  WeightTermOrder termOrder(termorderWeight(ridgeIdealOld));

	  log2 cerr<<"Lifting"<<endl;
	  PolynomialSet temp=ridgeIdealOld;
	  temp.markAndScale(T);
	  temp=temp.markedTermIdeal();
	  //  P<<temp;
	  checkSameLeadingTerms(ridgeIdealOld,idealGroebnerBasis);
	  for(PolynomialSet::const_iterator j=ridgeIdeal.begin();j!=ridgeIdeal.end();j++)
	    {
	      /*      {
		Term m=j->getMarked();
		P.printPolynomial(m);
	       	if(division(m,temp,T).isZero()){cerr<<"YES";}
		cerr<<endl;
		}*/

	      g2.push_back(divisionLift(*j, ridgeIdealOld, idealGroebnerBasis, termOrder));
	     // cerr<<"*";
	    }
	  assert(g2.isMarked());
	  //printMarkedTermIdeal(g2,"g2");
	  log2 cerr<<"Autoreducing"<<endl;

	  //  autoReduce(&g2,LexicographicTermOrder());
	  //PolynomialSet g2Old=g2;
	  int oldSize=g2.size();

	  // SINGULAR DOES NOT AUTO REDUCE AT THE MOMENT!
	  //idealGroebnerBasis=GE_autoReduce(g2);
	  {
	    idealGroebnerBasis=g2;
	    IntegerVectorList M;
            M.push_back(ridgeVector);
            M.push_back(rayVector);
            MatrixTermOrder T(M);
            assert(idealGroebnerBasis.checkMarkings(T));
            autoReduce(&idealGroebnerBasis,T);
	  }

	    //printMarkedTermIdeal(g2,"g2");
	  /*  if(g2.size()!=oldSize)
	    {
	      P<<g2Old;
	      P<<g2;
	      }*/
	  assert(idealGroebnerBasis.size()==oldSize);
	  //  idealGroebnerBasis=g2;
//	  assert(!containsMonomial(coneGroebnerBasis));
	  log2 cerr<<"Done changing cone"<<endl<<endl;

	  //  P<<coneGroebnerBasis;
	  //  P<<idealGroebnerBasis;

	  /*  log0{
	    WeightReverseLexicographicTermOrder T(ridgeVector);//    WeightTermOrder
	    PolynomialSet A=idealGroebnerBasis;
	    buchberger(&A,T);
	    //cerr<<A;
	    PolynomialSet B=initialFormsAssumeMarked(A,ridgeVector);
	    WeightReverseLexicographicTermOrder T2(rayVector);
	    buchberger(&B,T2);
	    cerr<<"RIGHT";
	    P<<initialFormsAssumeMarked(B,rayVector);
	    cerr<<"RIGHT?";
	    P<<coneGroebnerBasis;
	    }*/

	  log2 cerr << "Number of terms in new basis: "<< g2.totalNumberOfTerms()<<endl;


	updatePolyhedralCone();
}

IntegerVectorList TropicalTraverser::link(IntegerVector const &ridgeVector)
{
	assert(idealGroebnerBasis.containsInClosedGroebnerCone(ridgeVector));
	PolynomialSet tempIdeal=initialFormsAssumeMarked(idealGroebnerBasis,ridgeVector);

	//P<<tempIdeal;
	tempIdeal=saturatedIdeal(tempIdeal);//TODO: figure out if it is an advantage to saturate the ideal
	IntegerVectorList rays;

	PolyhedralCone saturatedHomogeneitySpace=homogeneitySpace(tempIdeal);
	if(saturatedHomogeneitySpace.dimensionOfLinealitySpace()==d)//if saturation of the initial ideal changed the homogeneity space things are easy
	  {
		IntegerVector v=theCone.link(ridgeVector).getUniquePoint();
		rays.push_back(v);
		rays.push_back(-v);
//	    rays.push_back(top.parentRay);
//	    rays.push_back(-top.parentRay);
	  }
	else
	  {
	    //P<<tempIdeal;
	    BergmanFan b=bergmanRayIntersection(tempIdeal);
	    //	  BergmanFan b=bergmanRayIntersection(initialIdeal(idealGroebnerBasis,facetStack.front().ridges.front()));


	    bool trouble=false;
	    for(BergmanFan::MaximalConeList::const_iterator i=b.cones.begin();i!=b.cones.end();i++)
		{
		  PolyhedralCone rayCone=i->theCone;
		  rayCone.canonicalize();
		  {
		    if(rayCone.getUniquePoint().isZero())trouble=true;
		  }
		  rays.push_back(rayCone.getUniquePoint());
		}
	    if(trouble)
		{
//		  b.print(P);
//		  P<<tempIdeal;
		  assert(!trouble);
	    }
	  }

	return rays;
}

PolyhedralCone & TropicalTraverser::refToPolyhedralCone()
{
	return theCone;
}

