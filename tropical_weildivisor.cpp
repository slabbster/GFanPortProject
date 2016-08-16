#include "tropical_weildivisor.h"
#include "tropical2.h"
#include "log.h"
#include "printer.h"
#include <iostream>

PolyhedralFan weilDivisor(PolyhedralFan const &F, PolyhedralFan const &G)//, Polynomial const &g)
{
  //  PolynomialRing R=g.getRing();
  //  PolyhedralFan G=PolyhedralFan::bergmanOfPrincipalIdeal(g);
  int n=G.getAmbientDimension();
  int d=F.getMaxDimension();

  PolyhedralFan retTemp(n);

  log1 cerr<<"Computing refinement"<<endl;
  for(PolyhedralConeList::const_iterator i=F.conesBegin();i!=F.conesEnd();i++)
    {
      log1 cerr<<"*";

#if 1
      bool found=false;
      {
        IntegerVector v=i->getRelativeInteriorPoint();
        PolyhedralCone c=G.coneContaining(v);
        if(!(c!=*i))
          {
            retTemp.insert(c);
            found=true;
          }
      }
      if(!found)
#endif
      for(PolyhedralConeList::const_iterator j=G.conesBegin();j!=G.conesEnd();j++)
	{
	  PolyhedralCone c=intersection(*i,*j);
	  c.canonicalize();
	  retTemp.insert(c);
	}
    }

  log1 cerr<<"Computing full complex"<<endl;
  retTemp=retTemp.fullComplex();
  PolyhedralFan ret(n);

  log1 cerr<<"Computing divisor"<<endl;
  for(PolyhedralConeList::iterator i=retTemp.conesBegin();i!=retTemp.conesEnd();i++)
    {
      log1 cerr<<"*";
      if(i->dimension()==d-1)
	{
	  IntegerVector v=i->getRelativeInteriorPoint();

	  AsciiPrinter P(Stderr);

	  log2 P<<v<<v<<"\n";

	  int multiplicity=0;
	  IntegerVector evaluationVector(n);

	  //	  Polynomial localg=initialForm(g,v);
	  //	  PolyhedralFan localG=PolyhedralFan::normalFanOfNewtonPolytope(localg);
	  PolyhedralFan localG=G.link(v);

	  for(PolyhedralConeList::const_iterator j=F.conesBegin();j!=F.conesEnd();j++)
	    {
	      if(j->contains(v))
		{
		  IntegerVectorList equations=j->getEquations();
		  IntegerVectorList inequalities1=j->getHalfSpaces();
		  IntegerVectorList inequalities2;
		  for(IntegerVectorList::const_iterator i=inequalities1.begin();i!=inequalities1.end();i++)
		    if(dotLong(v,*i)==0)inequalities2.push_back(*i);
		  PolyhedralCone localJ(inequalities2,equations,n);
		  localJ.canonicalize();

		  PolyhedralFan refinement(n);
		  for(PolyhedralConeList::const_iterator k=localG.conesBegin();k!=localG.conesEnd();k++)
		    {
		      PolyhedralCone sigma=intersection(localJ,*k);
		      sigma.canonicalize();
		      refinement.insert(sigma);
		    }

		  for(PolyhedralConeList::const_iterator k=refinement.conesBegin();k!=refinement.conesEnd();k++)
		    {
		      PolyhedralCone const &sigma(*k);
		      if(sigma.dimension()==d)
			{
			  /*IntegerVectorList rays=sigma.extremeRays();
			  assert(rays.size()==1);//SHOULD ALWAYS BE TRUE????
			  {
			    IntegerVector ray=*rays.begin();
			*/
			  IntegerVector ray=sigma.semiGroupGeneratorOfRay();
			  evaluationVector+=j->getMultiplicity()*ray;
			  //multiplicity+=j->getMultiplicity()*localg.degree(ray);
			  multiplicity+=j->getMultiplicity()*localG.evaluatePiecewiseLinearFunction(ray);
			}
		    }

		}
	    }
	  //multiplicity-=localg.degree(evaluationVector);
	  multiplicity-=localG.evaluatePiecewiseLinearFunction(evaluationVector);
	  if(multiplicity!=0)
	    {
	      PolyhedralCone c=*i;
	      c.setMultiplicity(multiplicity);
	      ret.insert(c);
	    }
	}
    }

  return ret;
}
