#include "tropicalmap.h"

#include <iostream>

using namespace std;

#include "log.h"

static IntegerVectorList vectorImages(PolynomialSet const &polynomialMap, IntegerVectorList const &l)
{
  int d=polynomialMap.size();
  IntegerVectorList ret;
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      IntegerVector v(d);
      int j=0;
      for(PolynomialSet::const_iterator J=polynomialMap.begin();J!=polynomialMap.end();J++,j++)
	v[j]=J->degree(*i);
      ret.push_back(v);
    }
  return ret;
}


PolyhedralCone imageOfConeUnderTropicalMap(PolynomialSet const &polynomialMap, PolyhedralCone const &cone)
{
	  int d=polynomialMap.size();
    IntegerVectorList linealitySpace=cone.linealitySpace().dualCone().getEquations();
    IntegerVectorList rays=cone.extremeRays();

    PolyhedralCone C=PolyhedralCone::givenByRays(vectorImages(polynomialMap,rays),vectorImages(polynomialMap,linealitySpace),d);
    C.canonicalize();
    return C;
}


PolyhedralFan imageOfTropicalMap(PolynomialSet const &polynomialMap, PolyhedralFan domain)
{
  for(PolynomialSet::const_iterator i=polynomialMap.begin();i!=polynomialMap.end();i++)
    {
      domain=refinement(domain,PolyhedralFan::normalFanOfNewtonPolytope(*i));
    }
  return imageOfTropicalMapGivenLinearityCones(polynomialMap,domain);
}

PolyhedralFan imageOfTropicalMapGivenLinearityCones(PolynomialSet const &polynomialMap, PolyhedralFan linearityCones)
{
	  int d=polynomialMap.size();
	  PolyhedralFan ret(d);
  /*
    We first put the images of the cones of the domain into the
    imageCones set. This automatically removes duplicates.
   */
  set<PolyhedralCone> imageCones;
  for(PolyhedralFan::coneIterator i=linearityCones.conesBegin();i!=linearityCones.conesEnd();i++)
	  imageCones.insert(imageOfConeUnderTropicalMap(polynomialMap,*i));
  log0 cerr << "Number of image cones "<<imageCones.size()<<endl;

  /*
    We then remove any cone which is contained in another.
   */
  for(set<PolyhedralCone>::iterator i=imageCones.begin();i!=imageCones.end();)
    {
      bool isContainedInOther=false;
      IntegerVectorList linealitySpace=i->linealitySpace().dualCone().getEquations();
      IntegerVectorList rays=i->extremeRays();
      for(set<PolyhedralCone>::iterator j=imageCones.begin();j!=imageCones.end();j++)
	{
	  if(j!=i)
	    if(j->contains(rays) && j->contains(linealitySpace))
	      {
		isContainedInOther=true;break;
	      }
	}
      if(isContainedInOther)
	{
	  set<PolyhedralCone>::iterator k=i;
	  i++;
	  imageCones.erase(k);
	}
      else
	i++;
    }
  log0 cerr << "Number of image cones "<<imageCones.size()<<endl;
  int I=0;
  for(set<PolyhedralCone>::const_iterator i=imageCones.begin();i!=imageCones.end();i++,I++)
    {
     //       if(I<31)continue;
     // if(I==32)break;
      //      AsciiPrinter P(Stderr);
      //      P<<C;



      ret.merge(*i);
      //ret.insert(*i);

      cerr<<"After merge"<< ret.size()<<endl;
      ret.removeNonMaximal();
      cerr<<"After simplification"<< ret.size()<<endl;
    }
  return ret;
}
