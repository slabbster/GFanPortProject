#include "mixedvolume.h"

#include "tropical_weildivisor.h"
#include "multiplicity.h"

int64 mixedVolume(PolynomialSet const &g_)
{
  PolynomialSet g=idealWithSameMultiplicity(g_);
  g.simplestPolynomialsFirst();
  int n=g.getRing().getNumberOfVariables();

	assert(g.size()<=n);

	PolyhedralFan f=PolyhedralFan::fullSpace(n);

	for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
	{
	    PolyhedralFan f2=PolyhedralFan::normalFanOfNewtonPolytope(*i);
	    if(f.size()==0)return 0;
	    f=weilDivisor(f,f2);
	}
	if(f.size()==0)return 0;
	assert(f.size()==1);
	return f.conesBegin()->getMultiplicity();
}

#include "halfopencone.h"
bool mixedVolumePositive(PolynomialSet const &g)
{
  return nonEmptyStableIntersection(g);
}
