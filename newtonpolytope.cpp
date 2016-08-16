#include "newtonpolytope.h"

#include "lp.h"

void removeNonExtremeVerticesOfPolytope(IntegerVectorList &polytope)
{
  if(polytope.empty())return;

  int n=polytope.begin()->size();

  for(IntegerVectorList::iterator j=polytope.begin();j!=polytope.end();j++)
    {
      j->resize(n+1);
      (*j)[n]=1;
    }

  for(IntegerVectorList::iterator j=polytope.begin();j!=polytope.end();j++)
    if(!isFacet(polytope,j))
      {
	IntegerVectorList::iterator k=j;
	j++;
	polytope.erase(k);
	j--;
      }
  
  for(IntegerVectorList::iterator j=polytope.begin();j!=polytope.end();j++)
    j->resize(n);
}


IntegerVectorList newtonPolytope(Polynomial const &p)
{
  IntegerVectorList polytope;

  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      polytope.push_back(i->first.exponent);      
    }

  removeNonExtremeVerticesOfPolytope(polytope);
  
  return polytope;
}


/*
  This routine should be split into several so that it can also handle closed cones and tropical hypersurfaces.
 */
HalfOpenConeList normalFan(int dimension, IntegerVectorList l, TermOrder const &t)
{
  HalfOpenConeList ret;

  removeNonExtremeVerticesOfPolytope(l);
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      IntegerVectorList strictList,nonStrictList;
      for(IntegerVectorList::const_iterator j=l.begin();j!=l.end();j++)
	if(j!=i)
	  {
	    IntegerVector v=*i-*j;
	    if(t(*i,*j))//we don't need to find the facets, right?
	      strictList.push_back(v);
	    else
	      nonStrictList.push_back(v);
	    
	    IntegerVectorList empty;
	    ret.push_back(HalfOpenCone(dimension,empty,nonStrictList,strictList,true));
	  }
    }

  return ret;
}
