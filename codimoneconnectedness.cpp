#include "codimoneconnectedness.h"
#include "graph.h"
#include "log.h"

#include <set>
#include <iostream>


void CodimOneConnectednessTester::insertFacetOrbit(IntegerVectorList const &ridges)
{
  data.push_back(ridges);
}


bool CodimOneConnectednessTester::isConnected()const
{
  set<IntegerVector> allRidges;
  for(list<IntegerVectorList>::const_iterator i=data.begin();i!=data.end();i++)
    for(IntegerVectorList::const_iterator j=i->begin();j!=i->end();j++)
      allRidges.insert(*j);

  vector<IntegerVector> allRidges2;//(allRidges.size());
  for(set<IntegerVector>::const_iterator i=allRidges.begin();i!=allRidges.end();i++)allRidges2.push_back(*i);

  int nRidges=allRidges2.size();
  int nFacets=data.size();

  Graph g(nFacets+nRidges);

  int facetIndex=0;
  for(list<IntegerVectorList>::const_iterator i=data.begin();i!=data.end();i++,facetIndex++)
    for(IntegerVectorList::const_iterator j=i->begin();j!=i->end();j++)
      g.addEdge(facetIndex,nFacets+(lower_bound(allRidges2.begin(),allRidges2.end(),*j)-allRidges2.begin()));

  int diameter=g.diameter();
  log2 cerr << "Diameter " << diameter << " nFacets " << nFacets << " nRidges " << nRidges << endl;
  return diameter<nFacets+nRidges;
}
