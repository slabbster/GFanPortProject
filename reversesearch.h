#ifndef REVERSESEARCH_H_INCLUDED
#define REVERSESEARCH_H_INCLUDED

#include "enumeration.h"
#include "termorder.h"

class ReverseSearch: public EnumerationAlgorithm
{
  int numberOfEdges;
  int numberOfVertices;
  bool isKnownToBeHomogeneous;
  bool broken;
  int treeSize(PolynomialSet &groebnerBasis);
  const TermOrder &termOrder;
 protected:
  bool computeSearchEdge(PolynomialSet &groebnerBasis, IntegerVector *edge);
 public:
  PolynomialSet findRoot(PolynomialSet groebnerBasis);
  ReverseSearch(const TermOrder &termOrder_);
  void enumerate(const PolynomialSet &groebnerBasis);
};


#endif
