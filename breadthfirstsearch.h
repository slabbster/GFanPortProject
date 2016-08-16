#ifndef BREADTHFIRSTSEARCH_H_INCLUDED
#define BREADTHFIRSTSEARCH_H_INCLUDED

#include "enumeration.h"
#include "termorder.h"
#include "symmetry.h"
#include "polyhedralcone.h"

PolyhedralCone coneFromMarkedBasis(PolynomialSet const &g);

class BreadthFirstSearch: public EnumerationAlgorithm
{
  bool minkowski;
  int numberOfEdges;
  int numberOfVertices;
  const SymmetryGroup &symmetryGroup;
  IntegerVectorList subspacePerp;
  /*  IntegerVectorList blocks;
  void markBoundary(IntegerVectorList &blocks, PolynomialSet const &g);
  IntegerVectorList computeLocalBlocks(PolyhedralCone const &c, bool copy);
  */
 public:
  BreadthFirstSearch(const SymmetryGroup &symmetryGroup_, bool minkowski_=false);
  void enumerate(const PolynomialSet &groebnerBasis);
  void setSubspace(IntegerVectorList const &subspacePerp_);
};


#endif
