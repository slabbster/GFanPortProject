#ifndef TRAVERSER_SPHERE_H_INCLUDED
#define TRAVERSER_SPHERE_H_INCLUDED

#include "symmetrictraversal.h"
#include "polynomial.h"

/**
 */

class SphereTraverser: public ConeTraverser
{
  int currentConeIndex;
  vector<PolyhedralCone> const &cones;
  map<IntegerVector,list<int> > adjacency;
          PolyhedralCone theCone;
        PolyhedralCone theRestrictingCone;
//        int n,d;
        void updatePolyhedralCone();
public:
  IntegerVector currentNormal;
        SphereTraverser(vector<PolyhedralCone> const &cones_, map<IntegerVector,list<int> > const &ajacency_, IntegerVector const &startCone, IntegerVector const &normal);
        virtual void changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector);
        virtual IntegerVectorList link(IntegerVector const &ridgeVector);
        PolyhedralCone & refToPolyhedralCone();
};

#endif
