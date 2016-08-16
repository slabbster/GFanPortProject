#include "vektor.h"

class CodimOneConnectednessTester
{
  list<IntegerVectorList> data;
 public:
  /**
     Adds a facet orbit. The list of vectors is the canonoical
     symmetry (w.r.t. symmetries of the ridge) invariant interior
     points of its ridges. Futhermore, these should be unique orbit
     representatives rather than just vectors.
   */
  void insertFacetOrbit(IntegerVectorList const &ridges);
  /*
    BE CAREFUL WHEN MAKING CONCLUSIONS - THIS ROUINE DOES NOT SUPPORT SYMMETRY.
   */
  bool isConnected()const;
};
