#ifndef LP_SOPLEX_H_INCLUDED
#define LP_SOPLEX_H_INCLUDED

#include "lp.h"
#include "lp_cdd.h"

class LpSolverSoPlexCddGmp : public LpSolverCddGmp
{
 public:
  const char *name(){return "SoPlexCddGmp";}
  bool isFacet(const IntegerVectorList &g, IntegerVectorList::const_iterator i);
  bool hasHomogeneousSolution(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations);
};

#endif
