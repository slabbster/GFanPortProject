#ifndef LP_CDD_H_INCLUDED
#define LP_CDD_H_INCLUDED

#include "lp.h"

#include <cstdio> /* Always include cstdio before gmp.h.*/
#include <gmp.h>  //SHOULD BE REMOVED

IntegerVector arrayToIntegerVector(mpq_t *point, int n);
void scaleToIntegerVector(mpq_t *point, int n);

class LpSolverCdd : public LpSolver
{
 public:
  const char *name(){return "cdd";}
  bool isFacet(const IntegerVectorList &g, IntegerVectorList::const_iterator i);
};

class LpSolverCddGmp : public LpSolver
{
 public:
  virtual const char *name(){return "cddgmp";}
  bool isFacet(const IntegerVectorList &g, IntegerVectorList::const_iterator i);
  bool interiorPoint(const IntegerVectorList &g, IntegerVector &result, bool strictlyPositive, IntegerVector const *equalitySet=0);
  bool hasInteriorPoint(const IntegerVectorList &g, bool strictlyPositive, IntegerVector const *equalitySet);
  IntegerVectorList::const_iterator shoot(const IntegerVectorList &g);
  bool positiveVectorInKernel(const IntegerVectorList &g, IntegerVector *result);
  int rankOfMatrix(const IntegerVectorList &g);
  IntegerVectorList extremeRaysInequalityIndices(const IntegerVectorList &inequalityList);
  void removeRedundantRows(IntegerVectorList *inequalities, IntegerVectorList *equalities, bool removeInequalityRedundancies);
  IntegerVector relativeInteriorPoint(int n, const IntegerVectorList &g, IntegerVector const *equalitySet);
  void dual(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations, IntegerVectorList *dualInequalities, IntegerVectorList *dualEquations);//lineality space have been computed
  virtual bool hasHomogeneousSolution(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations);
};

#endif
