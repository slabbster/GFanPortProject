#ifndef LP_H_INCLUDED
#define LP_H_INCLUDED

#include <stdio.h>

#include "vektor.h"

bool isFacet(const IntegerVectorList &g, IntegerVectorList::const_iterator i);
bool interiorPoint(const IntegerVectorList &g, IntegerVector &result, bool strictlyPositive, IntegerVector const *equalitySet=0);//The program terminates if the interior point could not be cast to an IntegerVector
bool hasInteriorPoint(const IntegerVectorList &g, bool strictlyPositive, IntegerVector const *equalitySet=0);
IntegerVectorList::const_iterator shootRay(const IntegerVectorList &g);
bool positiveVectorInKernel(const IntegerVectorList &g, IntegerVector *result);
int rankOfMatrix(const IntegerVectorList &g);
bool lpSetSolver(const char *name);
IntegerVectorList extremeRaysInequalityIndices(const IntegerVectorList &inequalityList);
void removeRedundantRows(IntegerVectorList *inequalities, IntegerVectorList *equalities, bool removeInequalityRedundancies);
IntegerVector relativeInteriorPoint(int n, const IntegerVectorList &g, IntegerVector const *equalitySet);
void dual(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations, IntegerVectorList *destInequalities, IntegerVectorList *destEquations);
bool hasHomogeneousSolution(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations);

/**
   Uses the hasHomogeneousSolution function to check if v is a linear
   combination of the rays and the linealitySpace where the
   coeefficents of the rays are forced to be non-negative.
 */
bool isInNonNegativeSpan(IntegerVector const &v, IntegerVectorList const &rays, IntegerVectorList const &linealitySpace);

class LpSolver
{
  static class LpSolver *list;
  class LpSolver *next;
 public:
  LpSolver();
  static LpSolver *find(const char *name);
  static void printList(FILE *f);
  virtual bool isFacet(const IntegerVectorList &g, IntegerVectorList::const_iterator i)=0;
  virtual bool interiorPoint(const IntegerVectorList &g, IntegerVector &result, bool strictlyPositive, IntegerVector const *equalitySet=0);
  virtual bool hasInteriorPoint(const IntegerVectorList &g, bool strictlyPositive, IntegerVector const *equalitySet=0);
  virtual const char *name()=0;
  virtual IntegerVectorList::const_iterator shoot(const IntegerVectorList &g);
  virtual bool positiveVectorInKernel(const IntegerVectorList &g, IntegerVector *result);
  virtual int rankOfMatrix(const IntegerVectorList &g);// in this header file because we use cdd for this
  virtual IntegerVectorList extremeRaysInequalityIndices(const IntegerVectorList &inequalityList);
  virtual void removeRedundantRows(IntegerVectorList *inequalities, IntegerVectorList *equalities, bool removeInequalityRedundancies);
  virtual IntegerVector relativeInteriorPoint(int n, const IntegerVectorList &g, IntegerVector const *equalitySet);
  virtual void dual(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations, IntegerVectorList *dualInequalities, IntegerVectorList *dualEquations);
  virtual bool hasHomogeneousSolution(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations);
  /* hasHomogeneousSolution()
    Description:
     n is the dimension of all vectors in the lists.
     Let A be the matrix with rows being vectors listed in "inequalities".
     This routine returns true iff there exists a vector x in R^n with Ax>=0 and with first coordinate =1 (eqivalently >0).
     If "equations" is non-empty its vectors are added into the inequality systems as rows with equality.
    About the choices for the parameter interface of this function:
     signs and order of coordinates were chosen to be compatible with cddlib.
  */
};


#endif
