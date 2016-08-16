#ifndef GROEBNERENGINE_H_INCLUDED
#define GROEBNERENGINE_H_INCLUDED

#include "polynomial.h"

class GroebnerEngine
{
  static class GroebnerEngine *list;
  class GroebnerEngine *next;
 public:
  GroebnerEngine();
  //  static LpSolver *find(const char *name);
  //  static void printList(FILE *f);
  virtual const char *name()=0;
  static GroebnerEngine *find(const char *name);

  /**
     Computes a Groebner basis of the input ideal with respect to the
     degree reverse lexicographic term order induced by the weight.
     The input generators of the ideal must be homogeneous with
     respect to a positive vector or all entries of weight must be
     strictly positive. The returned basis might not be minimal and
     might not be reduced.??????????????????
   */
  //  virtual PolynomialSet groebnerBasisDegreeRevLex(PolynomialSet const &ideal, IntegerVector const &weight)=0;
  virtual PolynomialSet groebnerBasis(bool &success, PolynomialSet const &idealGenerators, TermOrder const &termOrder, bool autoreduce)=0;
  /**
     
   */
  virtual PolynomialSet autoReduce(bool &success, PolynomialSet const &idealGenerators)=0;  
};


PolynomialSet GE_groebnerBasis(PolynomialSet const &idealGenerators, TermOrder const &termOrder, bool autoreduce);
PolynomialSet GE_autoReduce(PolynomialSet const &idealGenerators);  
bool GE_SetEngine(const char *name);


#endif
