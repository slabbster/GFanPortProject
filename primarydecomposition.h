#ifndef PRIMARYDECOMPOSITION_H_INCLUDED
#define PRIMARYDECOMPOSITION_H_INCLUDED

#include "polynomial.h"

class PrimaryDecompositionEngine
{
	  static class PrimaryDecompositionEngine *list;
	  class PrimaryDecompositionEngine *next;
 public:
  PrimaryDecompositionEngine();
  virtual const char *name()=0;
  static PrimaryDecompositionEngine *find(const char *name);

  virtual PolynomialSetList minimalAssociatedPrimes(PolynomialSet const &idealGenerators)=0;
};


list<PolynomialSet> minimalAssociatedPrimes(PolynomialSet const &idealGenerators);

#endif
