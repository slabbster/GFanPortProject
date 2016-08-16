#include "primarydecomposition.h"

PrimaryDecompositionEngine *PrimaryDecompositionEngine::list;

static bool initialized;

PrimaryDecompositionEngine::PrimaryDecompositionEngine()
{
   next=list;
   list=this;
}


PrimaryDecompositionEngine *PrimaryDecompositionEngine::find(const char *name)
{
   PrimaryDecompositionEngine *l=list;
   while(l)
      {
         if(std::string(l->name())==std::string(name))break;
         l=l->next;
      }
   return l;
}


PolynomialSetList minimalAssociatedPrimes(PolynomialSet const &idealGenerators)
{
	PrimaryDecompositionEngine *p=PrimaryDecompositionEngine::find("sagesingular");
	assert(p);
	return p->minimalAssociatedPrimes(idealGenerators);
}
