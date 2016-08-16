#include "groebnerengine.h"
#include "buchberger.h"
#include "printer.h"

class GroebnerEngineGfan : public GroebnerEngine
{
  virtual PolynomialSet groebnerBasis(bool &success, PolynomialSet const &idealGenerators, TermOrder const &termOrder, bool autoreduce)
  {
    success=true;
    PolynomialSet ret=idealGenerators;
    buchberger(&ret,termOrder);
    return ret;
  }
  virtual PolynomialSet autoReduce(bool &success, PolynomialSet const &idealGenerators)
  {
    success=true;
    PolynomialSet ret=idealGenerators;
    autoReduce_(&ret,LexicographicTermOrder());//term order is ignored in autoReduce_() but this could change
    return ret;
  }
  virtual const char* name()
  {
    return "gfan";
  }
};

static GroebnerEngineGfan groebnerEngineGfan;
