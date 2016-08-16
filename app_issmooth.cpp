#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "buchberger.h"
#include "wallideal.h"
#include "dimension.h"
#include "saturation.h"
#include "determinantpoly.h"
#include "tropical2.h"
#include "log.h"

class IsSmoothApplication : public GFanApplication
{
  SimpleOption optionInitialForms;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "Checks if the input ideal defines a smooth variety in $(C^*)^n$.\n";
  }
  IsSmoothApplication():
    optionInitialForms("--initialideals","Check the condition for the initial ideals of the ideals with respect to a given list of weight vectors. THE INPUT IDEAL MUST BE HOMOGENEOUS.")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_issmooth";
  }
  bool processSingle(PolynomialSet g)const
  {
    g=g.multiDeHomogenization();

    PolynomialRing theRing=g.getRing();
    PolynomialSet G=nonHomogeneousSaturation(g);
    log1 debug<<G;
    PolynomialSet G2=G;
    IntegerVector v=IntegerVector::allOnes(G2.getRing().getNumberOfVariables());
    buchberger(&G2,WeightReverseLexicographicTermOrder(v));
    log1 debug<<G2;
    int d=krullDimension(G2);
    PolynomialSet A=jacobiMinors(G,theRing.getNumberOfVariables()-d);
    for(PolynomialSet::const_iterator i=G.begin();i!=G.end();i++)
      A.push_back(*i);
    A.removeZeros();
    A.removeDuplicates();
 //   debug<<A;
    A=nonHomogeneousSaturation(A);
//    debug<<A;
    return (A.isUnitIdeal());
  }
  int main()
  {
    FileParser P(Stdin);
    PolynomialSet g=P.parsePolynomialSetWithRing();
    PolynomialRing theRing=g.getRing();

    if(optionInitialForms.getValue())
      {
        IntegerVectorList l=P.parseIntegerVectorList();
        bool isSchoen=true;
        for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
          {
            bool smooth=processSingle(initialIdeal(g,*i));
            if(!smooth)isSchoen=false;
            log1 debug<<smooth<<"\n";
          }
        pout<<int(isSchoen)<<"\n";
      }
    else
      pout<<int(processSingle(g))<<"\n";

    return 0;
  }
};

static IsSmoothApplication theApplication;
