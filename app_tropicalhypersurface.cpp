#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "termorder.h"
#include "ep_standard.h"
#include "ep_xfig.h"
#include "gfanapplication.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "tropical.h"
#include "tropical2.h"
#include "symmetry.h"
#include "halfopencone.h"
#include "log.h"


class TropicalHypersurfaceApplication : public GFanApplication
{

public:
  const char *helpText()
  {
    return "This program computes the tropical hypersurface defined by a principal"
      " ideal. The input is the polynomial ring followed by a set containing"
      " just a generator of the ideal.";
  }
  TropicalHypersurfaceApplication()
  {
    registerOptions();
  }
  const char *name()
  {
    return "_tropicalhypersurface";
  }
  int main()
  {
    FileParser P(Stdin);

    PolynomialSet f=P.parsePolynomialSetWithRing();
    int n=f.numberOfVariablesInRing();

    assert(f.size()==1);

    PolyhedralFan F=PolyhedralFan::bergmanOfPrincipalIdeal(*f.begin());

    {
      AsciiPrinter p(Stdout);
      PolyhedralFan a=F;
      a.printWithIndices(&p,FPF_default|FPF_values|FPF_multiplicities);
    }

    return 0;
  }
};

static TropicalHypersurfaceApplication theApplication;
