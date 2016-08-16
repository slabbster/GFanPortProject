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


class TropicalFunctionApplication : public GFanApplication
{

public:
  const char *helpText()
  {
    return "This program takes a polynomial and tropicalizes it. The output is piecewise linear function represented by a fan whose cones are the linear regions. Each ray of the fan gets the value of the tropical function assigned to it. In other words this program computes the normal fan of the Newton polytope of the input polynomial with additional information.";
  }
  TropicalFunctionApplication()
  {
    registerOptions();
  }
  const char *name()
  {
    return "_tropicalfunction";
  }
  int main()
  {
    FileParser P(Stdin);

    PolynomialSet f=P.parsePolynomialSetWithRing();
    int n=f.numberOfVariablesInRing();

    PolyhedralFan F=PolyhedralFan::normalFanOfNewtonPolytope(*f.begin());

    {
      AsciiPrinter p(Stdout);
      PolyhedralFan a=F;
      a.printWithIndices(&p,FPF_default|FPF_values);
    }

    return 0;
  }
};

static TropicalFunctionApplication theApplication;
