#include "parser.h"
#include "printer.h"
#include "lp.h"
#include "gfanapplication.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "polymakefile.h"
#include "tropicalmap.h"

class TropicalImageApplication : public GFanApplication
{
  StringOption inputOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program computes the image of the tropicalization of a polynomial map. The output is a polyhedral fan with support equal to the image. The input is the polynomial ring, followed by a list of coordinate polynomials. A domain different from $R^n$ can be chosen with the option -i.\n";
  }
  TropicalImageApplication():
    inputOption("-i","Specify the name of the file containing a polyhedral fan whose support is the domain of the function.",0)
  {
    registerOptions();
  }

  const char *name()
  {
    return "_tropicalimage";
  }

  int main()
  {
    PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();

    int n=g.getRing().getNumberOfVariables();



    PolyhedralFan domain=PolyhedralFan::fullSpace(n);

    if(inputOption.getValue())
      {
	domain=PolyhedralFan::readFan(inputOption.getValue());
      }

    PolyhedralFan f=imageOfTropicalMap(g,domain);

    AsciiPrinter P(Stdout);

    f.printWithIndices(&P,FPF_default|FPF_multiplicities/*|FPF_values*/);

    return 0;
  }
};

static TropicalImageApplication theApplication;

