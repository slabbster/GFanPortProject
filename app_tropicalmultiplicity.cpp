#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "gfanapplication.h"
#include "multiplicity.h"
#include "lp.h"
#include "saturation.h"
#include "log.h"

class TropicalMultiplicityApplication : public GFanApplication
{
public:
  const char *helpText()
  {
    return "This program computes the multiplicity of a tropical cone given a marked reduced Groebner basis for its initial ideal.\n";
  }
  TropicalMultiplicityApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_tropicalmultiplicity";
  }

  int main()
  {
    PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();

    log1 AsciiPrinter(Stdout).printPolynomialSet(nonHomogeneousSaturation(g));
    AsciiPrinter(Stdout).printInteger(multiplicity(g));
    AsciiPrinter(Stdout).printNewLine();

    return 0;
  }
};

static TropicalMultiplicityApplication theApplication;

