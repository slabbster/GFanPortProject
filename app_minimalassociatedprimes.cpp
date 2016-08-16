#include "parser.h"
#include "printer.h"
#include "primarydecomposition.h"
#include "gfanapplication.h"

class MinimalAssociatedPrimesApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program computes the minimal associated primes of an ideal. It attempts to call Singular through a Sage interface. The only reason for having this program is to illustrate how such communication can be done.\n";
  }
  MinimalAssociatedPrimesApplication()  {
    registerOptions();
  }

  const char *name()
  {
    return "_minimalassociatedprimes";
  }

  int main()
  {
    FileParser P(Stdin);

    PolynomialSet a=P.parsePolynomialSetWithRing();

    AsciiPrinter(Stdout).printPolynomialSetList(minimalAssociatedPrimes(a));

    return 0;
  }
};

static MinimalAssociatedPrimesApplication theApplication;
