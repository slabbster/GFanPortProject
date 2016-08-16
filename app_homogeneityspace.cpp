#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "termorder.h"
#include "gfanapplication.h"
#include "wallideal.h"
#include "polyhedralcone.h"

class HomogeneitySpaceApplication : public GFanApplication
{
public:
  const char *helpText()
  {
    return "This program computes the homogeneity space of a list of polynomials - as a cone. Thus generators for the homogeneity space are found in the section LINEALITY_SPACE. If you wish the homogeneity space of an ideal you should first compute a set of homogeneous generators and call the program on these. A reduced Groebner basis will always suffice for this purpose.\n";
  }
  HomogeneitySpaceApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_homogeneityspace";
  }

  int main()
  {
    FileParser P(Stdin);

    PolynomialSet g=P.parsePolynomialSetWithRing();

    IntegerVectorList l=wallInequalities(g);
    IntegerVectorList a;
    PolyhedralCone c(a,l);

    c.canonicalize();
    AsciiPrinter(Stdout).printPolyhedralCone(c);

    return 0;
  }
};

static HomogeneitySpaceApplication theApplication;

