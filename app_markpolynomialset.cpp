#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "minkowskisum.h"
#include "newtonpolytope.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"

class MarkPolynomialSetApplication : public GFanApplication
{
public:
  const char *helpText()
  {
    return "This program marks a set of polynomials with respect to the vector given at the end of the input, meaning that the largest terms are moved to the front. In case of a tie the lexicographic term order with $a>b>c...$ is used to break it.\n";
  }
  MarkPolynomialSetApplication()
  {
    registerOptions();
  }
  const char *name()
  {
    return "_markpolynomialset";
  }
  int main()
  {
    FileParser P(Stdin);
    PolynomialSet g=P.parsePolynomialSetWithRing();
    IntegerVector v=P.parseIntegerVector();
    WeightTermOrder t(v);
    g.markAndScale(t);
    AsciiPrinter(Stdout).printPolynomialRing(g.getRing());
    AsciiPrinter(Stdout).printPolynomialSet(g);

    return 0;
  }
};

static MarkPolynomialSetApplication theApplication;
