#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "newtonpolytope.h"
#include "buchberger.h"
#include "wallideal.h"
#include "log.h"

class IsMarkedGroebnerBasisApplication : public GFanApplication
{
public:
  const char *helpText()
  {
    return "This program checks if a set of marked polynomials is a Groebner basis with respect to its marking. First it is checked if the markings are consistent with respect to a positive vector. Then Buchberger's S-criterion is checked. The output is boolean value.\n";
  }
  IsMarkedGroebnerBasisApplication()
  {
    registerOptions();
  }
  const char *name()
  {
    return "_ismarkedgroebnerbasis";
  }
  int main()
  {
    FileParser P(Stdin);
    PolynomialSet g=P.parsePolynomialSetWithRing();

    bool isConsistent=isMarkingConsistent(g);

    log1 fprintf(Stderr,"Is marking consistent:%i\n",isConsistent);
    bool isGroebnerBasis=isConsistent;
    if(isGroebnerBasis)
      isGroebnerBasis=isMarkedGroebnerBasis(g);
    fprintf(Stdout,isGroebnerBasis?"true\n":"false\n");

    return 0;
  }
};

static IsMarkedGroebnerBasisApplication theApplication;
