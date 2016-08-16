#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "wallideal.h"
#include "lp.h"
#include "polyhedralcone.h"
#include "gfanapplication.h"
#include "polyhedralfan.h"
#include "halfopencone.h"
#include "matrix.h"
#include "regularsubdivision.h"
#include "log.h"

class RegularSubdivisionApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program takes a point configuration and a lifting vector a computes the corresponding regular subdivision PROJECTIVELY?.\n";
  }
  RegularSubdivisionApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_regularsubdivision";
  }

  int main()
  {
    FileParser P(Stdin);

    IntegerMatrix m=rowsToIntegerMatrix(P.parseIntegerVectorList());

    IntegerVector w=P.parseIntegerVector();

    printSetSetInt(stdout,regularSubdivision(m,w));

    return 0;
  }
};

static RegularSubdivisionApplication theApplication;
