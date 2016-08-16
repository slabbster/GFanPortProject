#include "parser.h"
#include "printer.h"
#include "matrix.h"
#include "symmetry.h"
#include "gfanapplication.h"


class ComposePermutationApplication : public GFanApplication
{
  FieldOption theFieldOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program takes a list of permutations and composes them in the given order.\n";
  }
  ComposePermutationApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_composepermutations";
  }

  int main()
  {
    IntegerMatrix m=rowsToIntegerMatrix(FileParser(Stdin).parseIntegerVectorList());

    int n=m.getWidth();
    IntegerVector p=SymmetryGroup::identity(n);
    for(int i=0;i<m.getHeight();i++)
      p=SymmetryGroup::compose(p,m[i]);

    AsciiPrinter(Stderr).printVector(p);

    return 0;
  }
};

static ComposePermutationApplication theApplication;
