#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "lp.h"

class TransposeMatrixApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation() // Not included since the program has no relation to the main programs
  {
    return false;
  }
  TransposeMatrixApplication()
  {
    registerOptions();
  }
  const char *name()
  {
    return "_transposematrix";
  }
  int main()
  {
    lpSetSolver("cddgmp");
    FileParser P(Stdin);
    IntegerVectorList v=P.parseIntegerVectorList();
    fprintf(Stderr,"Rank:%i\n",rankOfMatrix(v));
    AsciiPrinter(Stdout).printVectorList(transposeIntegerVectorList(v));
    return 0;
  }
  const char *helpText()
  {
    return "Takes a matrix and transposes it.\n";
  }
};

static TransposeMatrixApplication theApplication;
