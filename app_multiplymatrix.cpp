#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"

class MultiplyMatrixApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation() // Not included since the program has no relation to the main programs
  {
    return false;
  }
  MultiplyMatrixApplication()
  {
    registerOptions();
  }
  const char *name()
  {
    return "_multiplymatrix";
  }
  int main()
  {
    FileParser P(Stdin);
    IntegerVectorList a=P.parseIntegerVectorList();
    IntegerVectorList b=P.parseIntegerVectorList();
    AsciiPrinter(Stdout).printVectorList(multiplyIntegerVectorList(a,b));
    return 0;
  }
  const char *helpText()
  {
    return "Takes two matrices and multiplies them.\n";
  }
};

static MultiplyMatrixApplication theApplication;
