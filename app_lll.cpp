#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "lll.h"

class LLLApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  LLLApplication()
  {
    registerOptions();
  }
  const char *name()
  {
    return "_lll";
  }
  int main()
  {
    FileParser P(Stdin);

    IntegerVectorList ivl=P.parseIntegerVectorList();
    IntegerMatrix A=rowsToIntegerMatrix(ivl);


    AsciiPrinter Q(Stdout);

    Q.printString("Input generators:\n");
    AsciiPrinter(Stdout).printVectorList(ivl);

    IntegerMatrix Minv;
    IntegerMatrix A2=A;
    IntegerMatrix M=mlll(A,&Minv);
    Q.printString("LLL basis:\n");
    AsciiPrinter(Stdout).printVectorList(A.getRows());
    Q.printString("Transformation:\n");
    AsciiPrinter(Stdout).printVectorList(M.getRows());
    Q.printString("Inverse transformation:\n");
    AsciiPrinter(Stdout).printVectorList(Minv.getRows());

    Q.printString("Basis for integer kernal:\n");
    AsciiPrinter(Stdout).printVectorList(latticeKernelOfTransposed(A2).getRows());

    return 0;
  }
  const char *helpText()
  {
    return "Tests the LLL algorithm. Input is a list of lattice generators.\n";
  }
};

static LLLApplication theApplication;
