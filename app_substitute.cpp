#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"

typedef list<string> StringList;

class SubstituteApplication : public GFanApplication
{
public:
  const char *helpText()
  {
    return "This program changes the variable names of a polynomial ring. The input is a polynomial ring, a polynomial set in the ring and a new polynomial ring with the same coefficient field but different variable names. The output is the polynomial set written with the variable names of the second polynomial ring.\n"
      "Example:\n"
      "Input:\n"
      "Q[a,b,c,d]{2a-3b,c+d}Q[b,a,c,x]\n"
      "Output:\n"
      "Q[b,a,c,x]{2*b-3*a,c+x}\n";
  }
  SubstituteApplication()
  {
    registerOptions();
  }
  const char *name()
  {
    return "_substitute";
  }
  int main()
  {
    FileParser P(Stdin);
    PolynomialRing r=P.parsePolynomialRing();
    PolynomialSet s=P.parsePolynomialSet(r);
    PolynomialRing r2=P.parsePolynomialRing();
    AsciiPrinter(Stdout).printPolynomialRing(r2);
    AsciiPrinter(Stdout).printNewLine();
    AsciiPrinter(Stdout).printPolynomialSet(s.embeddedInto(r2));

    return 0;
  }
};

static SubstituteApplication theApplication;
