#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "binomial.h"
#include "latticeideal.h"
#include "field_rationals.h"
#include "lll.h"

class LatticeIdealApplication : public GFanApplication
{
  SimpleOption toricOption;
  SimpleOption justConvertOption;
public:
  bool includeInDefaultInstallation()
  {
    return true;
  }
  LatticeIdealApplication():
    toricOption("-t","Compute the toric ideal of the matrix whose rows are given on the input instead."),
	justConvertOption("--convert","Does not do any computation, but just converts the vectors to binomials.")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_latticeideal";
  }
  int main()
  {
    FileParser P(Stdin);

    IntegerVectorList ivl=P.parseIntegerVectorList();
    IntegerMatrix A=rowsToIntegerMatrix(ivl).transposed();

    if(toricOption.getValue())
      A=latticeKernelOfTransposed(A).transposed();

    IntegerVectorList b;
    if(justConvertOption.getValue())
    	b=A.transposed().getRows();
    else
    	b=latticeIdealRevLex(A);

    PolynomialRing theRing(Q,A.getHeight());
    PolynomialSet g(theRing);
    for(IntegerVectorList::const_iterator i=b.begin();i!=b.end();i++)
      g.push_back(integerVectorToBinomial(*i,theRing));

    AsciiPrinter(Stdout).printPolynomialRing(g.getRing());
    AsciiPrinter(Stdout).printPolynomialSet(g);
    return 0;
  }
  const char *helpText()
  {
    return "This program computes the lattice ideal of a lattice. The input is a list of generators for the lattice.\n";
  }
};

static LatticeIdealApplication theApplication;
