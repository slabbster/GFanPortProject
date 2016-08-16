#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "termorder.h"
#include "ep_standard.h"
#include "ep_xfig.h"
#include "genericwalk.h"
#include "gfanapplication.h"
#include "timer.h"

class WalkApplication : public GFanApplication
{
  FieldOption theFieldOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  WalkApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_walk";
  }

  int main()
  {
    lpSetSolver("cddgmp");


    FileParser p(Stdin);
    PolynomialSet g=p.parsePolynomialSetWithRing();

    //  StandardGradedLexicographicTermOrder t1;
    WeightReverseLexicographicTermOrder t1(StringParser("(1,1,1,1,1,1)").parseIntegerVector());
    AsciiPrinter P(Stderr);

    assert(g.checkMarkings(t1));
      //    g.mark(t1);
      // AsciiPrinter(Stdout).printPolynomialSet(g);

    t1.printMatrix(P,6);
    LexicographicTermOrder t2;
    //    LexicographicTermOrder t1(1);
    //    LexicographicTermOrder t2;

    fprintf(Stderr,"Computing starting Groebner basis\n");

    autoReduce(&g,t1);
    //    buchberger(&g,t1);

    AsciiPrinter(Stdout).printPolynomialSet(g);

    fprintf(Stderr,"Walking with perturbation\n");


    //    int n=g.numberOfVariablesInRing();
    //    for(int i=n;i>=0;i--)
    //    for(int i=n-4;i>=0;i--)
    {
      PolynomialSet g2(g.getRing());

      //      g2=genericWalkPerturbation(g,t1,t2,i,i);
      g2=genericWalkPerturbation(g,t1,t2,6,6);

      AsciiPrinter(Stdout).printPolynomialSet(g2);
    }

    return 0;
  }
  const char *helpText()
  {
    return "This program walks. THIS PROGRAM NEEDS TO BE FIXED.\n";
  }
};

static WalkApplication theApplication;

