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
#include "gfanapplication.h"
#include "symmetry.h"

class PolynomialSetUnionApplication : public GFanApplication
{
  SimpleOption optionSpecialMode;
  class PolynomialCmp
  {
  public:
    bool operator()(const Polynomial &a, const Polynomial &b)const
    {
      if(a.getMarked().m.exponent.sum()<b.getMarked().m.exponent.sum())return true;
      if(a.getMarked().m.exponent.sum()>b.getMarked().m.exponent.sum())return false;
      return LexicographicTermOrder()(a.getMarked().m.exponent,b.getMarked().m.exponent);
    }
  };
public:
  const char *helpText()
  {
    return "This program computes the union of a list of polynomial sets given as input. The polynomials must all belong to the same ring. The ring is specified on the input. After this follows the list of polynomial sets.\n";
  }
  PolynomialSetUnionApplication():
    optionSpecialMode("-s","Sort output by degree.\n")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_polynomialsetunion";
  }
  int main()
  {
    FileParser P(Stdin);
    PolynomialRing r=P.parsePolynomialRing();
    int c=P.nextNonBlank();
    bool first=true;
    assert(P.isLeftBracket(c));
    PolynomialSet s=P.parsePolynomialSet(r);
    c=P.nextNonBlank();
    while(','==c || first)
      {
	PolynomialSet temp=P.parsePolynomialSet(r);
	//	if(optionSpecialMode.getValue())temp.markAndScale(LexicographicTermOrder());
	s.unionSet(temp);
	c=P.nextNonBlank();
	first=false;
      }
    assert(P.isRightBracket(c));

    if(optionSpecialMode.getValue())s.sort(PolynomialCmp());

    AsciiPrinter(Stdout).printPolynomialRing(r);
    AsciiPrinter(Stdout).printNewLine();
    AsciiPrinter(Stdout).printPolynomialSet(s);
    return 0;
  }
};

static PolynomialSetUnionApplication theApplication;
