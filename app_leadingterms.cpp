#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "termorder.h"
#include "gfanapplication.h"
#include "wallideal.h"

class LeadingTermsApplication : public GFanApplication
{
  SimpleOption optionMultipleSets;
public:
  const char *helpText()
  {
    return "This program converts a list of polynomials to a list of their leading terms.\n";
  }
  LeadingTermsApplication():
    optionMultipleSets("-m","Do the same thing for a list of polynomial sets. That is, output the set of sets of leading terms.\n")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_leadingterms";
  }

  void p(PolynomialSet const &g)
  {
    PolynomialSet LT(g.getRing());

    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
      {
	LT.push_back(Polynomial(i->getMarked()));
      }

    AsciiPrinter(Stdout).printPolynomialSet(LT);
  }

  int main()
  {
    FileParser P(Stdin);

    if(!optionMultipleSets.getValue())
      {
	PolynomialSet s=P.parsePolynomialSetWithRing();
	AsciiPrinter(Stdout).printPolynomialRing(s.getRing());
	AsciiPrinter(Stdout).printNewLine();
	p(s);
      }
    else
      {
	PolynomialSetList l=P.parsePolynomialSetListWithRing();
	assert(l.size()!=0);
	AsciiPrinter(Stdout).printPolynomialRing(l.begin()->getRing());
	AsciiPrinter(Stdout).printNewLine();
	fprintf(Stdout,"{\n");
	for(PolynomialSetList::const_iterator i=l.begin();i!=l.end();i++)
	  {
	    if(i!=l.begin())fprintf(Stdout,",\n");
	    p(*i);
	  }
	fprintf(Stdout,"}\n");
      }

    return 0;
  }
};

static LeadingTermsApplication theApplication;
