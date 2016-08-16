#include <iostream>
#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "buchberger.h"
#include "reversesearch.h"
#include "termorder.h"
#include "genericwalk.h"
#include "gfanapplication.h"
#include "timer.h"
#include "tropical2.h"
#include "log.h"
#include "field_rationals.h"

class SPolynomialApplication : public GFanApplication
{
  SimpleOption optionTakeInitialForms;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program takes a list of polynomials and a weight vector. The output is the set of all...\n";
  }
  SPolynomialApplication():
    optionTakeInitialForms("-f",".\n")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_spolynomial";
  }

  int main()
  {
    PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();
    IntegerVector w=FileParser(Stdin).parseIntegerVector();

    //    WeightReverseLexicographicTermOrder T(w);
    WeightTermOrder T(w);

    g.markAndScale(T);
    PolynomialSet out(g.getRing());

    int I=0;
    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++,I++)
      {
	int J=0;
	for(PolynomialSet::const_iterator j=g.begin();j!=i;j++,J++)
	{
	  Polynomial s=sPolynomial(*i,*j);
	  IntegerVector cancelingTerm=max(i->getMarked().m.exponent,j->getMarked().m.exponent);
	  int cancelDegree=dotLong(cancelingTerm,w);
	  if(optionTakeInitialForms.getValue())
	    {
	      if(cancelDegree!=s.degree(w))
		{
		  out.push_back(initialForm(s,w));
		  if(initialForm(s,w).numberOfTerms()<=8)cerr<<"S("<<I<<","<<J<<")"<<endl;
		}
	    }
	  else
	    out.push_back(s);
	}
      }
    //    AsciiPrinter(Stdout).printPolynomialSet(out);

    PolynomialSet out2(out.getRing());
    for(PolynomialSet::const_iterator i=out.begin();i!=out.end();i++)
      {
	cerr<< i->numberOfTerms() << "    ";
	if(i->numberOfTerms()<=8)
	  {
	    AsciiPrinter(Stdout)<<"\n"<< *i<<"\n";
	    out2.push_back(*i);
	  }
      }
    out2.saturate();
    out2=out2.embeddedInto(PolynomialRing(Q,15));

    AsciiPrinter(Stdout).printPolynomialRing(out2.getRing());
    AsciiPrinter(Stdout).printPolynomialSet(out2);

    return 0;
  }
};

static SPolynomialApplication theApplication;

