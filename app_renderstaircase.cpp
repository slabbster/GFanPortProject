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
#include "renderer.h"

class RenderStaircaseApplication : public GFanApplication
{
  SimpleOption optionListOfPolynomialSets;
  IntegerOption optionMaxEntry;
  IntegerOption optionWidth;
public:
  const char *helpText()
  {
    return "This program renders a staircase diagram of a monomial initial ideal to an xfig file. The input is a Groebner basis of a (not necessarily monomial) polynomial ideal. The initial ideal is given by the leading terms in the Groebner basis. Using the -m option it is possible to render more than one staircase diagram. The program only works for ideals in a polynomial ring with three variables.\n";
  }
  RenderStaircaseApplication():
    optionListOfPolynomialSets("-m","Read multiple ideals from the input. The ideals are given as a list of lists of polynomials. For each polynomial list in the list a staircase diagram is drawn.\n"),
    optionMaxEntry("-d","Specifies the number of boxes being shown along each axis. Be sure that this number is large enough to give a correct picture of the standard monomials. The default value is 8.\n",8),
    optionWidth("-w","Width. Specifies the number of staircase diagrams per row in the xfig file. The default value is 5.\n",5)
  {
    registerOptions();
  }

  const char *name()
  {
    return "_renderstaircase";
  }

  int main()
  {
    FileParser P(Stdin);
    StandardMonomialRenderer r(Stdout);
    r.setMaxEntry(optionMaxEntry.getValue());
    r.setNumberOfDrawingsPerLine(optionWidth.getValue());

    if(optionListOfPolynomialSets.getValue())
      {
	PolynomialSetList l=P.parsePolynomialSetListWithRing();
	for(PolynomialSetList::const_iterator i=l.begin();i!=l.end();i++)
	  r.render(*i);
      }
    else
      {
	PolynomialSet p=P.parsePolynomialSetWithRing();
	r.render(p);
      }

    return 0;
  }
};

static RenderStaircaseApplication theApplication;
