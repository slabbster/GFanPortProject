#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "polyhedralfan.h"
#include "breadthfirstsearch.h"
#include "termorder.h"
#include "ep_standard.h"
#include "ep_xfig.h"
#include "gfanapplication.h"
#include "timer.h"

class StatsApplication : public GFanApplication
{
public:
  const char *helpText()
  {
    return "This program takes a list of reduced Groebner bases for the same ideal and computes various statistics. The following information is listed: the number of bases in the input, the number of variables, the dimension of the homogeneity space, the maximal total degree of any polynomial in the input and the minimal total degree of any basis in the input, the maximal number of polynomials and terms in a basis in the input.\n";
  }
  StatsApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_stats";
  }

  int main()
  {
    FileParser p(Stdin);

    PolynomialRing theRing=p.parsePolynomialRing();

    int c=p.nextNonBlank();

    int dmin=-1;
    int dmax=-1;
    int homog=-1;
    int n=-1;
    int maxNumberOfTerms=-1;
    int maxNumberOfPolynomials=-1;

    int counter=0;
    assert(p.isLeftBracket(c));
    do
      {
	PolynomialSet g=p.parsePolynomialSet(theRing);
	if(homog==-1)
	  {
	    n=g.numberOfVariablesInRing();
	    homog=dimensionOfHomogeneitySpace(g);
	  }
	int d=g.totalDegree();

	if(dmin==-1 || d<dmin)dmin=d;
	if(dmax==-1 || d>dmax)dmax=d;
	if(maxNumberOfTerms<(int)g.totalNumberOfTerms())maxNumberOfTerms=(int)g.totalNumberOfTerms();
	if(maxNumberOfPolynomials<(int)g.size())maxNumberOfPolynomials=(int)g.size();
	counter++;
      }
    while((c=p.nextNonBlank())==',');
    assert(p.isRightBracket(c));

    fprintf(Stdout,"Number of reduced Groebner bases: %i\n",counter);
    fprintf(Stdout,"Number of variables: %i\n",n);
    fprintf(Stdout,"Dimension of homogeneity space: %i\n",homog);
    fprintf(Stdout,"Maximal total degree of a Groebner basis: %i\n",dmax);
    fprintf(Stdout,"Minimal total degree of a Groebner basis: %i\n",dmin);
    pout<<"Maximal number of polynomials in Groebner basis: "<<maxNumberOfPolynomials<<"\n";
    pout<<"Maximal number of terms in Groebner basis: "<<maxNumberOfTerms<<"\n";
    return 0;
  }
};

static StatsApplication theApplication;

