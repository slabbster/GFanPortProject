#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "minkowskisum.h"
#include "newtonpolytope.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"

class IsGroebnerBasisApplication : public GFanApplication
{
  FieldOption theFieldOption;
  SimpleOption optionUGBTest;
public:
  bool includeInDefaultInstallation() // Not included since it requires Christophe's Minkowski sum program
  {
    return false;
  }
  const char *helpText()
  {
    return "Checks if an unmarked polynomial set is a Groebner basis with respect to some termorder.\n";
  }
  IsGroebnerBasisApplication():
    optionUGBTest("--ugb","Universal Groebner Basis test. Check if something is a Groebner basis with respect to any termorder.")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_isgroebnerbasis";
  }
  int main()
  {
    lpSetSolver("cddgmp");
    FileParser P(Stdin);
    PolynomialSet g=P.parsePolynomialSetWithRing();
    PolynomialRing theRing=g.getRing();

    IntegerVectorListList polytopes;

    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
      {
	polytopes.push_back(newtonPolytope(*i));
      }

    IntegerVectorListList sums;
    IntegerVectorList polytope=minkowski(polytopes,&sums);

    fprintf(Stderr,"Number of extreme vertices in Minkowski sum: %i\n",polytope.size());

    bool isGroebnerBasis=false;
    int counter=0;
    for(IntegerVectorListList::const_iterator i=sums.begin();i!=sums.end();i++)
      {
	fprintf(Stderr,"Checking vertex %i\n",counter++);
	IntegerVectorList::const_iterator j=i->begin();

	for(PolynomialSet::iterator k=g.begin();k!=g.end();k++)
	  {
	    k->mark(Monomial(theRing,*j));
	    j++;
	  }
	IntegerVectorList normals=wallInequalities(g);
	if(hasInteriorPoint(normals,true))
	  {
	    if(isMarkedGroebnerBasis(g))
	      {
		AsciiPrinter(Stderr).printPolynomialSet(g);
		isGroebnerBasis=true;
		if(!optionUGBTest.getValue())break;
	      }
	    else
	      {
		fprintf(Stderr,"Marking is not a Groebner basis\n");
		if(optionUGBTest.getValue())
		  {
		    AsciiPrinter(Stderr).printPolynomialSet(g);
		    assert(0);
		  }
	      }
	  }
	else
	  fprintf(Stderr,"Cone has no positive interior point\n");
      }

    fprintf(Stdout,isGroebnerBasis?"true\n":"false\n");
    return 0;
  }
};

static IsGroebnerBasisApplication theApplication;
