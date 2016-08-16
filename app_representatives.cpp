#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "minkowskisum.h"
#include "newtonpolytope.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "tropical.h"
#include "division.h"
#include "bergman.h"
#include "tropical2.h"
#include "dimension.h"
#include "timer.h"

class RepresentativesApplication : public GFanApplication
{
  FieldOption theFieldOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program computes takes generaters for a subgroup of S_n and a list of n-dimensional integer vectors. The output is a list of vectors, one from each orbit of elements of the list. \n";
  }
  RepresentativesApplication()
  {
    registerOptions();
  }
  const char *name()
  {
    return "_representatives";
  }
  int main()
  {
    PolynomialSetList tropical;

    lpSetSolver("cddgmp");
    FileParser P(Stdin);

    AsciiPrinter p(Stdout);


    IntegerVectorList generators=P.parseIntegerVectorList();
    assert(generators.size()!=0);
    int n=generators.begin()->size();
    SymmetryGroup s(n);
    s.computeClosure(generators);
    s.print(Stderr);
    fprintf(Stderr,"\n");

    IntegerVectorList vList=P.parseIntegerVectorList();

    IntegerVectorList rep;

    for(IntegerVectorList::const_iterator i=vList.begin();i!=vList.end();i++)
      {
	bool found=false;
	for(IntegerVectorList::const_iterator j=rep.begin();j!=rep.end();j++)
	  {
	    if(i->sum()==j->sum())
	      {
		for(SymmetryGroup::ElementContainer::const_iterator k=s.elements.begin();k!=s.elements.end();k++)
		  if(SymmetryGroup::compose(*k,*j)==*i)
		    {
		      found=true;
		      break;
		    }
	      }
	    if(found)break;
	  }
	if(!found)rep.push_back(*i);
      }

    p.printVectorList(rep);

    return 0;
  }
};

static RepresentativesApplication theApplication;
