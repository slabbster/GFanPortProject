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
#include "dimension.h"
#include "tropical.h"
#include "tropical2.h"
#include "log.h"

class TropicalBruteForceApplication : public GFanApplication
{
  //SimpleOption OptionTPlane;
public:
  const char *helpText()
  {
    return "This program takes a marked reduced Groebner basis for a homogeneous ideal and computes the tropical "
      "variety of the ideal as a subfan of the Groebner fan. "
      "The program is slow but works for any homogeneous ideal. If "
      "you know that your ideal is prime over the complex numbers or you "
      "simply know that its tropical variety is pure and connected in "
      "codimension one then use gfan_tropicalstartingcone and "
      "gfan_tropicaltraverse instead.\n";
  }
  TropicalBruteForceApplication()
  {
    registerOptions();
  }

  const char *name()//:
  //    OptionTPlane("--tplane","This option intersect the resulting
  {
    return "_tropicalbruteforce";
  }

  int main()
  {
    FileParser p(Stdin);

    PolynomialRing r=p.parsePolynomialRing();
    PolynomialSet G=p.parsePolynomialSet(r);

    if(!isMarkingConsistent(G))
      {
	fprintf(Stderr,"Input polynomial set is not marked consistently.\n");
	assert(0);
      }
    if(!isMarkedGroebnerBasis(G))
      {
	fprintf(Stderr,"Input polynomial set is not a marked Groebner basis.\n");
	assert(0);
      }

    autoReduce(&G, LexicographicTermOrder());

    int homog=dimensionOfHomogeneitySpace(G);
    int n=r.getNumberOfVariables();
    IntegerVector f=IntegerVector(n-homog+1);
    IntegerVector fTrop=IntegerVector(n-homog+1);

    //    ReverseLexicographicTermOrder A;
    LexicographicTermOrder A;
    ReverseSearch rs(A);
    EnumerationTargetCollector C;
    rs.setEnumerationTarget(&C);
    rs.enumerate(G);
    PolynomialSetList glist=C.getList();

    PolyhedralFan result(n);

    for(PolynomialSetList::const_iterator j=glist.begin();j!=glist.end();j++)
      {
	PolynomialSet g=*j;
	PolyhedralCone c=groebnerCone(g,true);
	PolyhedralFan F(n);
	F.insert(c);
	for(int i=0;i<n-homog+1;i++)
	  {
	    F.removeAllLowerDimensional();
	    PolyhedralFan F2(n);
	    while(!F.isEmpty())
	      {
		PolyhedralCone K=F.highestDimensionalCone();
		F.remove(K);

		IntegerVector w=K.getRelativeInteriorPoint();
		WeightTermOrder myOrder(w);
		if(g.checkMarkings(myOrder))
		  {
		    f[n-homog-i]++;
		    F2.insert(K);
		    if(!containsMonomial(initialFormsAssumeMarked(g,w)))
		      {
			result.insert(K);
			fTrop[n-homog-i]++;
		      }
		  }
	      }
	    F=F2.facetComplex();
	  }
      }

    AsciiPrinter P(Stdout);
    result.printWithIndices(&P, FPF_default);
    //    result.printWithIndices(&P, false, 0, false);

    AsciiPrinter Q(Stderr);
    log1 Q.printString("F-vector of Groebner fan:\n");
    log1 Q.printVector(f);
    log1 Q.printString("\nF-vector of tropical variety:\n");
    log1 Q.printVector(fTrop);
    log1 Q.printNewLine();

    return 0;
  }
};

static TropicalBruteForceApplication theApplication;

