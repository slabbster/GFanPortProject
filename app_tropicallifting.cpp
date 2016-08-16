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
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "tropical.h"
#include "tropical2.h"
#include "halfopencone.h"
#include "multiplicity.h"
#include "tropicalbasis.h"
#include "log.h"

class TropicalLiftingApplication : public GFanApplication
{
  SimpleOption optionNoMultiplicity;
  //  FieldOption theFieldOption;
  IntegerOption optionNumberOfNegativeVariables;
  SimpleOption optionOnlyOutputChoices;
public:
  const char *helpText()
  {
    return "This program is part of the Puiseux lifting algorithm implemented in Gfan and Singular. The Singular part of the implementation can be found in:\n\n"
      "Anders Nedergaard Jensen, Hannah Markwig, Thomas Markwig:\n tropical.lib. A SINGULAR 3.0 library for computations in tropical geometry, 2007 \n\n"
      "See also\n\nhttp://www.mathematik.uni-kl.de/~keilen/de/tropical.html\n\n"
 "and the paper\n\n Jensen, Markwig, Markwig: \"An algorithm for lifting points in a tropical variety\".\n\n"
      "Example:\n\n"
      "Run Singular from the directory where tropical.lib is located.\n"
      "Give the following sequence of commands to Singular:\n\n"
      "LIB \"tropical.lib\";\n"
      "ring R=0,(t,x,y,z),dp;\n"
      "ideal i=-y2t4+x2,yt3+xz+y;\n"
      "intvec w=1,-2,0,2;\n"
      "list L=tropicallifting(i,w,3);\n"
      "displaytropicallifting(L,\"subst\");\n"
      "This produces a Puiseux series solution to i with valuation (2,0,-2)\n";

  }
  TropicalLiftingApplication():
    optionNoMultiplicity("--noMult","Disable the multiplicity computation."),
    optionNumberOfNegativeVariables("-n","Number of variables that should have negative weight."),
    optionOnlyOutputChoices("-c","Only output a list of vectors being the possible choices.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_tropicallifting";
  }
  int main()
  {
    FileParser P(Stdin);

    PolynomialSet theInput=P.parsePolynomialSetWithRing();
    PolynomialRing theRing=theInput.getRing();

    int n=theInput.numberOfVariablesInRing();

    int homog=-1;
    {
      log1 debug.printString("Computing homogeneity space\n");
      PolynomialSet idealGroebnerBasis=theInput;
      IntegerVector grading=IntegerVector::allOnes(theInput.numberOfVariablesInRing());
      WeightReverseLexicographicTermOrder t(grading);
      buchberger(&idealGroebnerBasis,t);
      PolyhedralCone hspace=homogeneitySpace(idealGroebnerBasis);
      IntegerVectorList hv=hspace.dualCone().getEquations();
      homog=hv.size();
      log1 debug.printString("..done homogeneity space.\n");
    }

    PolynomialSet theOutput=tropicalBasisOfCurve(n,theInput,0,homog);

    /*    if(optionHomogenize.getValue())
      {
	theOutput=theOutput.deHomogenization();
      }
    */

    if(!optionOnlyOutputChoices.getValue())
      AsciiPrinter(Stdout).printPolynomialSet(theOutput);



    //    PolyhedralFan F=tropicalHyperSurfaceIntersectionClosed(n, theOutput);
    PolyhedralFan F=tropicalPrincipalIntersection(n, theOutput,homog); // dimension of lineality space could be computed to speed up computations


    //    fprintf(Stderr,"---------------------------\n");
    //AsciiPrinter(Stdout).printPolyhedralFan(F);
    //   fprintf(Stderr,"---------------------------\n");

    if(!optionOnlyOutputChoices.getValue())
      fprintf(Stdout,"\nA list of relative interior points:\n");

    IntegerVectorList rays=F.getRelativeInteriorPoints();
    if(!optionOnlyOutputChoices.getValue())
      AsciiPrinter(Stdout).printVectorList(rays);



    if(!optionOnlyOutputChoices.getValue())
      if(!optionNoMultiplicity.getValue())
	{
	  fprintf(Stdout,"Multiplicities:\n");
	  for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
	    {
	      if(!(i->isZero()))
		{

		  AsciiPrinter(Stdout).printVector(*i);
		  fprintf(Stdout,"   :");
		  AsciiPrinter(Stdout).printInteger(multiplicity(initialIdeal(theInput,*i)));
		  fprintf(Stdout,"\n");
		}
	    }
	}

    int t=optionNumberOfNegativeVariables.getValue();
    assert(t<=n);

    IntegerVectorList equations,inequalities;
    for(int i=0;i<t;i++)
      inequalities.push_back(-IntegerVector::standardVector(n,i));
    for(int i=t;i<n;i++)
      equations.push_back(IntegerVector::standardVector(n,i));
    PolyhedralCone C(inequalities,equations,n);



    if(!optionOnlyOutputChoices.getValue())
      fprintf(Stdout,"Possible choices:\n");
    for(PolyhedralConeList::const_iterator i=F.conesBegin();i!=F.conesEnd();i++)

    //    for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
      {
	PolyhedralCone C2=intersection(C,*i);
	IntegerVector v=C2.getRelativeInteriorPoint();

	//	    AsciiPrinter(Stdout).printVector(v);

	bool isNeg=true;
	for(int j=0;j<t;j++)
	  {
	    if(v[j]>=0)isNeg=false;
	  }
	if(isNeg)
	  {
	    AsciiPrinter(Stdout).printVector(v);fprintf(Stdout,"\n");
	  }
      }

    return 0;
  }
};

static TropicalLiftingApplication theApplication;
