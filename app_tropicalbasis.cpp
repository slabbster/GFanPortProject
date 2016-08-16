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
#include "tropicalbasis.h"
#include "log.h"

class TropicalBasisApplication : public GFanApplication
{

  SimpleOption optionHomogenize;
public:
  const char *helpText()
  {
    return "This program computes a tropical basis for an ideal defining a tropical curve. Defining a tropical curve means that the Krull dimension of R/I is at most 1 + the dimension of the homogeneity space of I where R is the polynomial ring. The input is a generating set for the ideal. If the input is not homogeneous option -h must be used.\n";
  }
  TropicalBasisApplication():
    optionHomogenize("-h","Homogenise the input before computing a tropical basis and dehomogenise the output. This is needed if the input generators are not already homogeneous.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_tropicalbasis";
  }
  int main()
  {
    FileParser P(Stdin);

    /*    StringParser P1("{1+76*a2*b*c+6*a2*c4+72*a*b3*c3,"
		    "46*b3*c+68*b3*a2*c+42*a3*c4+47*b6*a2}");
    StringParser P2("{1-79ab^4-27b^3c^2+32a^4b^2c,"
		    "-36c-35b^3c^3+11a^4b^2c^2+90a^3b^6}");
    */
    PolynomialSet theInput=P.parsePolynomialSetWithRing();
    PolynomialRing originalRing=theInput.getRing();

    if(optionHomogenize.getValue())
      {
	log1 fprintf(Stderr,"Homogenizing...\n");
	IntegerVector grading=IntegerVector::allOnes(theInput.numberOfVariablesInRing());

	/*	{  // It is not necessary to compute the true homogenization
	  WeightReverseLexicographicTermOrder t(grading);
	  buchberger(&theInput,t);
	  }*/

	PolynomialSet h=theInput.homogenization(theInput.getRing().withVariablesAppended("H"));
wallInequalities(h);

	log1 fprintf(Stderr,"The homogenized ideal:\n");
	log1 AsciiPrinter(Stderr).printPolynomialSet(h);

	theInput=h;
      }


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

    if(optionHomogenize.getValue())
      {
	theOutput=theOutput.embeddedInto(originalRing);
      }
    AsciiPrinter(Stdout).printPolynomialRing(theOutput.getRing());
    AsciiPrinter(Stdout).printNewLine();
    AsciiPrinter(Stdout).printPolynomialSet(theOutput);

    return 0;
  }
};

static TropicalBasisApplication theApplication;
