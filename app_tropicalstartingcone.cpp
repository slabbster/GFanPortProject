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
#include "log.h"

class TropicalStartingConeApplication : public GFanApplication
{
  SimpleOption useThisConeOption;
  SimpleOption dimensionOption;
  SimpleOption optionStableIntersection;
public:
  bool includeInDefaultInstallation()
  {
    return true;
  }
  const char *helpText()
  {
    return "This program attempts to compute a starting pair of marked reduced Groebner bases to be used as input for gfan_tropicaltraverse. If unsuccessful the program will say so. The input is a homogeneous ideal whose tropical variety is a pure d-dimensional polyhedral complex.\n";
  }
  TropicalStartingConeApplication():
    useThisConeOption("-g","Tell the program that the input is already a reduced Groebner basis."),
    dimensionOption("-d","Output dimension information to standard error."),
    optionStableIntersection("--stable","Find starting cone in the stable intersection or, equivalently, pretend that the coefficients are genereric.")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_tropicalstartingcone";
  }
  int main()
  {
    lpSetSolver("cddgmp");
    FileParser P(Stdin);

    PolynomialSet g=P.parsePolynomialSetWithRing();
    if(!useThisConeOption.getValue() && !optionStableIntersection.getValue())
      {
	buchberger(&g,StandardGradedLexicographicTermOrder());
      }
    //    assert(!containsMonomial(g));
    if(dimensionOption.getValue())
      {
	fprintf(Stderr,"Krull dimension of input ideal: %i\n",krullDimension(g));
	fprintf(Stderr,"Dimension of homogeneity space for full ideal: %i\n",dimensionOfHomogeneitySpace(g));
      }
    PolynomialSet fullBasis(g.getRing());
    log0 fprintf(Stderr,"guessing..\n");
    PolynomialSet g2=(optionStableIntersection.getValue())?
          guessInitialIdealWithoutMonomialStably(g,&fullBasis,false):
          guessInitialIdealWithoutMonomial(g,&fullBasis,false);

        //    PolynomialSet g2=guessInitialIdealWithoutMonomial(g,&fullBasis,true); //CHANGETHIS to false

    if(dimensionOption.getValue())
      {
	fprintf(Stderr,"Krull dimension of input ideal: %i\n",krullDimension(g));
	fprintf(Stderr,"Dimension of homogeneity space for full ideal: %i\n",dimensionOfHomogeneitySpace(g));
	fprintf(Stderr,"Dimension of homogeneity space for initial ideal: %i\n",dimensionOfHomogeneitySpace(g2));
      }
    AsciiPrinter(Stdout).printPolynomialRing(g2.getRing());
    AsciiPrinter(Stdout).printNewLine();
    AsciiPrinter(Stdout).printPolynomialSet(g2);
    AsciiPrinter(Stdout).printPolynomialSet(fullBasis);


    return 0;
  }
};

static TropicalStartingConeApplication theApplication;
