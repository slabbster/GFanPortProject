/*
  TO DO:
  Remove Timer
 */

#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "buchberger.h"
#include "reversesearch.h"
#include "termorder.h"
#include "genericwalk.h"
#include "gfanapplication.h"
#include "field_rationalfunctions2.h"
#include "timer.h"
#include "log.h"
//#include "singular.h"

static Timer globalTimer("Global-timer",1);

class BuchbergerApplication : public GFanApplication
{
  SimpleOption optionReadWeightVector;
  SimpleOption optionReverse;
  SimpleOption optionWalk;
  SimpleOption optionGenericWalk;
  IntegerOption optionParameters;
  //  SimpleOption optionPerformanceTimer;
public:
  const char *helpText()
  {
    return "This program computes a reduced lexicographic Groebner basis of the polynomial ideal given as input. The default behavior is to use Buchberger\'s algorithm. The ordering of the variables is $a>b>c...$ (assuming that the ring is Q[a,b,c,...]).\n";
  }
  BuchbergerApplication():
    optionReadWeightVector("-w","Compute a Groebner basis with respect to a degree lexicographic order with $a>b>c...$ instead. The degrees are given by a weight vector which is read from the input after the generating set has been read.\n"),
    optionReverse("-r","Use the reverse lexicographic order (or the reverse lexicographic order as a tie breaker if -w is used). The input must be homogeneous if the pure reverse lexicographic order is chosen. Ignored if -W is used.\n"),
    optionWalk("-W","Do a Groebner walk. The input must be a minimal Groebner basis. If -W is used -w is ignored.\n"),
    optionGenericWalk("-g","Do a generic Groebner walk. The input must be homogeneous and must be a minimal Groebner basis with respect to the reverse lexicographic term order. The target term order is always lexicographic. The -W option must be used.\n"),
    optionParameters("--parameters","With this option you can specify how many variables to treat as parameters instead of variables. This makes it possible to do computations where the coefficient field is the field of rational functions in the parameters. This does not work well at the moment.",0)
    //,optionPerformanceTimer("-T","enable performance timer")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_buchberger";
  }

  int main()
  {
    TermOrder *myOrder;

    PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();
    if(optionParameters.getValue())g=makeVariablesParameters(makeVariablesParameters(g.getRing(),optionParameters.getValue()),g);

    if(optionReadWeightVector.getValue())
      {
	IntegerVector w=FileParser(Stdin).parseIntegerVector();
	if(optionReverse.getValue())
	  myOrder = new WeightReverseLexicographicTermOrder(w);
	else
	  myOrder = new WeightTermOrder(w);
      }
    else
      if(optionReverse.getValue())
	myOrder =  new ReverseLexicographicTermOrder;
      else
	myOrder = new LexicographicTermOrder;

    globalTimer.on();

    if(optionWalk.getValue())
      {
	g.scaleMarkedCoefficientsToOne();
	log1 fprintf(Stderr,"Auto-reducing input...\n");
	autoReduce(&g,LexicographicTermOrder());
	log1 fprintf(Stderr,"Walking...\n");
	if(optionGenericWalk.getValue())
	  {
	    g=genericWalk(g, ReverseLexicographicTermOrder(), LexicographicTermOrder());
	  }
	else
	  {
	    g=ReverseSearch(LexicographicTermOrder()).findRoot(g);
	  }
      }
    else
      {
	//fprintf(Stderr,"TESTING Singular interface\n");
	//singularBuchberger(&g,*myOrder);
	buchberger(&g,*myOrder);


	//    AsciiPrinter(Stderr).printPolynomialSet(g);
	// fprintf(Stderr,"\n\n\n");

	g.scaleMarkedCoefficientsToOne();//	buchberger(&g,*myOrder);
	autoReduce(&g,LexicographicTermOrder());
      }

    globalTimer.off();

    AsciiPrinter(Stdout).printPolynomialRing(g.getRing());
    printf("\n");
    AsciiPrinter(Stdout).printPolynomialSet(g);

    //    if(optionPerformanceTimer.getValue())Timer::printList();

    return 0;
  }
};

static BuchbergerApplication theApplication;

