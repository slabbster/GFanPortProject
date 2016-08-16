#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "buchberger.h"
#include "wallideal.h"
#include "termorder.h"
#include "gfanapplication.h"
#include "tropical2.h"
#include "matrix.h"
#include <iostream>

class EvaluationApplication : public GFanApplication
{
	SimpleOption optionComplex;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program evaluates a list of polynomials in a point.\n";
  }
  EvaluationApplication():
	  optionComplex("--complex","Read and evaluate using complex numbers rather than real numbers.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_evaluation";
  }

  int main()
  {
    FileParser P(Stdin);

    PolynomialSet g=P.parsePolynomialSetWithRing();


    if(0){
    	for(int i1=1;i1<20;i1++)
        	for(int i2=1;i2<20;i2++)
            	for(int i3=1;i3<20;i3++)
            	{
					FloatVector v(4);
					v[0]=i1/20.0;
					v[1]=i2/20.0;
					v[2]=i3/20.0;
					v[3]=1-v[1]-v[2]-v[0];
                	pout.printFloatVector(g.evaluateFloat(v));
            	}
    }

    if(optionComplex.getValue())
    {
    	ComplexVector x=P.parseComplexVector();
    	pout.printComplexVector(g.evaluateComplex(x));
    }
    else
    {
    	FloatVector x=P.parseFloatVector();
    	pout.printFloatVector(g.evaluateFloat(x));
    }

    pout<<"\n";
    return 0;
  }
};

static EvaluationApplication theApplication;
