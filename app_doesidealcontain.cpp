#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "division.h"
#include "log.h"

class DoesIdealContainApplication : public GFanApplication
{
	SimpleOption remainderOption;
	SimpleOption multiplierOption;
public:
  const char *helpText()
  {
    return "This program takes a marked Groebner basis of an ideal I and a set of polynomials on its input and tests if the polynomial set is contained in I by applying the division algorithm for each element. The output is 1 for true and 0 for false.\n";
  }
  DoesIdealContainApplication():
	  remainderOption("--remainder","Tell the program to output the remainders of the divisions rather than outputting 0 or 1."),
	  multiplierOption("--multiplier","Reads in a polynomial that will be multiplied to the polynomial to be divided before doing the division.")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_doesidealcontain";
  }
  int main()
  {
    FileParser P(Stdin);
    PolynomialSet a=P.parsePolynomialSetWithRing();
    PolynomialRing R=a.getRing();
    PolynomialSet b=P.parsePolynomialSet(R);
    Polynomial multiplier=R.one();

    if(multiplierOption.getValue())
    {
    	multiplier=P.parsePolynomial(R);
    }

    if(remainderOption.getValue())
    {
    	PolynomialSet s(a.getRing());
    	for(PolynomialSet::const_iterator i=b.begin();i!=b.end();i++)
    	{
    	  WeightReverseLexicographicTermOrder T(termorderWeight(a));
    	  s.push_back(division(multiplier* *i,a,T/*LexicographicTermOrder()*/));
    	}
    	pout<<s.getRing()<<s;
    }
    else
    {
    	bool c=true;
    	for(PolynomialSet::const_iterator i=b.begin();i!=b.end();i++)
    	{
    		Polynomial remainder=division(multiplier* *i,a,LexicographicTermOrder());
	log2 AsciiPrinter(Stderr).printString("Remainder: ");
	log2 AsciiPrinter(Stderr).printPolynomial(remainder);
	log2 AsciiPrinter(Stderr).printNewLine();
	if(!remainder.isZero())
	  {
	    log1 AsciiPrinter(Stderr).printString("Polynomial not in ideal: ");
	    log1 AsciiPrinter(Stderr).printPolynomial(multiplier* *i);
	    log1 AsciiPrinter(Stderr).printNewLine();
	    c=false;
	    break;
	  }
      }
    AsciiPrinter(Stdout).printInteger(c);
    AsciiPrinter(Stdout).printNewLine();
    }
    return 0;
  }
};

static DoesIdealContainApplication theApplication;
