#include "dimension.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "division.h"
#include "field_rationals.h"

class SturmSequenceApplication : public GFanApplication
{
  SimpleOption evaluateOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "Given a polynomial this program computes the Sturm sequence of polynomials as defined in Sturm's Theorem.\n";
  }
  SturmSequenceApplication():
    evaluateOption("-e","Evaluate the sequence in a set of points given in a vector. (Only integer vectors are supported at the moment).")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_sturmsequence";
  }
  int main()
  {
    FileParser P(Stdin);
    Polynomial f1=P.parsePolynomialWithRing();

    Polynomial f2=f1.derivative();

    PolynomialRing theRing=f1.getRing();
    PolynomialSet result(theRing);
    result.push_back(f1);
    while(!f2.isZero())
      {
	result.push_back(f2);
	PolynomialSet g(theRing);
	Polynomial temp=f2;
	g.push_back(f2);
	g.markAndScale(LexicographicTermOrder());
	f2=(f1-f1)-division(f1,g,LexicographicTermOrder());
	f1=temp;
      }
    AsciiPrinter(Stdout).printPolynomialSet(result);

    if(evaluateOption.getValue())
      {
	IntegerVector v=P.parseIntegerVector();
	for(int i=0;i<v.size();i++)
	  {
	    FieldElement x=Q.zHomomorphism(v[i]);
	    AsciiPrinter(Stdout).printString("Evaluating in ");
	    AsciiPrinter(Stdout).printFieldElement(x);
	    AsciiPrinter(Stdout).printNewLine();
	    for(PolynomialSet::const_iterator j=result.begin();j!=result.end();j++)
	      {
		AsciiPrinter(Stdout).printFieldElement(j->evaluate(x));
		AsciiPrinter(Stdout).printNewLine();
	      }
	  }
      }

    return 0;
  }
};

static SturmSequenceApplication theApplication;
