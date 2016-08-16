#include "parser.h"
#include "printer.h"
#include "saturation.h"
#include "tropical.h"
#include "gfanapplication.h"

class SaturationApplication : public GFanApplication
{
  SimpleOption optionIsHomogeneous;
  SimpleOption optionNoIdeal;
public:
  const char *helpText()
  {
    return "This program computes the saturation of the input ideal with the product of the variables x_1,...,x_n. The ideal does not have to be homogeneous.\n";
  }
  SaturationApplication():
    optionIsHomogeneous("-h","Tell the program that the input is a homogeneous ideal (with homogeneous generators).\n"),
    optionNoIdeal("--noideal","Do not treat input as an ideal but just factor out common monomial factors of the input polynomials.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_saturation";
  }

  int main()
  {
    FileParser P(Stdin);

    PolynomialSet a=P.parsePolynomialSetWithRing();

    AsciiPrinter(Stdout).printPolynomialRing(a.getRing());
    AsciiPrinter(Stdout).printNewLine();
    if(optionNoIdeal.getValue())
    {
    	a.saturate();
    	pout<<a;
    }
   	else
    {
    	AsciiPrinter(Stdout).printPolynomialSet(optionIsHomogeneous.getValue()?saturatedIdeal(a):nonHomogeneousSaturation(a));
    }
    return 0;
  }
};

static SaturationApplication theApplication;
