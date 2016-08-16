#include "dimension.h"
#include "printer.h"
#include "parser.h"
#include "buchberger.h"
#include "gfanapplication.h"
#include "multiplicity.h"

class KrullDimensionApplication : public GFanApplication
{
  SimpleOption isGroebnerBasisOption;
  SimpleOption optionNumberOfSolutions;
public:
  KrullDimensionApplication():
    isGroebnerBasisOption("-g","Tell the program that the input is already a reduced Groebner basis."),
    optionNumberOfSolutions("--nsolutions","When the ideal is zero-dimensional this option will give the number of solutions over the algebraic closure of the coefficient field.")
  {
    optionNumberOfSolutions.hide();
    registerOptions();
  }
  const char *name()
  {
    return "_krulldimension";
  }
  int main()
  {
    FileParser P(Stdin);
    PolynomialSet g=P.parsePolynomialSetWithRing();
    if(!isGroebnerBasisOption.getValue())buchberger(&g,StandardGradedLexicographicTermOrder());
    int kd=krullDimension(g);
    pout<<kd<<"\n";
    if(optionNumberOfSolutions.getValue())
      {
        if(kd!=0){debug<<"Ideal is not zero dimensional!\n";assert(0);}
        pout<<numberOfStandardMonomials(g)<<"\n";
      }
    return 0;
  }
  const char *helpText()
  {
    return "Takes an ideal $I$ and computes the Krull dimension of R/I where R is the polynomial ring. This is done by first computing a Groebner basis.\n"
;
  }
};

static KrullDimensionApplication theApplication;
