#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "buchberger.h"
#include "wallideal.h"
#include "termorder.h"
#include "gfanapplication.h"
#include "tropical2.h"
#include "matrix.h"

class TropicalEvaluationApplication : public GFanApplication
{
public:
  const char *helpText()
  {
    return "This program evaluates a tropical polynomial function in a given set of points.\n";
  }
  TropicalEvaluationApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_tropicalevaluation";
  }

  int main()
  {
    FileParser P(Stdin);

    PolynomialSet g=P.parsePolynomialSetWithRing();
    IntegerVectorList w=P.parseIntegerVectorList();
    IntegerMatrix m=rowsToIntegerMatrix(w,g.getRing().getNumberOfVariables());

    IntegerMatrix result(m.getHeight(),g.size());

    for(int i=0;i<result.getHeight();i++)
      {
	PolynomialSet::const_iterator J=g.begin();
	for(int j=0;j<result.getWidth();j++,J++)
	  result[i][j]=J->degree(m[i]);
      }

    AsciiPrinter(Stdout).printVectorList(result.getRows());
    return 0;
  }
};

static TropicalEvaluationApplication theApplication;
