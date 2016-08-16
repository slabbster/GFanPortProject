#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "minors.h"
#include "field_rationals.h"
#include "tropical2.h"
#include "matrix.h"
#include "buchberger.h"
#include "tropicaldeterminant.h"

class TropicalRankApplication : public GFanApplication
{
  SimpleOption optionKapranov;
  SimpleOption optionDeterminant;
public:
  const char *helpText()
  {
    return "This program will compute the tropical rank of matrix given as input. Tropical addition is MAXIMUM.\n";
  }
  TropicalRankApplication():
    optionKapranov("--kapranov","Compute Kapranov rank instead of tropical rank."),
    optionDeterminant("--determinant","Compute the tropical determinant instead.")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_tropicalrank";
  }
  int main()
  {
    FileParser P(Stdin);

    IntegerMatrix m=rowsToIntegerMatrix(P.parseIntegerVectorList());

    if(optionDeterminant.getValue())
    {
//    	tropicalDeterminantTest();
    	pout<<tropicalDeterminant(m);
    	pout<<"\n";
    	return 0;
    }

    IntegerVector w=flattenMatrix(m);
    int d=m.getHeight();
    int n=m.getWidth();
    int theRank=-1;
    int min=d>n?n:d;
    for(int r=1;r<=min;r++)
      {
	PolynomialRing R(Q,matrixVariableNames("m",d,n));

	PolynomialSet p=minors(R,r,d,n,false,false);

	if(optionKapranov.getValue())
	  {
	    WeightReverseLexicographicTermOrder T(w);
	    buchberger(&p,T);
	  }

	PolynomialSet q=initialForms(p,w);
	bool containsMonomial=false;
	for(PolynomialSet::const_iterator i=q.begin();i!=q.end();i++)
	  if(i->isMonomial())
	    {
	      containsMonomial=true;
	      break;
	    }
	fprintf(Stderr,"%ix%i picks monomial: %s %i\n",r,r,containsMonomial?"true, meaning that rank is >":"false, meaning that rank is <=",r-1);
	if(!containsMonomial)
	  {
	    theRank=r-1;
	    break;
	  }
      }
    if(theRank==-1)
      theRank=min;
    AsciiPrinter(Stdout).printInteger(theRank);
    AsciiPrinter(Stdout).printNewLine();

   return 0;
  }
};

static TropicalRankApplication theApplication;
