#include "latticeideal.h"
#include "binomial.h"
#include "buchberger.h"
#include "printer.h"
#include "field_rationals.h"

IntegerVectorList latticeIdealRevLex(IntegerMatrix const &lattice)
{
  PolynomialRing theRing(Q,lattice.getHeight());
  IntegerMatrix L=lattice.transposed();
  PolynomialSet g(theRing);

  for(int i=0;i<L.getHeight();i++)
    g.push_back(integerVectorToBinomial(L[i],theRing));

  //  AsciiPrinter(Stderr).printPolynomialSet(g);


  for(int i=0;i<L.getWidth();i++)
    {
      ReverseLexicographicTermOrder t(i);
      buchberger(&g,t);
      for(PolynomialSet::iterator j=g.begin();j!=g.end();j++)
	j->saturate();
    }
  IntegerVectorList ret;
  for(PolynomialSet::const_iterator j=g.begin();j!=g.end();j++)
    ret.push_back(binomialToIntegerVector(*j));

  return ret;
}
