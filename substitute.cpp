#include "substitute.h"

Polynomial multiplicativeChange(Polynomial const &p, IntegerMatrix const &mat)
{
  PolynomialRing theRing=p.getRing();
  Polynomial ret(theRing);
  if(!p.isZero())
    {
      IntegerVector rel=p.terms.begin()->first.exponent;
      for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
	ret+=Term(i->second,Monomial(theRing,mat.vectormultiply(i->first.exponent)));
    }
  return ret;
}


PolynomialSet multiplicativeChange(PolynomialSet const &g, IntegerMatrix const &mat)
{
  PolynomialRing theRing=g.getRing();
  PolynomialSet ret(theRing);
  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    ret.push_back(multiplicativeChange(*i,mat));

  return ret;
}
