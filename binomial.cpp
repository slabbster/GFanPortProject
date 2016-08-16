#include "binomial.h"
#include "field_zmodpz.h"

IntegerVector binomialToIntegerVector(Polynomial const &p)
{
  if(p.isZero())
    {
      assert(0); // There is no way of finding the dimension!!
      return IntegerVector(0);
    }
  assert(p.numberOfTerms()==2);

  assert(p.isMarked());

  IntegerVector markedExponent=p.getMarked().m.exponent;

  IntegerVector ret;
  for(TermMap::const_iterator j=p.terms.begin();j!=p.terms.end();j++)
    {
      IntegerVector dif=markedExponent-j->first.exponent;
      
      if(!dif.isZero())
	{
	  ret=dif;
	}
    }
  return ret;
}


Polynomial integerVectorToBinomial(IntegerVector const &v, PolynomialRing const &r)
{
  IntegerVector pos=max(v,IntegerVector(v.size()));
  IntegerVector neg=pos-v;
  Field field=r.getField();
  //  if(!field)field=theZMod2ZField();//should this be changed to the global field??
  Term a(field.zHomomorphism(1),Monomial(r,pos));
  Term b(field.zHomomorphism(-1),Monomial(r,neg));
  Polynomial p=Polynomial(a)+Polynomial(b);
  p.mark(Monomial(r,pos));

  return p;
}
