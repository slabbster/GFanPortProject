#include "saturation.h"
#include "matrix.h"
#include "multiplicity.h"
#include "buchberger.h"
#include "division.h"
#include "tropical.h"
#include "printer.h"

/*
#include "substitution.h"

PolynomialSet colonIdeal(PolynomialSet const &ideal, polynomial f)
{
  int n=ideal.numberOfVariablesInRing();

  IntegerMatrix mat(n+1,n);
  for(int i=0;i<n;i++)mat[i][i]=1;

  PolynomialSet ret=multiplicativeChange(mat,ideal);
  f=multiplicativeChange(mat,f);

  f*=Monomial(IntegerVector::standardVector(n+1,n));

  ret.push_back();
}
*/

FieldElement getFieldElement(PolynomialSet const &s)
{
  FieldElement a;
  for(PolynomialSet::const_iterator i=s.begin();i!=s.end();i++)
    for(TermMap::const_iterator j=i->terms.begin();j!=i->terms.end();j++)
      return j->second;
  assert(0);
  return a;
}

/*
 * This procedure does not work
 */
static PolynomialSet singleSaturation(PolynomialSet const &s)//TO DO: call idealIntersection to do the intersection
{
	assert(0);
	PolynomialRing theRing=s.getRing();
  int n=s.numberOfVariablesInRing();

  FieldElement one=getFieldElement(s).one();

  // add a variable

  //  PolynomialSet s2=s;s2.changeNumberOfVariables(n+1);

  IntegerVector temp(theRing.getNumberOfVariables());
  PolynomialSet s2=s.homogenization(theRing.withVariablesAppended("T"),&temp);

  /*  IntegerMatrix mat(n+1,n);
  for(int i=0;i<n;i++)mat[i][i]=1;
  PolynomialSet s2=multiplicativeChangeInv(s,mat);*/

  // multiply by (1-t)
  Monomial m(s2.getRing(),IntegerVector::standardVector(n+1,n));

  PolynomialSet s3(s2.getRing());
  for(PolynomialSet::const_iterator i=s2.begin();i!=s2.end();i++)
    {
      s3.push_back(*i-((*i)*m));
    }

  for(int i=0;i<n;i++)
    s3.push_back(Polynomial(Term(one,Monomial(s3.getRing(),IntegerVector::standardVector(n+1,n)+IntegerVector::standardVector(n+1,i)))));

  //  WeightReverseLexicographicTermOrder T(IntegerVector::standardVector(n+1,n));
  LexicographicTermOrder T(n);
  buchberger(&s3,T);

  PolynomialSet s4=s3.polynomialRingIntersection(theRing);
  s4.saturate();

  return s4;
}

/*
 *
 */
PolynomialSet nonHomogeneousSaturation(PolynomialSet const &s)
{
	PolynomialSet s2=s.homogenization(s.getRing().withVariablesAppended("T"));
	PolynomialSet s3=saturatedIdeal(s2);
	PolynomialSet s4=s3.deHomogenizationInSameRing();
	return s4.polynomialRingIntersection(s.getRing());
	/*  PolynomialRing theRing=s.getRing();
  StandardGradedLexicographicTermOrder T;
  PolynomialSet s1=s;
  PolynomialSet s2(theRing);
  buchberger(&s1,T);
  s1.sort_();

  do
    {
      s2=s1;
      s1=singleSaturation(s2);

      buchberger(&s1,T);
      s1.sort_();
    }
  while(!s1.isEqualTo(s2));

  return s1;
*/
}


PolynomialSet idealIntersection(PolynomialSet const &a, PolynomialSet const &b)
{
  PolynomialRing theRing=a.getRing();
  int n=a.numberOfVariablesInRing();
  assert(n==b.numberOfVariablesInRing());

  FieldElement one=getFieldElement(a).one();

  // add a variable
  PolynomialRing theRing2=theRing.withVariablesAppended("T");
  IntegerVector temp(theRing.getNumberOfVariables());
  PolynomialSet A=a.homogenization(theRing2,&temp);
  PolynomialSet B=b.homogenization(theRing2,&temp);
  //  PolynomialSet A=a;A.changeNumberOfVariables(n+1);
  //  PolynomialSet B=b;B.changeNumberOfVariables(n+1);


  // multiply by (1-t)
  Monomial m(theRing2,IntegerVector::standardVector(n+1,n));

  PolynomialSet C(theRing2);
  for(PolynomialSet::const_iterator i=B.begin();i!=B.end();i++)
    C.push_back(*i-((*i)*m));

  for(PolynomialSet::const_iterator i=A.begin();i!=A.end();i++)
    C.push_back(((*i)*m));

  LexicographicTermOrder T(n);
  buchberger(&C,T);

  return C.polynomialRingIntersection(theRing);
}
