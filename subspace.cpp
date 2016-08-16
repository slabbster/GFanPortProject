#include "subspace.h"

#include "division.h"
#include "buchberger.h"
#include "field_rationals.h"


static Polynomial vectorToPolynomial(PolynomialRing const &r, IntegerVector const &v)
{
  /*static Field* field;

  if(!field)field=Field::find("GmpRationals"); // this is a bit stupid. We should add a field_rationals header file instead
  assert(field);
  */
  Polynomial ret(r);

  for(int i=0;i<v.size();i++)
    if(v[i])
      {
	ret+=Term(r.getField().zHomomorphism(v[i]),Monomial(r,IntegerVector::standardVector(v.size(),i)));
      }

  return ret;
}


static IntegerVector polynomialToVector(Polynomial const &f)
{
  /*  static Field* field;

  if(!field)field=Field::find("GmpRationals"); // this is a bit stupid. We should add a field_rationals header file instead
  assert(field);
  */
  int n=f.getNumberOfVariables();

  //  fprintf(Stderr,"%i\n",n);
  
  vector<FieldElement> r(n);

  for(int i=0;i<n;i++)
    r[i]=f.getRing().getField().zHomomorphism(0);

  for(TermMap::const_iterator i=f.terms.begin();i!=f.terms.end();i++)
    {
      for(int j=0;j<n;j++)
	if(i->first.exponent==IntegerVector::standardVector(n,j))
	  {
	    r[j]=i->second;
	  }
    }

  return primitiveVector(r);
}

Subspace::Subspace(IntegerVectorList const &generators, int ambientDimension):
  basis2(integerMatrixToFieldMatrix(rowsToIntegerMatrix(generators,ambientDimension),Q))
		  //  ,  basis(PolynomialRing(Q,ambientDimension))
{
  //  n=ambientDimension;
  n=basis2.getWidth();
  /*  {
    PolynomialRing theRing=basis.getRing();
    if(n==-1)
      {
	assert(!generators.empty());
	n=generators.begin()->size();
      }
    for(IntegerVectorList::const_iterator i=generators.begin();i!=generators.end();i++)
      {
	assert(i->size()==n);
	basis.push_back(vectorToPolynomial(theRing,*i));
      }
    buchberger(&basis,LexicographicTermOrder());
    }*/
  
  basis2.reduce();
  basis2.REformToRREform();
  basis2.scaleFirstNonZeroEntryToOne();
  basis2.removeZeroRows();
  integerRep=getRepresentation();
  integerRep.sort();
}

bool Subspace::contains(IntegerVector const &v)const
{
  if(v.size()!=n)
    {
      fprintf(Stderr,"v.size()=%i, n=%i",v.size(),n);
    }
  assert(v.size()==n);

  //  bool oldRet=division(vectorToPolynomial(PolynomialRing(Q,v.size()),v),basis,LexicographicTermOrder()).isZero();

  bool ret=basis2.canonicalize(integerVectorToFieldVector(v,Q)).isZero();

  //  assert(oldRet==ret);

  return ret;
}


int Subspace::dimension()
{
  //    assert(basis2.reduceAndComputeRank()==basis.size());
  return basis2.reduceAndComputeRank();
  //  return basis.size();
}


int Subspace::ambientDimension()const
{
  return n;
}


Subspace sum(Subspace const &a, Subspace const &b)
{
  /*Subspace ret=a;

      ret.basis.insert(ret.basis.end(),b.basis.begin(),b.basis.end());
  buchberger(&ret.basis,LexicographicTermOrder());
  
    return ret;
  */
   Subspace ret=a;
  ret.basis2=combineOnTop(a.basis2,b.basis2);
  ret.basis2.reduce();
  ret.basis2.REformToRREform();
  ret.basis2.scaleFirstNonZeroEntryToOne();
  ret.basis2.removeZeroRows();

  ret.integerRep=ret.getRepresentation();
  ret.integerRep.sort();

  return ret;
}


IntegerVectorList Subspace::getRepresentation()const
{
  /*IntegerVectorList ret;

      for(PolynomialSet::const_iterator i=basis.begin();i!=basis.end();i++)
    {
      ret.push_back(polynomialToVector(*i));
    }
  */
  IntegerVectorList ret2;
  for(int i=0;i<basis2.getHeight();i++)
    {
      if(!basis2[i].isZero())
	{
	  /*	  int sign=0;*/
	  IntegerVector temp=basis2[i].primitive();
	  /*  for(int j=0;j<n;j++)if(temp[j]!=0){sign=temp[j];break;}
	  if(sign>0)*/
	    ret2.push_back(temp);
	    /*	  else
	    ret2.push_back(-temp);
	    */	}
    }

  /*  AsciiPrinter(Stderr).printVectorList(ret);
  AsciiPrinter(Stderr).printVectorList(ret2);

  fprintf(Stderr,"------------------\n");
  */
  return ret2;
}

IntegerVector Subspace::canonicalizeVector(IntegerVector const &v)const
{
  /*  Polynomial f=vectorToPolynomial(PolynomialRing(Q,v.size()),v);
  f=division(f,basis,LexicographicTermOrder());

  IntegerVector ret=polynomialToVector(f);
  */
  IntegerVector ret2=basis2.canonicalize(integerVectorToFieldVector(v,Q)).primitive();

  //  assert((ret-ret2).isZero());

  return ret2;
}


bool Subspace::operator<(Subspace const &b)const
{
  return integerRep<b.integerRep;
}
