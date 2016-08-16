#include "singularconversion.h"
#include <assert.h>

const  char *singular_date=__DATE__" "__TIME__;

ring singularRing(PolynomialRing const &r)
{
  //  FieldRationalsImplementation *q=dynamic_cast<FieldRationalsImplementation>r.getField().implementingObject;

  ring ret=(ring)omAlloc0(sizeof(sip_sring));

  if(r.getField().isRationals())
    {
      ret->ch=0;
    }
  else
    {
      assert(0);
    }

  ret->N=r.getNumberOfVariables();

  ret->names=(char**) omAlloc(ret->N*sizeof(char*));
  for(int i=0;i<ret->N;i++)
    ret->names[i]=omStrDup(r.getVariableName(i).c_str());

  ret->order=(int*)omAlloc0(3*sizeof(int));
  ret->block0=(int*)omAlloc0(3*sizeof(int));
  ret->block1=(int*)omAlloc0(3*sizeof(int));

  //  ret->order[0]=ringorder_wp;//degree revlex
  ret->order[0]=ringorder_Wp;//degree lex
  ret->block0[0]=1;ret->block1[0]=ret->N;
  ret->order[1]=ringorder_C;

  ret->wvhdl=(int**)omAlloc0(3*sizeof(int*));
  ret->wvhdl[0]=(int*)omAlloc(ret->N*sizeof(int));
  for(int i=0;i<ret->N;i++)
    ret->wvhdl[0][i]=1;

  rComplete(ret);
  rChangeCurrRing(ret);

  return ret;
}


void freeSingularRing(ring R)
{
  rKill(R);
}


poly singularPolynomial(Polynomial const &p)
{
  int n=p.getRing().getNumberOfVariables();
  poly r=NULL;

  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      poly p=pInit();

      FieldElement c=i->second;

      mpq_t *C=fieldElementToGmp(c);
      number n2=nlRInit(0);
      mpz_init_set(&(n2)->z,mpq_numref(*C));
      mpz_init_set(&(n2)->n,mpq_denref(*C));
      n2->s=0;
      nNormalize(n2);
      pSetCoeff0(p,n2);

      for(int j=0;j<n;j++)
	pSetExp(p,j+1,i->first.exponent[j]);
      pSetComp(p,0);
      pSetm(p);
      r=pAdd(r,p);
    }

  return r;
}


ideal singularPolynomialSet(PolynomialSet const &g)
{
  int m=g.size();


  ideal i=idInit(m,1);

  int J=0;
  for(PolynomialSet::const_iterator j=g.begin();j!=g.end();j++,J++)
    {
      i->m[J]=singularPolynomial(*j);
    }

  return i;
}


FieldElement fromSingularCoefficient(PolynomialRing const &r, number c)
{
  FieldElement C(r.getField());

  //  if((nInt(c)==0) && ((SR_HDL(c) & SR_INT)==0))
  if(((SR_HDL(c) & SR_INT)==0))
    {
      mpq_t value;

      mpq_init(value);

      mpz_set(mpq_numref(value), &(c->z));
      if (c->s <3)
      mpz_set(mpq_denref(value), &(c->n));
      else
	 mpz_set_si(mpq_denref(value),1);

      C=fieldElementFromGmp(&value);
    }
  else
    {
      C=r.getField().zHomomorphism(nInt(c));
    }
  return C;
}


Polynomial fromSingularPolynomial(PolynomialRing const &r, poly &p)
{
  int n=r.getNumberOfVariables();
  Polynomial ret(r);

  poly q=p;
  while(p)
    {
      IntegerVector v(n);
      for(int i=0;i<n;i++)
	v[i]=pGetExp(p,i+1);

      ret+=Term(fromSingularCoefficient(r,pGetCoeff(p)),Monomial(r,v));
      p=pNext(p);
    }
  if(q)
    {
      IntegerVector v(n);
      for(int i=0;i<n;i++)
	v[i]=pGetExp(q,i+1);

      ret.mark(Monomial(r,v));
    }

  return ret;
}


PolynomialSet fromSingularIdeal(PolynomialRing const &r, ideal i)
{
  PolynomialSet ret(r);

  for(int j=0;j<IDELEMS(i);j++)
    if(i->m[j]!=NULL)
      ret.push_back(fromSingularPolynomial(r,i->m[j]));

  //  AsciiPrinter(Stderr).printPolynomialSet(ret);

  return ret;
}
