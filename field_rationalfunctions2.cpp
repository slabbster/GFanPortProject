#include "field_rationalfunctions2.h"

#include <memory>
#include <assert.h>
#include <gmp.h>
#include <sstream>
#include <iostream>

#include "termorder.h"
#include "division.h"
#include "buchberger.h"
#include "saturation.h"
#include "printer.h"

#include "log.h"

int FieldElementRationalFunctions2Living;

class FieldElementRationalFunction2 : public FieldElementImplementation
{
 public:
  Polynomial p,q;
  FieldElementRationalFunction2(FieldImplementation &a):
    FieldElementImplementation(a),
    p(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing()),
    q(Term(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing().getField().zHomomorphism(1),Monomial(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing())))
  {
    FieldElementRationalFunctions2Living++;
  }
  FieldElementRationalFunction2(FieldImplementation &a,int n_):
    FieldElementImplementation(a),
    p(Term(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing().getField().zHomomorphism(n_),Monomial(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing()))),
    q(Term(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing().getField().zHomomorphism(1),Monomial(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing())))
    {
      if(n_==0)p=Polynomial(((FieldRationalFunctions2Implementation*)&a)->getPolynomialRing());
      FieldElementRationalFunctions2Living++;
    }
  void normalize()
  {
/*    PolynomialSet g(p.getRing());
    g.push_back(p);
    g.push_back(q);
//    LexicographicTermOrder T;
    LexicographicTermOrder T;
    buchberger(&g,T);
    minimize(&g);
//    assert(g.size()==1);
    if(g.size()==1)
    {
    	Polynomial r=g.front();

    PolynomialSet Q(p.getRing());
    PolynomialSet P(p.getRing());
    division(p,g,T,&P);
    division(q,g,T,&Q);
    p=*P.begin();
    q=*Q.begin();
    assert(!q.isZero());
    FieldElement a=q.terms.rbegin()->second;
    p*=a;
    q*=a.inverse();
    }*/
    }
  FieldElementRationalFunction2(FieldImplementation &a, Polynomial const &p_, Polynomial const &q_):
    FieldElementImplementation(a),
    p(p_),
    q(q_)
    {

      FieldElementRationalFunctions2Living++;
    }
  virtual ~FieldElementRationalFunction2()
    {
      FieldElementRationalFunctions2Living--;
    }
  FieldElementRationalFunction2& operator=(const FieldElementRationalFunction2& a)
    {
      assert(0);
      return *this;
    }

  void operator*=(const FieldElementImplementation &a)
    {
      const FieldElementRationalFunction2 *A=(const FieldElementRationalFunction2*)&a;
      assert(A);

      p*=A->p;
      q*=A->q;
      normalize();
    }
  void operator+=(const FieldElementImplementation &a)
    {
      const FieldElementRationalFunction2 *A=(const FieldElementRationalFunction2*)&a;
      assert(A);

      p=p*A->q+A->p*q;
      q=A->q*q;
      normalize();
    }
  void madd(const FieldElementImplementation &a,const FieldElementImplementation &b)
    {
      const FieldElementRationalFunction2 *A=(const FieldElementRationalFunction2*)&a;
      const FieldElementRationalFunction2 *B=(const FieldElementRationalFunction2*)&b;
      assert(A);
      assert(B);

      p=p*(A->q*B->q)+(A->p*B->p)*q;
      q=A->q*B->q*q;
      normalize();
    }
  FieldElementRationalFunction2 *one() const;
  bool isZero() const
    {
      return p.isZero();
    }

  FieldElementRationalFunction2 *sum(const FieldElementImplementation &b)const
    {
      const FieldElementRationalFunction2 *B=(const FieldElementRationalFunction2*)&b;

      FieldElementRationalFunction2 *r= new FieldElementRationalFunction2(*getField(),p*B->q+B->p*q,B->q*q);

      return r;
    }
  FieldElementRationalFunction2 *difference(const FieldElementImplementation &b)const
    {
      const FieldElementRationalFunction2 *B=(const FieldElementRationalFunction2*)&b;
      FieldElementRationalFunction2 *r= new FieldElementRationalFunction2(*getField(),p*B->q-B->p*q,B->q*q);
      return r;
    }
  FieldElementRationalFunction2 *negation()const
    {
      FieldElementRationalFunction2 *r= new FieldElementRationalFunction2(*getField(),p-p-p,q);

      return r;
    }
  FieldElementImplementation *inverse()const
  {
    if(isZero())
      {
	AsciiPrinter P(Stderr);
	P.printString("Error inverting FieldElement: ");
	//	P.printFieldElement(*this);
	P.printString("\n");
	assert(0);
      }

    FieldElementRationalFunction2 *r= new FieldElementRationalFunction2(*getField(),q,p);

    return r;
  }

  int sign()const
  {
	  assert(0);//not an ordered field (yet)
    if(isZero())return 0;
    return p.terms.rbegin()->second.sign();
  }

  static string LaTeXTranslator(const string &s)
  {
	  assert(0);//not supported yet
/*    int startIndex=0;
    string sign;
    if(s[0]=='-')
      {
	sign=string("-");
	startIndex=1;
      }
    int slashIndex=-1;
    for(int i=startIndex;i<s.length();i++)if(s[i]=='/')slashIndex=i;
    if(slashIndex==-1)
      return string(s);

    return sign+string("{").append(s,startIndex,slashIndex-startIndex)+string("\\over ").append(s,slashIndex+1,s.length()-slashIndex-1)+string("}");
*/
	  }

  std::string toString(bool writeIfOne=true, bool alwaysWriteSign=false, bool latexMode=false /*, bool mathMode=true*/) const
  {
    stringstream s;

    //    s << "(" <<p.toString(latexMode) << "/" << q.toString(latexMode) << ")";
    s << "(";
    s<<p.toString(latexMode);
    s << "/";
    s << q.toString(latexMode);
    s << ")";
    return s.str();
  }

  FieldElementRationalFunction2 *copy()const
  {
    FieldElementRationalFunction2 *r= new FieldElementRationalFunction2(*getField());
    r->p=p;
    r->q=q;

    return r;
  }
};

PolynomialRing FieldRationalFunctions2Implementation::getPolynomialRing()const
{
  return thePolynomialRing;
}

bool FieldRationalFunctions2Implementation::isRationals()const
{
  return false;
}

FieldRationalFunctions2Implementation::FieldRationalFunctions2Implementation(PolynomialRing const &r):
  thePolynomialRing(r)
{
}

std::string FieldRationalFunctions2Implementation::toString()const
{
  stringstream s;
  s<< thePolynomialRing.getField().toString() << "("<<thePolynomialRing.toStringVariableNames()<< ")";
  return s.str();
}

FieldElementImplementation *FieldRationalFunctions2Implementation::zHomomorphismImplementation(int n)
{
  FieldElementImplementation *ret=new FieldElementRationalFunction2(*this,n);
  return ret;
}

FieldElement FieldRationalFunctions2Implementation::zHomomorphism(int n)
{
  return FieldElement(zHomomorphismImplementation(n));
}

const char *FieldRationalFunctions2Implementation::name()
{
  return "Rational functions in n variables";
}

FieldElementRationalFunction2 *FieldElementRationalFunction2::one() const
{
  return new FieldElementRationalFunction2(*getField(),1);
}



FieldRationalFunctions2::FieldRationalFunctions2(PolynomialRing const &r):
  Field(new FieldRationalFunctions2Implementation(r))
{
}


FieldElement FieldRationalFunctions2::polynomialToFraction(Polynomial const &p)
{
  FieldRationalFunctions2Implementation *imp=dynamic_cast<FieldRationalFunctions2Implementation*>(implementingObject);
  Polynomial q=Term(imp->getPolynomialRing().getField().zHomomorphism(1),Monomial(imp->getPolynomialRing()));

  return new FieldElementRationalFunction2(*imp, p, q);
}

/*****************************************************
 * Conversion functions
 *****************************************************/
PolynomialRing makeVariablesParameters(PolynomialRing const &r, int numberOfParameters)
{
	assert(numberOfParameters>=0);
	assert(numberOfParameters<=r.getNumberOfVariables());
	vector<string> names(numberOfParameters);
	for(int i=0;i<numberOfParameters;i++)names[i]=r.getVariableName(i);
	PolynomialRing temp(r.getField(),names);
	vector<string> names2(r.getNumberOfVariables()-numberOfParameters);
	for(int i=0;i<names2.size();i++)names2[i]=r.getVariableName(i+numberOfParameters);

	return PolynomialRing(FieldRationalFunctions2(temp),names2);
}

Polynomial makeVariablesParameters(PolynomialRing const &genericRing, Polynomial const &p)
{
	Polynomial ret(genericRing);
	FieldRationalFunctions2Implementation const *coefficientField=dynamic_cast<FieldRationalFunctions2Implementation const*>(genericRing.getField().implementingObject);
	FieldRationalFunctions2 &DANGER=(FieldRationalFunctions2&)genericRing.getField();
	PolynomialRing coefRing=coefficientField->getPolynomialRing();

	for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
	{
		FieldElement c=i->second;
		IntegerVector v=i->first.exponent;
		IntegerVector coefficientExponent=v.subvector(0,p.getRing().getNumberOfVariables()-genericRing.getNumberOfVariables());
		IntegerVector monomialExponent=v.subvector(p.getRing().getNumberOfVariables()-genericRing.getNumberOfVariables(),v.size());
		FieldElement c2=DANGER.polynomialToFraction(Term( c,Monomial(coefRing, coefficientExponent)));//does the numerator not belong to a field?
		ret+=Polynomial(Term(c2,Monomial(genericRing,monomialExponent)));
	}
	return ret;
}

PolynomialSet makeVariablesParameters(PolynomialRing const &genericRing, PolynomialSet const &p)
{
	PolynomialSet ret(genericRing);
	for(PolynomialSet::const_iterator i=p.begin();i!=p.end();i++)
		ret.push_back(makeVariablesParameters(genericRing,*i));

	return ret;
}
