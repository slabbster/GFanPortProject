#include "field_rationalfunctions.h"

#include <memory>
#include <assert.h>
#include <cstdio> /* Always include cstdio before gmp.h.*/
#include <gmp.h>
#include <sstream>
#include <iostream>

#include "termorder.h"
#include "division.h"
#include "buchberger.h"
#include "printer.h"

#include "log.h"

int FieldElementRationalFunctionsLiving;

class FieldElementRationalFunction : public FieldElementImplementation
{
 public:
  Polynomial p,q;
  FieldElementRationalFunction(FieldImplementation &a):
    FieldElementImplementation(a),
    p(((FieldRationalFunctionsImplementation*)&a)->getPolynomialRing()),
    q(Term(((FieldRationalFunctionsImplementation*)&a)->getPolynomialRing().getField().zHomomorphism(1),Monomial(((FieldRationalFunctionsImplementation*)&a)->getPolynomialRing(),IntegerVector(1))))
  {
    FieldElementRationalFunctionsLiving++;
  }
  FieldElementRationalFunction(FieldImplementation &a,int n_):
    FieldElementImplementation(a),
    p(Term(((FieldRationalFunctionsImplementation*)&a)->getPolynomialRing().getField().zHomomorphism(n_),Monomial(((FieldRationalFunctionsImplementation*)&a)->getPolynomialRing(),IntegerVector(1)))),
    q(Term(((FieldRationalFunctionsImplementation*)&a)->getPolynomialRing().getField().zHomomorphism(1),Monomial(((FieldRationalFunctionsImplementation*)&a)->getPolynomialRing(),IntegerVector(1))))
    {
      if(n_==0)p=Polynomial(((FieldRationalFunctionsImplementation*)&a)->getPolynomialRing());
      /* fprintf(stderr,"constructing\n");
      AsciiPrinter(Stderr).printPolynomial(p);
      AsciiPrinter(Stderr).printPolynomial(q);
      fprintf(stderr,"\n");*/
      FieldElementRationalFunctionsLiving++;
    }
  void normalize()
  {
    PolynomialSet g(p.getRing());
    g.push_back(p);
    g.push_back(q);
    LexicographicTermOrder T;
    buchberger(&g,T);
    minimize(&g);
    assert(g.size()==1);
    Polynomial r=g.front();
    PolynomialSet Q(p.getRing());
    PolynomialSet P(p.getRing());
    division(p,g,T,&P);
    division(q,g,T,&Q);
    p=*P.begin();
    q=*Q.begin();
    assert(!q.isZero());
    //    FieldElement a=q.terms.begin()->second;
    //    TermMap temp=q.terms.end();
    FieldElement a=q.terms.rbegin()->second;
    p*=a;
    q*=a.inverse();
    //    assert(0);//make q's leading coefficient 1
  }
  FieldElementRationalFunction(FieldImplementation &a, Polynomial const &p_, Polynomial const &q_):
    FieldElementImplementation(a),
    p(p_),
    q(q_)
    {

      FieldElementRationalFunctionsLiving++;
    }
  /*  FieldElementRationalFunction(FieldImplementation &a, mpq_t *n_):FieldElementImplementation(a)
  {
    FieldElementRationalsLiving++;
    mpq_init(value);
    mpq_set(value,*n_);
    }*/
  virtual ~FieldElementRationalFunction()
    {
      FieldElementRationalsLiving--;
      //      mpq_clear(value);
    }
  FieldElementRationalFunction& operator=(const FieldElementRationalFunction& a)
    {
      assert(0);
      /*      const FieldElementRational *A=(const FieldElementRational*)&a;
      if (this != A) {
        mpq_clear(value);
        mpz_init_set(mpq_numref(value), mpq_numref(a.value));
        mpz_init_set(mpq_denref(value), mpq_denref(a.value));
      }
      */
      return *this;
    }

  void operator*=(const FieldElementImplementation &a)
    {
      const FieldElementRationalFunction *A=(const FieldElementRationalFunction*)&a;
      assert(A);

      //      TimerScope ts(&rationalTimer);
      p*=A->p;
      q*=A->q;
      normalize();
      //      mpq_mul(value,value,A->value);
    }
  void operator+=(const FieldElementImplementation &a)
    {
      const FieldElementRationalFunction *A=(const FieldElementRationalFunction*)&a;
      assert(A);

      p=p*A->q+A->p*q;
      q=A->q*q;
      normalize();
    }
  void madd(const FieldElementImplementation &a,const FieldElementImplementation &b)
    {
      const FieldElementRationalFunction *A=(const FieldElementRationalFunction*)&a;
      const FieldElementRationalFunction *B=(const FieldElementRationalFunction*)&b;
      assert(A);
      assert(B);

      p=p*(A->q*B->q)+(A->p*B->p)*q;
      q=A->q*B->q*q;
      normalize();
    }

  FieldElementRationalFunction *one() const;
  bool isZero() const
    {
      return p.isZero();
    }

  FieldElementRationalFunction *sum(const FieldElementImplementation &b)const
    {
      const FieldElementRationalFunction *B=(const FieldElementRationalFunction*)&b;

      FieldElementRationalFunction *r= new FieldElementRationalFunction(*getField(),p*B->q+B->p*q,B->q*q);

      return r;
    }
  FieldElementRationalFunction *difference(const FieldElementImplementation &b)const
    {
      const FieldElementRationalFunction *B=(const FieldElementRationalFunction*)&b;
      FieldElementRationalFunction *r= new FieldElementRationalFunction(*getField(),p*B->q-B->p*q,B->q*q);
      return r;
    }
  FieldElementRationalFunction *negation()const
    {
      FieldElementRationalFunction *r= new FieldElementRationalFunction(*getField(),p-p-p,q);

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

    FieldElementRationalFunction *r= new FieldElementRationalFunction(*getField(),q,p);

    return r;
  }

  int sign()const
  {
    if(isZero())return 0;
    return p.terms.rbegin()->second.sign();
  }

  static string LaTeXTranslator(const string &s)
  {
    int startIndex=0;
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
    /*
    s<<"(";
    FieldElementRational *tempOne=one();
    FieldElementRational *temp=difference(*tempOne);
    bool isOne=temp->isZero();
    //fprintf(Stderr,"DELETE\n");
    delete temp;
    temp=sum(*tempOne);
    bool isMinusOne=temp->isZero();
    //fprintf(Stderr,"DELETE\n");
    delete temp;
    //fprintf(Stderr,"DELETE\n");
    delete tempOne;

    if(!writeIfOne && isOne)
      {
	if(alwaysWriteSign)return std::string("+");
	return std::string("");
      }
    if(!writeIfOne && isMinusOne)
      return std::string("-");
    static char s[1290*1000];
    //    mpq_get_str(s,10,value); //// CHECK BUFFER SIZE!!!!
    // Changed to make code gmp 3.1.1 compatible
    mpz_get_str(s,10,mpq_numref(value)); //CHECK BUFFER SIZE!!!!
    string S(s);
    if(mpz_cmp_ui(mpq_denref(value),1)!=0)
      {
	mpz_get_str(s,10,mpq_denref(value)); //CHECK BUFFER SIZE!!!!
	S=S+string("/")+string(s);
      }

    if(latexMode)S=LaTeXTranslator(S);

    if(alwaysWriteSign && mpq_sgn(value)!=-1)
      return std::string("+")+S;
      return S;
    */
    return s.str();
  }

  FieldElementRationalFunction *copy()const
  {
    FieldElementRationalFunction *r= new FieldElementRationalFunction(*getField());
    //      fprintf(Stderr,"NEW\n");
    /*
    mpq_clear(r->value);
    mpz_init_set(mpq_numref(r->value), mpq_numref(value));
    mpz_init_set(mpq_denref(r->value), mpq_denref(value));
    */
    r->p=p;
    r->q=q;

    return r;
  }
};

PolynomialRing FieldRationalFunctionsImplementation::getPolynomialRing()const
{
  return thePolynomialRing;
}

bool FieldRationalFunctionsImplementation::isRationals()const
{
  return false;
}

FieldRationalFunctionsImplementation::FieldRationalFunctionsImplementation(Field const &f_, string const &parameterName_):
  coefficientField(f_),
  parameterName(parameterName_),
  thePolynomialRing(f_,1)
{
  vector<string> l;
  l.push_back(parameterName_);
  thePolynomialRing= PolynomialRing(f_,l);
}

std::string FieldRationalFunctionsImplementation::toString()const
{
  stringstream s;
  s<< coefficientField.toString() << "(" << parameterName << ")";
  return s.str();
}

FieldElementImplementation *FieldRationalFunctionsImplementation::zHomomorphismImplementation(int n)
{
  /*  if(n==0)
    {
      static FieldElementImplementation *p;
      if(p==0)p=new FieldElementRational(*this,0);
      p->refCount++;
      return p;
    }
  else
    if(n==1)
      {
	static FieldElementImplementation *p;
	if(p==0)p=new FieldElementRational(*this,1);
	p->refCount++;
	return p;
      }
    else
      if(n==-1)
	{
	  static FieldElementImplementation *p;
	  if(p==0)p=new FieldElementRational(*this,-1);
	  p->refCount++;
	  return p;
	}
  */
  FieldElementImplementation *ret=new FieldElementRationalFunction(*this,n);
  //      fprintf(Stderr,"NEW\n");
  //  ret->refCount++;
  return ret;
}

FieldElement FieldRationalFunctionsImplementation::zHomomorphism(int n)
{
  //      fprintf(Stderr,"NEW\n");
  //  return FieldElement(new FieldElementRational(*this,n));
  return FieldElement(zHomomorphismImplementation(n));
}

const char *FieldRationalFunctionsImplementation::name()
{
  return "Rational functions";
}

FieldElementRationalFunction *FieldElementRationalFunction::one() const
{
  //      fprintf(Stderr,"NEW\n");
  return new FieldElementRationalFunction(*getField(),1);
}



FieldRationalFunctions::FieldRationalFunctions(Field const &coefficientField, string const &parameterName):
  Field(new FieldRationalFunctionsImplementation(coefficientField,parameterName))
{
}


FieldElement FieldRationalFunctions::exponent(int power)
{
  FieldRationalFunctionsImplementation *imp=dynamic_cast<FieldRationalFunctionsImplementation*>(implementingObject);
  IntegerVector v(1);
  v[0]=power;
  Polynomial p=Term(imp->getPolynomialRing().getField().zHomomorphism(1),Monomial(imp->getPolynomialRing(),v));
  Polynomial q=Term(imp->getPolynomialRing().getField().zHomomorphism(1),Monomial(imp->getPolynomialRing(),IntegerVector(1)));

  return new FieldElementRationalFunction(*imp, p, q);
}


FieldElement FieldRationalFunctions::fromCoefficientField(FieldElement const &c)
{
  FieldRationalFunctionsImplementation *imp=dynamic_cast<FieldRationalFunctionsImplementation*>(implementingObject);
  IntegerVector v(1);
  Polynomial p=Term(c,Monomial(imp->getPolynomialRing(),v));
  Polynomial q=Term(imp->getPolynomialRing().getField().zHomomorphism(1),Monomial(imp->getPolynomialRing(),IntegerVector(1)));

  return new FieldElementRationalFunction(*imp, p, q);
}



FieldElement FieldRationalFunctions::substitute(FieldElement const &e, FieldElement const &tvalue)const
{
  FieldElementRationalFunction *eimp=dynamic_cast<FieldElementRationalFunction*>(e.implementingObject);
  FieldElement p=eimp->p.evaluate(tvalue);
  FieldElement q=eimp->q.evaluate(tvalue);

  assert(!q.isZero());

  return p*q.inverse();
}

#include "field_rationals.h"
#include <iostream>
void testRationalFunctionField()
{
    FieldRationalFunctions F(Q,string("t"));
    AsciiPrinter(Stderr).printField(F);


    FieldElement a=F.zHomomorphism(2);

    FieldElement b=F.exponent(3);

    a=b*a;
    AsciiPrinter(Stderr).printPolynomial(((dynamic_cast<FieldElementRationalFunction*>((a.implementingObject)))->p));

    /*    //    cerr<< ((FieldElementRationalFunction*)&(a.implementingObject))-> p.toString(false);
    cerr<< (dynamic_cast<FieldElementRationalFunction*>(a.implementingObject))-> p.toString(false);
    */


    cerr<<a.implementingObject->toString();
}
