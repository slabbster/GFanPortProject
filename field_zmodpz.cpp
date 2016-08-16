#include "field.h"

#include "field_zmodpz.h"
#include <memory>
#include <assert.h>

#include "printer.h"

static int modP(int a, int p)
{
  return ((a%p)+p)%p;
}

static void swap(int &a,int &b)
{
  int temp=a;
  a=b;
  b=temp;
}

static int gcd(int a, int b, int &r, int &s)
{
  int A=a;
  int B=b;
  assert(a>0);
  assert(b>0);
  int aA=1;
  int aB=0;
  int bA=0;
  int bB=1;
  while(0!=b)
    {
      //      fprintf(Stderr,"a:%i b:%i aA:%i aB:%i bA:%i bB:%i\n",a,b,aA,aB,bA,bB);
      assert(a==aA*A+aB*B);
      assert(b==bA*A+bB*B);
      if(a>b)
	{
	  swap(a,b);
	  swap(aA,bA);
	  swap(aB,bB);
	}
      int n=b/a;
      assert(n!=0);
      b-=n*a;
      bA-=n*aA;
      bB-=n*aB;
    }
  //  fprintf(Stderr,"a:%i b:%i aA:%i aB:%i bA:%i bB:%i\n",a,b,aA,aB,bA,bB);
  assert(a==aA*A+aB*B);
  assert(b==bA*A+bB*B);
  r=aA;
  s=aB;
  return a;
}

const char *FieldZModPZImplementation::name()
{
  static char s[20];
  sprintf(s,"Zmod%iZ",p);
  return s;
  // return "Zmod2Z";
}


static bool isPrime(int p)
{
  for(int i=2;i*i<=p;i++)
    {
      if(p%i ==0)return false;
    }
  return true;
}



FieldZModPZImplementation::FieldZModPZImplementation(int p_):
  p(p_)
{
  if(p>=32768 || p<2)
    {
      fprintf(Stderr,"Prime %i out of range!\n",p_);
      assert(0);
    }
  if(!isPrime(p_))
    {
      fprintf(Stderr,"%i is not a prime.\n",p_);
      assert(0);
    }
}

FieldZModPZ::FieldZModPZ(int p_):
  Field(new FieldZModPZImplementation(p_))
{
}


int FieldZModPZ::getP()
{
  return ((FieldZModPZImplementation*)(implementingObject))->getP();
}


int FieldZModPZImplementation::getP()const
{
  return p;
}

std::string FieldZModPZImplementation::toString()const
{
  //  return std::string("Z/")+string(p)+string("Z");
  //  return std::string("Z/")+string("Z");
  static char s[20];
  sprintf(s,"Z/%iZ",p);
  return std::string(s);
}

class FieldElementZModPZ : public FieldElementImplementation
{
  int value;
  //  class FieldZModPZImplementation &theFieldImplementation;
 public:
  FieldZModPZImplementation &getF()const
  {
    return *(FieldZModPZImplementation*)getField();
  }
  FieldElementZModPZ(class FieldZModPZImplementation &theField_):
    FieldElementImplementation(theField_),
    value(0)
    {
    }
  FieldElementZModPZ(class FieldZModPZImplementation &theField_, int n_):
    FieldElementImplementation(theField_)
    {
      value=modP(n_,getF().getP());
    }
  virtual ~FieldElementZModPZ() //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    {
    }
  FieldElementZModPZ& operator=(const FieldElementZModPZ& a)
    {
      assert(0);
      const FieldElementZModPZ *A=(const FieldElementZModPZ*)&a;
      if (this != A) {
	value=A->value;
      }
      return *this;
    }

  void operator*=(const FieldElementImplementation &a)
    {
      const FieldElementZModPZ *A=(const FieldElementZModPZ*)&a;
      assert(A);
      value=modP(value*A->value,getF().getP());
    }
  void operator+=(const FieldElementImplementation &a)
    {
      const FieldElementZModPZ *A=(const FieldElementZModPZ*)&a;
      assert(A);
      value=modP(value+(A->value),getF().getP());
    }
  void madd(const FieldElementImplementation &a,const FieldElementImplementation &b)
    {
      const FieldElementZModPZ *A=(const FieldElementZModPZ*)&a;
      const FieldElementZModPZ *B=(const FieldElementZModPZ*)&b;
      assert(A);
      assert(B);
      value=modP(value+modP(A->value*B->value,getF().getP()),getF().getP());
    }
  FieldElementZModPZ *one()const
  {
    return new FieldElementZModPZ(getF(),1);
  }
  bool isZero() const
  {
    return value==0;
  }

  FieldElementZModPZ *sum(const FieldElementImplementation &b)const
    {
      const FieldElementZModPZ *B=(const FieldElementZModPZ*)&b;
      FieldElementZModPZ *r= new FieldElementZModPZ(getF());
      r->value=modP(value+B->value,getF().getP());
      return r;
    }
  FieldElementZModPZ *difference(const FieldElementImplementation &b)const
    {
      const FieldElementZModPZ *B=(const FieldElementZModPZ*)&b;
      FieldElementZModPZ *r= new FieldElementZModPZ(getF());
      r->value=modP(value-B->value,getF().getP());
      return r;
    }
  FieldElementZModPZ *negation()const
    {
      FieldElementZModPZ *r= new FieldElementZModPZ(getF());
      r->value=modP(-value,getF().getP());

      return r;
    }
  FieldElementZModPZ *inverse()const
  {
    FieldElementZModPZ *r= new FieldElementZModPZ(getF());

    if(isZero())
      {
	AsciiPrinter P_(Stderr);
	P_.printString("Error inverting FieldElement: ");
	//	P.printFieldElement(*this);
	P_.printString("\n");
	assert(0);
      }

    //r->value=value; // Very simple when P=2

    int R,S;
    int d=gcd(value,getF().getP(),R,S);
    R=modP(R,getF().getP());
    assert(d==1);
    r->value=R;


    if(!(modP(R*value,getF().getP())==1))
      {
	fprintf(Stderr,"%i, %i\n",R,value);
      }
    assert(modP(R*value,getF().getP())==1);

    return r;
  }

  static string intToString(int v)
  {
    char s[20];
    sprintf(s,"%i",v);
    return string(s);
  }

  std::string toString(bool writeIfOne=true, bool alwaysWriteSign=false, bool latexMode=false) const
  {
    bool isOne=(value==1);
    bool isMinusOne=false;//...

    isMinusOne=(value==getF().getP()-1) && (value!=1);

    if(!writeIfOne && isOne)
      {
	if(alwaysWriteSign)return std::string("+");
	return std::string("");
      }
    if(!writeIfOne && isMinusOne)
      return std::string("-");

    if(alwaysWriteSign )
      return std::string("+")+intToString(value);
    return intToString(value);
  }

  FieldElementZModPZ *copy()const
  {
    FieldElementZModPZ *r= new FieldElementZModPZ(getF());

    r->value=value;
    return r;
  }
  int integerRepresentative()const
  {
    return value;
  }
};

FieldElementImplementation *FieldZModPZImplementation::zHomomorphismImplementation(int n)
{
  FieldElementImplementation *ret=new FieldElementZModPZ(*this,n);
  //  ret->refCount++; //!!!!!!!!!!!
  return ret;
}


FieldElement FieldZModPZImplementation::zHomomorphism(int n)
{
  return FieldElement(new FieldElementZModPZ(*this,n));
}

//static FieldZModPZ theField(2);
//static FieldZModPZ theField(32003);

/*Field *theZMod2ZField()
{
  return &theField;
  }*/
