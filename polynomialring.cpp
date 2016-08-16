#include "polynomialring.h"
#include <sstream>
#include <string>
#include <sstream>
#include "printer.h"
#include "log.h"
#include "polynomial.h"

#define plog1 if(0)

static int log10(int x)
{
  if(x>=10000)return 5;
  if(x>=1000)return 4;
  if(x>=100)return 3;
  if(x>=10)return 2;
  if(x>=1)return 1;
  return 0;
}

PolynomialRing::PolynomialRing(Field const &f, int numberOfVariables)
{
  vector<string> names;
  for(int i=0;i<numberOfVariables;i++)
    {
      stringstream s;
      s<<'x';
      for(int j=0;j<log10(numberOfVariables-1)-log10(i);j++)s<<'0';
      if(i!=0)s<<i;
      names.push_back(s.str());
    }
  implementingObject= new PolynomialRingImplementation(f,numberOfVariables,names);
  implementingObject->refCount++;
  plog1 fprintf(Stderr,"Constructing PolynomialRing\n");
}

PolynomialRing::PolynomialRing(Field const &f, vector<string> const &variables)
{
  implementingObject= new PolynomialRingImplementation(f,variables.size(),variables);
  implementingObject->refCount++;
  plog1 fprintf(Stderr,"Constructing PolynomialRing\n");
}

PolynomialRing PolynomialRing::withVariablesAppended(string variableNames)const
{
  vector<string> names2=implementingObject->variableNames;
  names2.push_back(variableNames);///!!!!!!!
  return PolynomialRing(implementingObject->theField,names2);
}


int PolynomialRing::variableIndex(string const &name)const
{
  for(int i=0;i<implementingObject->variableNames.size();i++)
    {
      if(implementingObject->variableNames[i]==name)return i;
    }
  return -1;
}


string const &PolynomialRing::getVariableName(int i)const
{
  assert(implementingObject);
  assert(i>=0);
  assert(i<implementingObject->variableNames.size());
  return implementingObject->variableNames[i];
}


vector<string> PolynomialRing::getVariableNames()const
{
  return implementingObject->variableNames;
}


string PolynomialRing::toStringVariableNames()const
{
	stringstream s;

	for(int i=0;i<getNumberOfVariables();i++)
	{
		if(i!=0)s<<",";
		s<<getVariableName(i);
	}

	return s.str();
}


PolynomialRing::~PolynomialRing()
{
  assert(implementingObject);
  implementingObject->refCount--;
  assert(implementingObject->refCount>=0);
  if(implementingObject->refCount==0)
    {
      plog1 fprintf(Stderr,"Deleting implementing object\n");
      delete implementingObject;
    }
  implementingObject=0;
  plog1 fprintf(Stderr,"Destructing PolynomialRing\n");
}


PolynomialRing::PolynomialRing(PolynomialRing const &a)
  :implementingObject(a.implementingObject)
{
  implementingObject->refCount++;
  plog1 fprintf(Stderr,"Copying PolynomialRing\n");
}


PolynomialRing& PolynomialRing::operator=(const PolynomialRing& a)
{
  plog1 fprintf(Stderr,"Assigning Field\n");
  if(this==&a)
    {
      plog1 fprintf(Stderr,"---selfassigning\n");
      return *this;
    }

  if(implementingObject&& 0==(--(implementingObject->refCount)))delete implementingObject;
  //assert(a.implementingObject);
  if(a.implementingObject)
    {
      implementingObject=a.implementingObject;
      implementingObject->refCount++;
    }
  else
    implementingObject=0;
  return *this;
}


vector<string> matrixVariableNames(string base, int height, int width)
{
  vector<string> ret;
  for(int i=0;i<height;i++)
    for(int j=0;j<width;j++)
      {
	stringstream s;
	// s<<base<<"["<<i<<"]["<<j<<"]";
	s<<base<<i<<j;
	ret.push_back(s.str());
      }
  return ret;
}


vector<string> vectorVariableNames(string base, int n)
{
  vector<string> ret;
  for(int i=0;i<n;i++)
    {
      stringstream s;
      // s<<base<<"["<<i<<"]["<<j<<"]";
      s<<base<<i;
      ret.push_back(s.str());
    }
  return ret;
}


vector<string> subsetVariableNames(string base, int n, int choose, bool M2Convention)
{
  vector<string> ret;
  vector<int> I;
  if(M2Convention)
    {
      for(int i=0;i<n-choose;i++)I.push_back(0);
      for(int i=0;i<choose;i++)I.push_back(1);
    }
  else
    {
      for(int i=0;i<choose;i++)I.push_back(1);
      for(int i=0;i<n-choose;i++)I.push_back(0);
    }
  do
  {
    stringstream s;
    s<<base;
    for(int i=0;i<n;i++)
      if(I[i])s<<i;
    ret.push_back(s.str());
  }
  while(M2Convention?next_permutation(I.begin(),I.end()):prev_permutation(I.begin(),I.end()));
  return ret;
}


int subsetToVariableIndex(set<int> const &s, int n, int choose, bool M2Convention)
{
  int ret=0;
  vector<int> I;
  if(M2Convention)
    {
      for(int i=0;i<n-choose;i++)I.push_back(0);
      for(int i=0;i<choose;i++)I.push_back(1);
    }
  else
    {
      for(int i=0;i<choose;i++)I.push_back(1);
      for(int i=0;i<n-choose;i++)I.push_back(0);
    }
  do
    {
      bool isRightSet=true;
      for(set<int>::const_iterator i=s.begin();i!=s.end();i++)
	{
	  if(I[*i]!=1)isRightSet=false;
	}
      if(isRightSet)break;
      ret++;
    }
  while(M2Convention?next_permutation(I.begin(),I.end()):prev_permutation(I.begin(),I.end()));

  return ret;
}


Polynomial PolynomialRing::one()const
{
  return Polynomial(Term(getField().zHomomorphism(1),Monomial(*this,IntegerVector(getNumberOfVariables()))));
}


Polynomial PolynomialRing::zero()const
{
  return Polynomial(*this);
}

Polynomial PolynomialRing::monomialFromExponentVector(IntegerVector const &v)const
{
  return Polynomial(Term(getField().zHomomorphism(1),Monomial(*this,v)));
}

Polynomial PolynomialRing::ithVariable(int i)const
{
	return Polynomial(Term(getField().zHomomorphism(1),Monomial(*this,IntegerVector::standardVector(getNumberOfVariables(),i))));
}

Polynomial PolynomialRing::polynomialFromField(FieldElement const &c)const
{
	return Polynomial(Term(c,Monomial(*this,IntegerVector(getNumberOfVariables()))));
}
