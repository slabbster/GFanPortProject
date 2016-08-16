#include "polynomial.h"
#include "field_rationals.h"
#include "printer.h"
#include "log.h"

#include <sstream>
#include <iostream>

static int rIndex(int i, int j, int N, bool withMasses)
{
  if(j<i){i^=j;j^=i;i^=j;}
  int r=withMasses*N;
  for(int I=0;I<N;I++)
    for(int J=I+1;J<N;J++)
      {
	if((I==i)&&(J==j))return r;
	r++;
      }
  assert(0);
}

static int sIndex(int i, int j, int N, bool withMasses)
{
  return rIndex(i,j,N,withMasses)+N*(N-1)/2;
}

static Polynomial S(PolynomialRing const &r, int N, bool withMasses, int i, int j, bool withSVariables)
{
  Polynomial ret(r);
  int n=r.getNumberOfVariables();

  if(withSVariables)
    {
      ret+=Term(r.getField().zHomomorphism(1),Monomial(r,IntegerVector::standardVector(n,sIndex(i,j,N,withMasses))));
    }
  else
    {
      ret+=Term(r.getField().zHomomorphism(-1),Monomial(r,IntegerVector(n)));
      ret+=Term(r.getField().zHomomorphism(1),Monomial(r,-3*IntegerVector::standardVector(n,rIndex(i,j,N,withMasses))));
    }

  return ret;
}


static Polynomial rPolynomial(PolynomialRing const &r, int N, bool withMasses, int i, int j, int exponent=2)
{
  if(i==j)return Polynomial(r);
  return Term(r.getField().zHomomorphism(1),Monomial(r,exponent*IntegerVector::standardVector(r.getNumberOfVariables(),rIndex(i,j,N,withMasses))));
}


static vector<string> variableNames(int N, bool withMasses, bool withSVariables)
{
  vector<string> ret;

  if(withMasses)
    {
      for(int i=0;i<N;i++)
	{
	  stringstream s;
	  s<<"m"<<i+1;
	  ret.push_back(s.str());
	}
    }

  for(int i=0;i<N;i++)
    for(int j=i+1;j<N;j++)
      {
	stringstream s;
	s<<"r"<<i+1<<j+1;
	ret.push_back(s.str());
      }

  if(withSVariables)
    {
      for(int i=0;i<N;i++)
	for(int j=i+1;j<N;j++)
	  {
	    stringstream s;
	    s<<"S"<<i+1<<j+1;
	    ret.push_back(s.str());
	  }
    }

  return ret;
}


Polynomial AlbouyChencinerEquation(PolynomialRing const &r, int N, bool withMasses, int i, int j, bool symmetric, bool withSVariables, bool saturate)
{
  int n=r.getNumberOfVariables();
  Polynomial ret(r);
  for(int k=0;k<N;k++)
    {
      Term massTerm=Term(r.getField().zHomomorphism(1),Monomial(r,withMasses*IntegerVector::standardVector(n,k)));
      if(i!=k)
	ret+=massTerm*
	  (
	   S(r,N,withMasses,i,k,withSVariables)*(
		       rPolynomial(r,N,withMasses,j,k)
		       -rPolynomial(r,N,withMasses,i,k)
		       -rPolynomial(r,N,withMasses,i,j)
		       )
	   );
      if(symmetric)
	if(j!=k)
	  ret+=massTerm*
	    (
	     S(r,N,withMasses,j,k,withSVariables)*(
			 rPolynomial(r,N,withMasses,i,k)
			 -rPolynomial(r,N,withMasses,j,k)
			 -rPolynomial(r,N,withMasses,i,j)
			 )
	     );
    }
  if(saturate)ret.saturate();
  return ret;
}


PolynomialSet AlbouyChencinerEquations(int N, bool withMasses, bool symmetric, bool withSVariables, bool saturate)
{
  PolynomialRing r(Q,variableNames(N,withMasses,withSVariables));

  PolynomialSet ret(r);
  for(int i=0;i<N;i++)
    {
      if(!symmetric)
	for(int j=0;j<i;j++)
	  ret.push_back(AlbouyChencinerEquation(r,N,withMasses,i,j,symmetric,withSVariables,saturate));

      for(int j=i+1;j<N;j++)
	ret.push_back(AlbouyChencinerEquation(r,N,withMasses,i,j,symmetric,withSVariables,saturate));
    }
  return ret;
}


PolynomialSet DziobekEquations(PolynomialRing const &r, int N, bool withMasses, bool withSVariables, bool saturate)
{
  PolynomialSet ret(r);

  list<int> a;
  for(int i=4;i<N;i++)
    a.push_back(0);
  a.push_back(1);
  a.push_back(1);
  a.push_back(1);
  a.push_back(1);

  do
    {
      list<int>::const_iterator i=a.begin();
      int I=0;
      while((*i)==0){i++;I++;}
      int first=I;
      do{i++;I++;}while((*i)==0);
      int second=I;
      do{i++;I++;}while((*i)==0);
      int third=I;
      do{i++;I++;}while((*i)==0);
      int fourth=I;

      //      cerr<<first<<second<<third<<fourth;
      Polynomial f1=S(r,N,withMasses,first,second,withSVariables)*S(r,N,withMasses,third,fourth,withSVariables)-S(r,N,withMasses,first,third,withSVariables)*S(r,N,withMasses,second,fourth,withSVariables);
      if(saturate)f1.saturate();
      Polynomial f2=S(r,N,withMasses,first,third,withSVariables)*S(r,N,withMasses,second,fourth,withSVariables)-S(r,N,withMasses,first,fourth,withSVariables)*S(r,N,withMasses,second,third,withSVariables);
      if(saturate)f2.saturate();
      Polynomial f3=S(r,N,withMasses,first,second,withSVariables)*S(r,N,withMasses,third,fourth,withSVariables)-S(r,N,withMasses,first,fourth,withSVariables)*S(r,N,withMasses,second,third,withSVariables);
      if(saturate)f3.saturate();
      ret.push_back(f1);
      ret.push_back(f2);
      ret.push_back(f3);
    }
  while(next_permutation(a.begin(),a.end()));

  return ret;
}


static Polynomial mlookup(PolynomialRing const &r, int N, bool withMasses,int i,int j)
{
  if(i==j)return r.zero();
  if(i==0 || j==0)return r.one();
  return Polynomial(Term(r.getField().zHomomorphism(1),Monomial(r,2*IntegerVector::standardVector(r.getNumberOfVariables(),rIndex(i-1,j-1,N,withMasses)))));
}


PolynomialSet nbodyDeterminants(PolynomialRing const &r, int N, bool withMasses, int determinantSize)
{
  vector<int> l;

  for(int i=0;i<N+1-determinantSize;i++)
    l.push_back(1);
  for(int i=0;i<determinantSize-1;i++)
    l.push_back(0);

  PolynomialSet ret(r);

  if(determinantSize==N+2)return ret;

  do
    {
      vector<int> indexList;
      indexList.push_back(0);
      for(int i=0;i<l.size();i++)if(l[i]==0)indexList.push_back(i+1);

      log1
      {
    	  for(int A=0;A<determinantSize;A++)
    	  {
    		  for(int B=0;B<determinantSize;B++)
    			  AsciiPrinter(Stderr)<<mlookup(r,N,withMasses,indexList[A],indexList[B])<<";";
				  cerr<<endl;
    	  }
      }

      vector<int> perm;
      for(int i=0;i<indexList.size();i++)
	perm.push_back(i);

      Polynomial p(r);

      do
	{
	  Polynomial prod=r.one();
	  for(int j=0;j<perm.size();j++)
	    {
	      prod*=mlookup(r,N,withMasses,indexList[j],indexList[perm[j]]);
	    }
	  int s=1;
	  for(int x=0;x<perm.size();x++)
	    for(int y=0;y<x;y++)
	      if(perm[y]>perm[x])s*=-1;
	  if(s==1)
	    p+=prod;
	  else
	    p-=prod;
	}
      while(next_permutation(perm.begin(),perm.end()));
      ret.push_back(p);
    }
  while(prev_permutation(l.begin(),l.end()));

  return ret;
}


Polynomial massEquation(PolynomialRing const &r, int N, bool withMasses, bool saturate)
{
  Polynomial ret(r);
  for(int i=0;i<N;i++)
    for(int j=0;j<i;j++)
      {
	Polynomial mm=r.one();
	if(withMasses)mm=Polynomial(Term(r.getField().zHomomorphism(1),Monomial(r,IntegerVector::standardVector(r.getNumberOfVariables(),i)+IntegerVector::standardVector(r.getNumberOfVariables(),j))));
	ret+=(rPolynomial(r,N,withMasses,i,j,2)-rPolynomial(r,N,withMasses,i,j,-1))*mm;
      }
  if(saturate)ret.saturate();
  return ret;
}


Polynomial SEquation(PolynomialRing const &r, int N, bool withMasses, int i, int j, bool withSVariables, bool saturate=true)
{
  Polynomial ret(r);
  int n=r.getNumberOfVariables();

  ret-=Term(r.getField().zHomomorphism(1),Monomial(r,IntegerVector::standardVector(n,sIndex(i,j,N,withMasses))));

  ret+=Term(r.getField().zHomomorphism(-1),Monomial(r,IntegerVector(n)));
  ret+=Term(r.getField().zHomomorphism(1),Monomial(r,-3*IntegerVector::standardVector(n,rIndex(i,j,N,withMasses))));

  if(saturate)ret.saturate();

  return ret;
}


PolynomialSet SEquations(PolynomialRing const &r, int N, bool withMasses, bool saturate)
{
  PolynomialSet ret(r);
  for(int i=0;i<N;i++)
    for(int j=0;j<i;j++)
      {
	ret.push_back(SEquation(r,N,withMasses,i,j,true,saturate));
      }
  return ret;
}
