#include "minors.h"
#include "printer.h"
#include <algorithm>

using namespace std;

static int lookup(vector<int> const &v, int i)
{
  for(int j=0;j<v.size();j++)
    {
      if(v[j])
	{
	  if(i==0)
	    return j;
	  i--;
	}
    }
  return 0;
}

static int sign(vector<int> const &v)
{
  int ret=1;
  for(int i=0;i<v.size();i++)
    for(int j=i+1;j<v.size();j++)
      if(v[i]>v[j])ret=-ret;

  return ret;
}

PolynomialSet minors(PolynomialRing const &R, int r, int d, int n, bool withNames, bool M2Convention)
{
  int namesOffset=d*n;
  int entriesOffset=0;
  PolynomialSet ret(R);
  vector<int> I;
  if(M2Convention)
    {
      for(int i=0;i<d-r;i++)I.push_back(0);
      for(int i=0;i<r;i++)I.push_back(1);
    }
  else
    {
      for(int i=0;i<r;i++)I.push_back(1);
      for(int i=0;i<d-r;i++)I.push_back(0);
    }

  do
    {
      vector<int> J;
      if(M2Convention)
	{
	  for(int i=0;i<n-r;i++)J.push_back(0);
	  for(int i=0;i<r;i++)J.push_back(1);
	}
      else
	{
	  for(int i=0;i<r;i++)J.push_back(1);
	  for(int i=0;i<n-r;i++)J.push_back(0);
	}
      do
	{
	  Polynomial p(R);
	  if(withNames)
	    p+=Term(R.getField().zHomomorphism(-1),Monomial(R,IntegerVector::standardVector(R.getNumberOfVariables(),namesOffset++)));
	    
	  vector<int> K;
	  for(int i=0;i<r;i++)K.push_back(i);
	  do
	    {
	      IntegerVector v(R.getNumberOfVariables());
	      vector<int>::const_iterator s=K.begin();
	      for(int i=0;i<r;i++,s++)
		{
		  int y=lookup(I,i);
		  int x=lookup(J,*s);
		  int index=y*n+x;
		  v+=IntegerVector::standardVector(R.getNumberOfVariables(),index+entriesOffset);
		}
	      p+=Term(R.getField().zHomomorphism(sign(K)),Monomial(R,v));
	    }
	  while(next_permutation(K.begin(),K.end()));

	  ret.push_back(p);
	}
      while(M2Convention?next_permutation(J.begin(),J.end()):prev_permutation(J.begin(),J.end()));
    }
  while(M2Convention?next_permutation(I.begin(),I.end()):prev_permutation(I.begin(),I.end()));
  return ret;
}
