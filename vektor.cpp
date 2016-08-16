#include "vektor.h"
#include <stdio.h>
#include "printer.h"

void outOfRange(int i, int n)
{
  fprintf(Stderr,"Index out of range: index=%i vector lenght=%i\n",i,n);
  assert(0);
}

IntegerVectorList transposeIntegerVectorList(IntegerVectorList const &l)
{
  int n=l.size();
  assert(n);
  int m=l.begin()->size();

  IntegerVectorList r;

  for(int i=0;i<m;i++)
    {
      IntegerVector column(n);

      int j=0;
      for(IntegerVectorList::const_iterator J=l.begin();J!=l.end();J++)
	{
	  column[j]=(*J)[i];
	  j++;
	}
      r.push_back(column);
    }
  return r;
}


IntegerVectorList multiplyIntegerVectorList(IntegerVectorList const &A, IntegerVectorList const &B)
{
  int s1=A.size();
  assert(s1);
  int r=A.begin()->size();
  int t=B.size();
  assert(t);
  int s2=B.begin()->size();

  assert(s1==s2);

  IntegerVectorList C;
  for(IntegerVectorList::const_iterator I=B.begin();I!=B.end();I++)
  {
    IntegerVector column(r);

    int j=0;
    for(IntegerVectorList::const_iterator J=A.begin();J!=A.end();J++)
      {
	column+=(*I)[j]*(*J);
	j++;
      }
    C.push_back(column);
  }
  return C;
}


/* Computes the greatest common divisor of two integers.
   The result is always positive.
   Asserts if both integers are zero.
 */
int gcdGFAN(int r, int s)
{
  if(r<0)r=-r;
  if(s<0)s=-s;

  if(s<r)
    {
      int t;
      t=r;
      r=s;
      s=t;
    }

  while(r!=0)
    {
      int t=s%r;
      s=r;
      r=t;
    }
  assert(s!=0);

  return s;
}


/* Returns positive gcd of elements in the vector.
   The result is always positive.
   Asserts if the vector is the zero-vector.
 */
int gcdOfVector(IntegerVector const &v)
{
  int ret=0;
  for(int i=0;i<v.size();i++)if(ret=v[i])break;
  if(ret<0)ret=-ret;
  assert(ret!=0);
  for(int i=0;i<v.size();i++)ret=gcdGFAN(ret,v[i]);

  return ret;
}

int gcdOfVectorCandidate(IntegerVector const &v, int ret)
{
  for(int i=0;i<v.size();i++)
	  {
	  ret=gcdGFAN(ret,v[i]);
	  if(ret==1)break;
	  }

  return ret;
}

/* Uses gcd to put (the direction given by) the vector in a unique form by scaling.
   Only positive scaling is applied.
   If the input is the zero-vector so is the output.
   Otherwise, the common divisor of the output is 1.
*/
void normalizedLowLevel(IntegerVector const &v, IntegerVector &dest)
{
	int a;
	int n=v.size();
	assert(n==dest.size());
	int smallest=0x7fffffff;
	int zeroTest=0;
	for(int i=0;i<n;i++)
	{
		int vv=v.UNCHECKEDACCESS(i);
		if(vv>0)
		{
			if(vv<smallest)smallest=vv;
		}
		else if(vv<0)
		{
			if(-vv<smallest)smallest=-vv;
		}
		zeroTest|=vv;
	}
	if(smallest==1)goto returnv;

  if(zeroTest==0)goto returnv;
  a=gcdOfVectorCandidate(v,smallest);
  for(int i=0;i<n;i++)dest[i]=v.UNCHECKEDACCESS(i)/a;

  return;
  returnv:
  for(int i=0;i<n;i++)dest[i]=v.UNCHECKEDACCESS(i);
}

IntegerVector normalized(IntegerVector const &v)
{
	int n=v.size();
	int smallest=0x7fffffff;
	for(int i=0;i<n;i++)
	{
		if(v[i]>0)
		{
			if(v[i]<smallest)smallest=v[i];
		}
		else if(v[i]<0)
		{
			if(-(v[i])<smallest)smallest=-(v[i]);
		}
	}
	if(smallest==1)return v;

  if(v.isZero())return v;
  int a=gcdOfVectorCandidate(v,smallest);
//  int a=gcdOfVector(v);
  IntegerVector V(n);
  for(int i=0;i<n;i++)V[i]=v[i]/a;

  return V;
}


IntegerVectorList subvectorsOfIntegerVectorList(IntegerVectorList const &l, list<int> const &chosen)
{
  IntegerVectorList ret;

  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    ret.push_back(i->subvector(chosen));

  return ret;
}


void removeDuplicates(IntegerVectorList &l)
{
  l.sort();
  l.unique();
}
