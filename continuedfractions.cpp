#include "continuedfractions.h"
#include "vektor.h"

#include <stdio.h>

void doubleToFraction(double f, int &numerator, int &denominator, int maksIter)
{
  bool changeSign=f<0;
  if(changeSign)f*=-1;

  if(f<0.000001 || maksIter<1)
    {
      numerator=0;
      denominator=1;
    }
  else
    {
      double r=1/f;
      int R=int(r);//truncates
      r-=R;
      int n2,d2;
      doubleToFraction(r,n2,d2,maksIter-1);
      /* 1/f-R=r=n2/d2; => 1/f=(n2+d2*R)/d2; => f=d2/(n2+d2*R); */
      numerator=d2;
      denominator=n2+d2*R;

      if(changeSign)
	numerator*=-1;
    }
}



void doubleVectorToFractions(vector<double> v, vector<int> &numerators, int &denominator)
{
  int n=v.size();
  numerators=vector<int>(n);
  denominator=1;
  if(n==0)return;
  for(int i=0;i<n;i++)
    {
      int num,den;
      doubleToFraction(v[i],num,den);
      if(den==0)goto error;
      denominator=(((signed long long)den)*denominator)/gcdGFAN(denominator,den);
    }
  for(int i=0;i<n;i++)
    {
      int num,den;
      doubleToFraction(v[i],num,den);
      if(den!=0)
    	  numerators[i]=num*(denominator/den);
      else
    	  numerators[i]=num;//if v[i] is large it can happen that we lift to 1/0. Then we just produce some random result.
    }
  return;
error:
	for(int i=0;i<n;i++)numerators[i]=0;denominator=1;
}
