#include "graph.h"

#include <sstream>
#include "matrix.h"

int Graph::diameter(int *timesAttained)const
{
  IntegerMatrix A(n,n);

  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      A[i][j]=n;
  for(set<Edge>::const_iterator i=edges.begin();i!=edges.end();i++)
    {
      A[i->a][i->b]=1;
      if(!isDirected)A[i->b][i->a]=1;
    }
  for(int i=0;i<n;i++)A[i][i]=0;

  bool wasChange=true;
  while(wasChange)
    {
      wasChange=false;
      for(int i=0;i<n;i++)
	for(int j=0;j<n;j++)
	  {
	    int best=A[i][j];
	    for(int k=0;k<n;k++)
	      {
		int s=A[i][k]+A[k][j];
		if(s<best)
		  {
		    best=s;
		    wasChange=true;
		  }
	      }
	    A[i][j]=best;
	  }
    }

  int timesAttained2=0;
  int ret=0;
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      {
	if(A[i][j]>ret)
	  {
	    ret=A[i][j];
	    timesAttained2=0;
	  }
	if(A[i][j]==ret)
	  timesAttained2++;
      }

  if(timesAttained)*timesAttained=timesAttained2;

  return ret;
}


string Graph::toString()const
{
  stringstream ret;

  ret<<"("<<n<<",{"<<endl;
  for(set<Edge>::const_iterator i=edges.begin();i!=edges.end();i++)
    {
      if(i!=edges.begin())ret<<","<<endl;
      ret<<"("<<i->a<<","<<i->b<<")"<<endl;
    }
  ret<<"}"<<endl;

  return ret.str();
}
