#include "scarf.h"

#include <vector>
#include <set>
#include <algorithm>

#include "latticeideal.h"
#include "subspace.h"
#include "printer.h"
#include "polyhedralcone.h"

IntegerVectorList neighbours(IntegerMatrix const &A)
{
  IntegerVectorList testSet=latticeIdealRevLex(A);
  

  IntegerMatrix T=A.transposed();
  IntegerMatrix B=IntegerMatrix::identity(T.getWidth());
  B.append(T);
  B=B.transposed();
  //  AsciiPrinter(Stderr).printVectorList(B.getRows());
  
  //  AsciiPrinter Q(Stderr);
  IntegerVectorList temp;
  PolyhedralCone c(temp,B.getRows(),B.getWidth());
  //  c.canonicalize();
  c=c.dualCone();
  Subspace U(c.getEquations(),c.ambientDimension());
  
  IntegerVectorList neighbours;
  for(IntegerVectorList::const_iterator i=testSet.begin();i!=testSet.end();i++)
    {
      IntegerVector v=*i;
      IntegerVector w=v;
      w.grow(c.ambientDimension());
      IntegerVector u=U.canonicalizeVector(w).subvector(A.getHeight(),w.size());
      neighbours.push_back(u);
      //  Q.printVector(v);
      //  Q.printVector(A.vectormultiply(u));
      assert(A.vectormultiply(u)==v); //the method depends on the implementation of Subspace
    }

  return neighbours;
}


bool satisfiesA1(IntegerMatrix const &A)
{
  IntegerMatrix T=A.transposed();
  IntegerVectorList temp;
  PolyhedralCone c(temp,T.getRows(),T.getWidth());


  return intersection(PolyhedralCone::positiveOrthant(T.getWidth()),c).containsPositiveVector();

  //PolyhedralCone d=intersection(PolyhedralCone::positiveOrthant(T.getWidth()),c);
  //AsciiPrinter Q(Stderr);
  //d.print(&Q);
  //return intersection(PolyhedralCone::positiveOrthant(T.getWidth()),c).containsPositiveVector();
}


bool satisfiesA2(IntegerMatrix const &A)
{
  int n=A.getWidth();
  int m=A.getHeight();
  assert(m>=n);
  IntegerVector v(m);
  for(int i=0;i<n;i++)v[i]=-1;

  do
    {
       IntegerVectorList l;
       IntegerVectorList k=A.transposed().getRows();
       
       for(int i=0;i<m;i++)
	 if(v[i])
	   l.push_back(IntegerVector::standardVector(m,i));
	 else
	   k.push_back(IntegerVector::standardVector(m,i));

       PolyhedralCone c(l,k,m);

       if(!c.isZero())return false;
    }
  while(v.nextPermutation());

  return true;
}


bool satisfiesA3(IntegerMatrix const &A, IntegerVectorList const *N)
{
  IntegerVectorList N2;
  if(!N)
    {
      N2=neighbours(A);
      N=&N2;
    }
  for(IntegerVectorList::const_iterator i=N->begin();i!=N->end();i++)
    for(int j=0;j<A.getHeight();j++)
      {
	if(dot(A[j],*i)==0)return false;
      }
  return true;
}


bool satisfiesA3i(IntegerMatrix const &A, int j, IntegerVectorList const *N)
{
  assert(j>=0 && j<A.getHeight());
  IntegerVectorList N2;
  if(!N)
    {
      N2=neighbours(A);
      N=&N2;
    }
  for(IntegerVectorList::const_iterator i=N->begin();i!=N->end();i++)if(dot(A[j],*i)==0)return false;

  return true;
}


IntegerVectorList orientedNeighbours(IntegerVectorList const& N, IntegerVector const &v)
{
  IntegerVectorList ret;
  for(IntegerVectorList::const_iterator i=N.begin();i!=N.end();i++)
    {
      if(dot(v,*i)>0)
	ret.push_back(*i);
      else
	ret.push_back(-*i);	
    }
  return ret;
}


static IntegerVector shift(IntegerMatrix const &N, IntegerVector simplex, IntegerVector const &offset)
{

  {
  for(int i=0;i<simplex.size();i++)
    {
      IntegerVector newNeighbour=N[simplex[i]]+offset;
      int newIndex=-1;
      for(int j=0;j<N.getHeight();j++)
	if(N[j]==newNeighbour)newIndex=j;
      assert(newIndex!=-1);
      simplex[i]=newIndex;
    }
  }
  return simplex;
}


IntegerVector kFlip(IntegerMatrix const &A, IntegerMatrix const &N, IntegerVector simplex, int vertex)
{
  int vertex2=-1;
  int best=-100000;
  for(int i=0;i<simplex.size();i++)
    {
      if(i!=vertex)
	{
	  int d=dot(N[simplex[i]],A[vertex]);
	  if(d>best)
	    {
	      best=d;
	      vertex2=i;
	    }
	}
    }
  assert(vertex2!=-1);
  simplex[vertex]=simplex[vertex2];

  simplex=shift(N,simplex,-N[simplex[0]]);

  best=10000;
  int newRow=-1;
  for(int i=0;i<N.getHeight();i++)
    {
      bool inside=true;
      for(int j=0;j<simplex.size();j++)
	{
	  if(j!=vertex2)
	    if(dot(N[i]-N[simplex[j]],A[j])>=0)inside=false;
	}
      if(inside)
	{
	  int d=dot(N[i],A[vertex2]);
	  if(d<best)
	    {
	      best=d;
	      newRow=i;
	    }
	}
    }
  assert(newRow!=-1);

  simplex[vertex2]=newRow;

  simplex=shift(N,simplex,-N[simplex[0]]);
  return simplex;
}

void traverseScarfComplex(IntegerMatrix const &A, IntegerMatrix const &N, IntegerVector simplex)
{
  set<IntegerVector> simplices;
  IntegerVectorList active;

  active.push_back(simplex);

  while(!active.empty())
    {
      IntegerVector s=active.front();

      if(simplices.count(s)==0)
	{
	  fprintf(Stderr,"processing:");
	  AsciiPrinter(Stderr).printVector(s);
	  simplices.insert(s);
	  fprintf(Stderr,"\n");
	    
	  for(int i=0;i<simplex.size();i++)
	    {
	      IntegerVector s2=kFlip(A,N,s,i);
	      
	      active.push_back(s2);
	    }
	}
      active.pop_front();
    }
  
}


int growPolytope(IntegerMatrix const &A, IntegerMatrix const &N, IntegerVector simplex, int inequality, int n)
{
  int best=10000;
  int newRow=-1;
  for(int i=0;i<N.getHeight();i++)
    {
      bool inside=true;
      for(int j=0;j<n;j++)
	{
	  if(j!=inequality)
	    if(dot(N[i]-N[simplex[j]],A[j])>=0)inside=false;
	}
      /*      fprintf(Stderr,"Testing ");
      AsciiPrinter(Stderr).printVector(N[i]);
      fprintf(Stderr," inside:%i\n",inside);
      */
      for(int j=0;j<n;j++)
	if(dot(N[i]-N[simplex[j]],A[inequality])<0)inside=false;
      //      fprintf(Stderr," inside:%i\n",inside);

      if(inside)
	{
	  int d=dot(N[i],A[inequality]);
	  if(d<best)
	    {
	      best=d;
	      newRow=i;
	    }
	}
    }
  /*  fprintf(Stderr,"Best ");
  AsciiPrinter(Stderr).printVector(N[newRow]);
  fprintf(Stderr,"\n");
  */
  assert(newRow!=-1);
  return newRow;
}

IntegerVector computeMaximalScarfSimplex(IntegerMatrix const &A, IntegerMatrix const &N)
{
  IntegerVector s(A.getHeight());

  for(int i=1;i<s.size();i++)
    s[i]=growPolytope(A,N,s,i,i);

  return s;
}
