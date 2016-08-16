#include "triangulation.h"
#include <set>
#include <iostream>
#include "polyhedralcone.h"
#include "determinant.h"
#include "lp.h"
#include "linalg.h"
#include "log.h"

//---------------------------------------
// Printing routines
//---------------------------------------

void printIntList(list<int> const &v)
{
  FILE *f=Stderr;
  fprintf(f,"{");
  for(list<int>::const_iterator i=v.begin();i!=v.end();i++)
    {
      if(i!=v.begin())fprintf(f," ");
      fprintf(f,"%i",*i);
    }
  fprintf(f,"}\n");
}

static void printCone(Triangulation::Cone const &c)
{
  FILE *f=Stderr;
  fprintf(f,"{");
  for(list<int>::const_iterator i=c.begin();i!=c.end();i++)
    {
      if(i!=c.begin())fprintf(f," ");
      fprintf(f,"%i",*i);
    }
  fprintf(f,"} %i",c.changeSign);
}

void printIntListList(list<list<int> > const &l)
{
  for(list<list<int> >::const_iterator i=l.begin();i!=l.end();i++)
    printIntList(*i);
}

void printIntListList(list<Triangulation::Cone > const &l)
{
  cerr<<"{";
  for(list<Triangulation::Cone >::const_iterator i=l.begin();i!=l.end();i++)
    {
      if(i!=l.begin())cerr<<",\n";
      printCone(*i);
    }
  cerr<<"}\n";
}


//---------------------------------------
//
//---------------------------------------

static int signSwaps(list<int> const &c)
{
  int ret=1;
  for(list<int>::const_iterator i=c.begin();i!=c.end();i++)
    for(list<int>::const_iterator j=i;j!=c.end();j++)
      if(*j<*i)ret=-ret;
  return ret;
}

void Triangulation::coneSort(Triangulation::Cone &c)
{
  c.changeSign*=signSwaps(c);

  c.sort();
}

IntegerVectorList Triangulation::coneToVectorList(Cone2 const &c, IntegerMatrix const &rays)
{
  IntegerVectorList ret;

  for(Cone::const_iterator i=c.begin();i!=c.end();i++)
    ret.push_back(rays[*i]);

  return ret;
}

int Triangulation::coneDim(Cone2 const &c, IntegerMatrix const &rays)
{
  return rankOfMatrix(coneToVectorList(c,rays));
}

Triangulation::Cone Triangulation::firstSimplex(Cone const &c, IntegerMatrix const &rays)
{
  Cone ret;
  int d=0;

  for(Cone::const_iterator i=c.begin();i!=c.end();i++)
    {
      ret.push_back(*i);
      if(coneDim(ret,rays)!=ret.size())ret.pop_back();
    }
  return ret;
}

IntegerVectorList Triangulation::coneComplement(Cone c, IntegerMatrix const &rays)//returns generators of orth. complement.
{
  IntegerVectorList equations=coneToVectorList(c,rays);
  IntegerVectorList empty;
  return PolyhedralCone(empty,equations,rays.getWidth()).dualCone().getEquations();
}

int Triangulation::signVisible(int v, Cone const &c, IntegerVectorList const &complement, IntegerMatrix const &rays)
{
  IntegerVectorList l1=coneToVectorList(c,rays);
  l1.push_back(rays[v]);
  for(IntegerVectorList::const_iterator i=complement.begin();i!=complement.end();i++)
    {
      l1.push_back(*i);
    }
  return determinantSign(l1)*c.changeSign;
}

list<Triangulation::Cone> Triangulation::triangulateRek(int d, Cone2 const &c, IntegerMatrix const &rays, bool revlex, bool ignoreContainedRays)
{
  int revlexSignChange=(revlex?-1:1);
  list<Cone> ret;

  if(d==0)
    {
      ret.push_back(Cone());
      return ret;
    }

  //  log0 cerr<<"recursing"<<endl;

  ret=triangulateRek(d-1,c,rays,revlex);

  /* Check that rays do span a d-dimensional subspace,  and let *i denote the first entry of c such that the span is d */
  Cone l;
  Cone::const_iterator i;
  for(i=c.begin();i!=c.end();i++)
    {
      l.push_back(*i);
      if(coneDim(l,rays)==d)break;
    }
  assert(i!=c.end());

  /* Compute the orthogonal complement of the cone
   */
  IntegerVectorList complement;
  {
    Cone temp=*ret.begin();
    temp.push_back(*i);
    complement=coneComplement(temp,rays);
  }

  list<Cone>    boundary;

  /* Make two copies of the lower dimensional triangulation one for each orientation */
  for(list<Cone>::iterator i2=ret.begin();i2!=ret.end();i2++)
    {
      coneSort(*i2);
      boundary.push_back(*i2);
      i2->changeSign*=-1;
      boundary.push_back(*i2);
    }

  /* Build up the triangulation in new dimension */
  ret = list<Cone>();

  for(;i!=c.end();i++)
    {
      //      log0 cerr << "progress\n";
      {//We are done if we leave the subspace we are working in:
	bool done=false;
	for(IntegerVectorList::const_iterator j=complement.begin();j!=complement.end();j++)
	  if(dotLong(*j,rays[*i])!=0)
	    {
	      done=true;
	      break;
	    }
	if(done)break;
      }

      if(revlex)ret = list<Cone>();
      //      list<Cone> newBoundary;
      set<Cone> newBoundary;
      /* Run through old boundary */
      for(list<Cone>::iterator j=boundary.begin();j!=boundary.end();)
	{
	  /* For every visible triangle in the boundary */
	  if(signVisible(*i,*j,complement,rays)==1*revlexSignChange)
	    {
	      Cone b=*j;
	      /* remove the boundary triangle */
	      list<Cone>::iterator tempj=j;
	      j++;
	      if(!revlex)boundary.erase(tempj);

	      /* xor all other facets of the new higher dimensional simplex into the new boundary set */
	      for(Cone::const_iterator k=b.begin();k!=b.end();k++)
		{
		  Cone temp;
		  for(Cone::const_iterator l=b.begin();l!=b.end();l++)
		    if(l!=k)
		      temp.push_back(*l);
		    else
		      temp.push_back(*i);
		  {
		    temp.changeSign=b.changeSign;
		    if(revlex)temp.changeSign=-b.changeSign;
		    coneSort(temp);
		    bool found=false;
		    /*		    for(list<Cone>::iterator l=newBoundary.begin();l!=newBoundary.end();l++)
		      if(*l==temp)
			{
			  newBoundary.erase(l);
			  found=true;
			  break;
			}
		    if(!found)newBoundary.push_back(temp);
		    */
		    if(newBoundary.count(temp))newBoundary.erase(temp);else newBoundary.insert(temp);
		  }
		}
	      /* build higher dimensional triangle and add it to the triangulation */
	      b.push_back(*i);
	      b.changeSign=-b.changeSign;
	      coneSort(b);
	      ret.push_back(b);
	    }
	  else
	    {
	      /* remove the boundary triangle */
	      list<Cone>::iterator tempj=j;
	      j++;
	      if(revlex)boundary.erase(tempj);
	    }
	}
      //      boundary.splice(boundary.begin(),newBoundary);
      for(set<Cone>::const_iterator i=newBoundary.begin();i!=newBoundary.end();i++)boundary.push_back(*i);
      //log0 cerr<< "ret.size:"<<ret.size()<<"boudary.size"<<boundary.size()<<endl;
      //printIntListList(ret);
      //printIntListList(boundary);
    }
  return ret;
}

list<Triangulation::Cone> Triangulation::triangulate(Cone2 c, IntegerMatrix const &rays, bool revlex) //computes a lexicographic triangulation
{
  //  coneSort(c);
  c.sort();
  return triangulateRek(coneDim(c,rays),c,rays,revlex);
}

list<Triangulation::Cone> Triangulation::triangulate(IntegerMatrix const &rays, bool revlex)
{
  /* Reduce to subspace. Could this change orientation??? */
  FieldMatrix R=integerMatrixToFieldMatrix(rays,Q).transposed();
  R.reduce();
  R.removeZeroRows();
  IntegerMatrix rays2=fieldMatrixToIntegerMatrixPrimitive(R.transposed());

  Cone c;
  for(int i=0;i<rays2.getHeight();i++)c.push_back(i);
  return triangulate(c,rays2,revlex);
}


list<Triangulation::Cone> Triangulation::boundary(list<Cone> cones)
{
  set<Cone> ret;

  //  printIntListList(cones);
  for(list<Cone>::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      Cone c=*i;
      //      coneSort(c);
      int doSwap=-1;
      for(Cone::iterator j=c.begin();j!=c.end();j++)
	{
	  doSwap=-doSwap;
	  Cone c2;
	  for(Cone::iterator k=c.begin();k!=c.end();k++)
	    if(k!=j)c2.push_back(*k);
	  c2.changeSign=doSwap*c.changeSign;
	  coneSort(c2);
	  for(set<Cone>::iterator l=ret.begin();l!=ret.end();l++)
	    {
	      Cone a=*l;
	      Cone b=c2;
	      coneSort(a);
	      coneSort(b);

	      if((!(a<b)&&(!(b<a)))){ret.erase(l);goto erased;};
	    }
	  ret.insert(c2);
	erased:
	  ;
	}
    }
  list<Cone> ret2;
  for(set<Cone>::const_iterator i=ret.begin();i!=ret.end();i++)ret2.push_back(*i);
  return ret2;
}


IntegerVectorList Triangulation::normals(IntegerMatrix &rays)
{
  list<Cone> b=boundary(triangulate(rays));



  FieldMatrix raySpan=integerMatrixToFieldMatrix(rays,Q);
  FieldMatrix rayKernel=raySpan.reduceAndComputeKernel();

  rayKernel.reduce();
  //  fprintf(Stderr,"RESULT\n");

  //  printIntListList(b);

  // AsciiPrinter P(Stderr);
  // fprintf(Stderr,"Reducer:\n");
  //  rayKernel.print(P);
  //P.printVectorList(rays.getRows());
  //raySpan.printMatrix(P);
  set<IntegerVector> ret2;

  for(list<Cone>::const_iterator i=b.begin();i!=b.end();i++)
    {
      //      fprintf(Stderr,"HEJ!!\n");

      //AsciiPrinter P(Stderr);
      IntegerVectorList temp;
      for(Cone::const_iterator j=i->begin();j!=i->end();j++)
	temp.push_back(rays[*j]);
      //printIntList(*i);
      //fprintf(Stderr,"test\n");
      FieldMatrix raySpan2=integerMatrixToFieldMatrix(rowsToIntegerMatrix(temp),Q);

      FieldMatrix raySpan3=combineOnTop(raySpan2,rayKernel);

      //      raySpan3.printMatrix(P);
      //fprintf(Stderr,"test2\n");
      FieldMatrix rayKernel2=raySpan3.reduceAndComputeKernel();
      //fprintf(Stderr,"RaySpan2:\n");
      //raySpan2.printMatrix(P);
      //rayKernel2.printMatrix(P);
      //fprintf(Stderr,"test3\n");

      assert(rayKernel2.getHeight()==1);
      FieldVector v=rayKernel2[0];


      int swaps=0;
      //fprintf(Stderr,"test4\n");
      for(int j=0;j<rays.getHeight();j++)
	{
	  //fprintf(Stderr,"test5\n");
	  FieldElement d=dot(v,integerVectorToFieldVector(rays[j],Q));
	  if(d.sign()<0)
	    {
	      v=(Q.zHomomorphism(-1))*v;
	      swaps++;
	    }
	}
      assert(swaps<2);


      //fprintf(Stderr,"HEJ!!\n");
      //rayKernel.printMatrix(P);
      //v.print(P);
      //      rayKernel.canonicalize(v).print(P);
      //AsciiPrinter(Stderr).printVector(rayKernel.canonicalize(v).primitive());
      //      ret2.insert(rayKernel.canonicalize(v).primitive());
      ret2.insert(rayKernel.canonicalize(v).primitive());

    }
  IntegerVectorList ret;
  for(set<IntegerVector>::const_iterator i=ret2.begin();i!=ret2.end();i++)
    ret.push_back(*i);

  return ret;
}


vector<list<int> > Triangulation::removeOrientation(list<Cone> const &triangulation)
{
  vector<list<int> > ret(triangulation.size());

  int I=0;
  for(list<Cone>::const_iterator i=triangulation.begin();i!=triangulation.end();i++,I++)
    ret[I]=*i;

  return ret;
}


vector<list<int> > Triangulation::removeOrientation(vector<Cone> const &triangulation)
{
  vector<list<int> > ret(triangulation.size());

  int I=0;
  for(vector<Cone>::const_iterator i=triangulation.begin();i!=triangulation.end();i++,I++)
    ret[I]=*i;

  return ret;
}
