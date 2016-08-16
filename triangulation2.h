#ifndef TRIANGULATION2_H_INCLUDED
#define TRIANGULATION2_H_INCLUDED

#include <set>
#include <iostream>
#include <algorithm>
using namespace std;
#include "polyhedralfan.h"
#include "matrix.h"
#include "field_rationals.h"
#include "graph.h"
#include "determinant.h"
#include "lp.h"
#include "wallideal.h"
#include "linalg.h"
#include "printer.h"
#include "log.h"
/**
     TODO: Merge this class with the Triangulation class.
   */

  class Triangulation2
  {
    //    IntegerVector circuit(IntegerMatrix const &A, IntegerVector const &basis, int j);

    /**
       A projective vector configuration. The columns are the vectors.
     */
    IntegerMatrix A;
    IntegerMatrix Atransposed;
    /**
       The dual of a subdivision of the cone over A is a polyhedron with facet normals among the columns.
       The subdivision being a triangulation is equivalent to the polyhedron being simple.
       Bases is either the sets of facet normals giving rise to a poiont of the polyhedron
       OR the list of triangles of the triangulation.
     */
  public:
    Triangulation2(IntegerMatrix const &A_):
      A(A_),
      Atransposed(A_.transposed())
      {
      }
    set<IntegerVector> bases;
    int getN()const
    {
      return A.getWidth();
    }
    int getD()const
    {
      return A.getHeight();
    }
    IntegerVector complement(IntegerVector const &v, int n)const
    {
      IntegerVector ret(0);
      int j=0;
      for(int i=0;i<n;i++)
	{
	  if(j>=v.size())
	    ret.push_back(i);
	  else
	    if(i==v[j])
	      j++;
	    else
	      {
		ret.push_back(i);
	      }
	}
      return ret;
    }
    /**
     * Test if the configuration
     */
    bool isHomogeneous()
    {
/*      FieldMatrix temp=integerMatrixToFieldMatrix(A.transposed(),Q);
      FieldMatrix s=temp.solver();
      FieldVector y=s.canonicalize(concatenation(FieldVector::allOnes(Q,temp.getHeight()),FieldVector(Q,temp.getWidth())));
      bool ret=y.subvector(0,temp.getHeight()).isZero();
*/
      IntegerVector temp(A.getWidth());
      bool ret=!positiveVectorInKernel(A.getRows(),&temp);
      //debug<<A.getRows()<<temp;
      log1 debug<<(ret?"VECTOR CONFIGURATION IS \"HOMOGENEOUS\"\n":"VECTOR CONFIGURATION IS NOT \"HOMOGENEOUS\"\n");
      return ret;
    }
    /**
     * The support of the secondary fan is given by {x:{y^TA<x}!=emptyset}.
     * Hence the support is the projection of {(x,y):y^TA<x} onto the X component.
     * This function computes the support which is a polyhedral cone.
     */
    PolyhedralCone secondaryFanSupport()
    {
      if(!isHomogeneous())
//    	if(!positiveVectorInKernel(A.getRows(),0))
    	{
    		FieldMatrix A2=integerMatrixToFieldMatrix(A,Q);
    		IntegerVectorList empty;
    		PolyhedralCone C(fieldMatrixToIntegerMatrix(combineLeftRight(FieldMatrix::identity(Q,A2.getWidth()),A2.transposed())).getRows(),empty,A.getHeight()+A.getWidth());
    		C.canonicalize();
    		return C.projection(A2.getWidth());
    	}
    	else
    		return PolyhedralCone(A.getWidth());
    	}
    PolyhedralCone secondaryCone()const
    {
      IntegerVectorList empty;
      return PolyhedralCone(facets(),empty,A.getWidth(),PCP_impliedEquationsKnown);
    }
    IntegerVector interiorPoint()const
    {
      //      IntegerVectorList empty;
      //return PolyhedralCone(inequalities(),empty,A.getWidth()).getRelativeInteriorPoint();
      //return PolyhedralCone(facets(),empty,A.getWidth()).getRelativeInteriorPoint();
      return secondaryCone().getRelativeInteriorPoint();
    }


    IntegerVector determinantInequality3(list<int> const &indices) const
    {
      IntegerMatrix AA(A.getHeight(),indices.size());
      int i=0;
      for(list<int>::const_iterator I=indices.begin();I!=indices.end();I++,i++)
        {
          for(int j=0;j<A.getHeight();j++)
            AA[j][i]=A[j][*I];
        }
      FieldMatrix AAA=integerMatrixToFieldMatrix(AA,Q);
      IntegerVector temp=AAA.reduceAndComputeVectorInKernel().primitive();
      IntegerVector u(A.getWidth());
      i=0;
      for(list<int>::const_iterator I=indices.begin();I!=indices.end();I++,i++)
        u[*I]=-temp[i];
      return u;
    }

    IntegerVector determinantInequality2(list<int> const &indices) const
    {
      IntegerMatrix AA(A.getHeight(),indices.size());
      int i=0;
      for(list<int>::const_iterator I=indices.begin();I!=indices.end();I++,i++)
        {
          for(int j=0;j<A.getHeight();j++)
            AA[j][i]=A[j][*I];
        }
      IntegerVector v=vectorInKernel(AA);
      IntegerVectorList V=AA.getRows();
      V.push_back(v);
      IntegerVector u(A.getWidth());
      i=0;
      for(list<int>::const_iterator I=indices.begin();I!=indices.end();I++,i++)
        u[*I]=v[i];
      if(determinantSign(V)==1)
        return -u;
      return u;
    }


    IntegerVector determinantInequality1(list<int> const &indices) const
    {
      //  fprintf(stderr,"SET:");
      //  for(list<int>::const_iterator i=indices.begin();i!=indices.end();i++)fprintf(stderr," %i",*i);
      //  fprintf(stderr,"\n");
      int m=indices.size();
      FieldVector ret(Q,A.getWidth());

      int I=m;//HER
      for(list<int>::const_iterator i=indices.begin();i!=indices.end();i++,I++)
	{
	  IntegerMatrix B(m-1,A.getHeight());
	  int J=0;
	  for(list<int>::const_iterator j=indices.begin();j!=indices.end();j++)
	    if(i!=j)
	      {
		for(int k=0;k<B.getHeight();k++)
		  {
		    B[k][J]=A[k][*j];
		  }
		J++;
	      }
	  FieldMatrix B2=integerMatrixToFieldMatrix(B,Q);
	  ret[*i]=(Q.zHomomorphism(1-2*(I&1)))*B2.reduceAndComputeDeterminant();
	}
      //      AsciiPrinter T(Stderr);      ret.print(T);
      return ret.primitive();
    }
    IntegerVector determinantInequality(list<int> const &indices) const
    {
      return determinantInequality3(indices);

      IntegerVector r1=determinantInequality1(indices);
      IntegerVector r2=determinantInequality2(indices);
      IntegerVector r3=determinantInequality3(indices);
      debug<<r1<<"\n"<<r2<<"\n"<<r3<<"\n";
      return r3;
    }
    IntegerMatrix subsetRows(IntegerMatrix const &ATransposed, IntegerVector const &cols)const
    {
      IntegerMatrix ret(cols.size(),ATransposed.getWidth());

      for(int i=0;i!=cols.size();i++)ret[i]=ATransposed[cols[i]];
      return ret;
    }
    FieldElement volume(IntegerVector const &v, IntegerMatrix const &ATransposed)const
    {
      IntegerMatrix B(A.getHeight(),A.getHeight());
      for(int j=0;j<v.size();j++)
	B[j]=ATransposed[v[j]];
      FieldMatrix B2=integerMatrixToFieldMatrix(B,Q);
      FieldElement det=B2.reduceAndComputeDeterminant();
      if(det.sign()<0)return -det;
      return det;
    }
    set<set<int> > coDimensionOneTriangles()const
    {
      set<set<int> > codimOne;
      for(set<IntegerVector>::const_iterator j=bases.begin();j!=bases.end();j++)
	{
	  set<int> temp;
	  for(int i=0;i<j->size();i++)
	    {
	      temp.insert((*j)[i]);
	    }
	  for(int i=0;i<j->size();i++)
	    {
	      temp.erase((*j)[i]);
	      codimOne.insert(temp);
	      temp.insert((*j)[i]);
	    }
	}
      return codimOne;
    }
    Graph edgeGraph()const
    {
      Graph ret(bases.size());

      set<set<int> > codimOne=coDimensionOneTriangles();

      for(set<set<int> >::const_iterator i=codimOne.begin();i!=codimOne.end();i++)
	{
	  list<int> triangIndices;
	  int J=0;
	  for(set<IntegerVector>::const_iterator j=bases.begin();j!=bases.end();j++,J++)
	    {
	      bool isSuperSet=true;
	      for(set<int>::const_iterator k=i->begin();k!=i->end();k++)
		{
		  bool hasBeenFound=false;
		  for(int l=0;l<j->size();l++)
		    if((*j)[l]==*k)
		      {
			hasBeenFound=true;
			break;
		      }
		  if(!hasBeenFound)
		    {
		      isSuperSet=false;
		      break;
		    }
		}
	      if(isSuperSet)
		triangIndices.push_back(J);
	    }
	  if(triangIndices.size()==2)
	    ret.addEdge(triangIndices.front(),triangIndices.back());
	}

      //      cerr << ret.toString();

      return ret;
    }
    set<int> difference(IntegerVector const &v, set<int> const &s)const
    {
      set<int> ret;

      for(int i=0;i<v.size();i++)if(s.count(v[i])==0)ret.insert(v[i]);
      return ret;
    }
/** Computes v setminus s - or actually it returns true if this set has size 1 and false other wise.
 *  In case of true the unique element is stored in theDifferencs.
 *  The routine assumes that the number of elements in v is one larger than that of s.
 */
    bool differenceOne(IntegerVector const &v, set<int> const &s, int &theDifference)const
    {
      int diffSize=0;
      for(int i=0;i<v.size();i++)
        if(s.count(v[i])==0)
          {
            diffSize++;
            if(diffSize==2)return false;
            theDifference=v[i];
          }
      return diffSize==1;
    }
    /**
       Seems to compute the outer normals of the secondary cone of the triangulation.
     */

    IntegerVectorList inequalitiesFast()const
    {
      IntegerVectorList ret;

//      cerr<<"bases.size"<<bases.size()<<"A.height"<<A.getHeight()<<endl;

      // we get an inequality for every A-column not a vertex, i.e. not appearing in the triangulation
      for(int i=0;i<A.getWidth();i++)
	{
	  bool appears=false;
	  for(set<IntegerVector>::const_iterator j=bases.begin();j!=bases.end();j++)
	    {
	      for(int k=0;k<j->size();k++)
		if(i==(*j)[k]){appears=true;goto done;}
	    }
	done:
	  if(!appears)
	    {
	      // We now find a triangle which covers the ith column
	      set<IntegerVector>::const_iterator j;
	      int lastSign=0;
	      for(j=bases.begin();j!=bases.end();j++)
		{
		  IntegerMatrix ATransposed=Atransposed;
		  bool containsI=true;
		  IntegerMatrix subMatrix=subsetRows(Atransposed,*j);
		  int sign=determinantSign(subMatrix.getRows());
		  for(int k=0;k<j->size();k++)
		    {
		      subMatrix[k]=ATransposed[i];
		      if(sign*determinantSign(subMatrix.getRows())<0)
			{
			  containsI=false;
			  break;
			}
		      subMatrix[k]=ATransposed[(*j)[k]];
		    }
		  lastSign=sign;
		  if(containsI)break;
		}
	      if(j==bases.end())
		{
		  AsciiPrinter(Stderr).printVector(Atransposed[i]);
		}
	      assert(j!=bases.end());
	      list<int> temp;
	      for(int k=0;k<j->size();k++)temp.push_back((*j)[k]);
	      temp.push_back(i);
	      ret.push_back(lastSign* determinantInequality(temp));//Do we need to use the sign here?//HERE
	    }
	}


      // we get an inequality for every codim one simplex of the triangulation, i.e. every edge of the polyhedron

// TODO: the call to determinantSign below can be avoided if we carefully keep track of orientation
// - that is, the bases should be stored with an orientation flag just as it is done in the trinaglation class.
//      debug<<">\n";
#if 0
  {
      set<set<int> > codimOne;
      for(set<IntegerVector>::const_iterator j=bases.begin();j!=bases.end();j++)
	{
	  set<int> temp;
	  for(int i=0;i<j->size();i++)
	    {
	      temp.insert((*j)[i]);
	    }
	  for(int i=0;i<j->size();i++)
	    {
	      temp.erase((*j)[i]);
	      codimOne.insert(temp);
	      temp.insert((*j)[i]);
	    }
	}
      for(set<set<int> >::const_iterator i=codimOne.begin();i!=codimOne.end();i++)
	{
	  list<int> additional;
	  for(set<IntegerVector>::const_iterator j=bases.begin();j!=bases.end();j++)
	    {
              int theDifference;
              if(differenceOne(*j,*i,theDifference))
                {
                  additional.push_back(theDifference);
                }
	    }
	  if(additional.size()==2)
	    {
	      list<int> temp;
	      for(set<int>::const_iterator j=i->begin();j!=i->end();j++)temp.push_back(*j);
	      list<int>::const_iterator a=additional.begin();

	      list<int> temp2=temp;
	      temp2.push_back(*a);
	      IntegerVector temp3(temp2.size());
	      list<int>::const_iterator K=temp2.begin();
	      for(int k=0;k<temp3.size();k++,K++)temp3[k]=*K;
	      if(determinantSign(subsetRows(Atransposed,temp3).getRows())>0)
		{
		  temp.push_back(*a);
		  a++;
		  temp.push_back(*a);
		}
	      else
		{
		  a++;
		  temp.push_back(*a);
		  a--;
		  temp.push_back(*a);
		}
	      ret.push_back(determinantInequality(temp));//HERE
	    }
	}
  }
#else
  {
      multimap<IntegerVector,pair<int,bool> > codimOne;
      for(set<IntegerVector>::const_iterator j=bases.begin();j!=bases.end();j++)
        {
          bool orientation=0>determinantSign(subsetRows(Atransposed,*j).getRows());
          IntegerVector temp=j->subvector(1,j->size());
          for(int i=0;i<j->size();i++)
            {
              if(i)temp[i-1]=(*j)[i-1];
              IntegerVector temp2=temp;
              int nswaps=mergeSort(temp2)+orientation;
              codimOne.insert(pair<IntegerVector,pair<int,bool> >(temp2,pair<int,bool>((*j)[i],(nswaps+i+j->size())&1)));
            }
        }
      for(multimap<IntegerVector,pair<int,bool> >::const_iterator i=codimOne.begin();i!=codimOne.end();i++)
        {
          multimap<IntegerVector,pair<int,bool> >::const_iterator next=i;next++;
          if(next==codimOne.end())break;
          if(i->first == next->first)
            {
              list<int> temp;
              IntegerVector const &v=i->first;
              for(int j=0;j<v.size();j++)temp.push_back(v[j]);
              if(i->second.second)
                {
                  temp.push_back(i->second.first);
                  temp.push_back(next->second.first);
                }
              else
                {
                  temp.push_back(next->second.first);
                  temp.push_back(i->second.first);
                }
              ret.push_back(determinantInequality(temp));
              i++;
            }
        }
  }
#endif
      return ret;
    }
    IntegerVectorList inequalities()const
    {
      int n=A.getWidth();
      IntegerVectorList ret;
      for(set<IntegerVector>::const_iterator i=bases.begin();i!=bases.end();i++)
	{
	  IntegerMatrix A2transposed(i->size()+1,A.getHeight());
	  for(int j=0;j<i->size();j++)
	    {
	      A2transposed[j]=A.column((*i)[j]);
	    }
	  IntegerVector iComplement=complement(*i,n);

	  for(int k=0;k<iComplement.size();k++)
	    {
	      IntegerVector nyC;

	      A2transposed[i->size()]=A.column(iComplement[k]);
	      IntegerMatrix temp=A2transposed.transposed();
	      nyC=vectorInKernel(temp);


	      IntegerVector c(n);
	      for(int j=0;j<i->size();j++)
		c[(*i)[j]]=nyC[j];
	      c[iComplement[k]]=nyC[i->size()];
	      if(c[iComplement[k]]>0)
		{
		  ret.push_back(c);
		}
	      else
		{
		  ret.push_back(-c);
		}
	    }
	}
      return ret;
    }
    int totalVolume()const
    {
      FieldElement s(Q);
      for(set<IntegerVector>::const_iterator i=bases.begin();i!=bases.end();i++)
	{
	  FieldElement vol=volume(*i,Atransposed);
	  s=s+vol;
	}
      cerr<<"retateawtat"<< s.toString()<<endl;
      return toInteger(s);
    }
    void print(Printer &p)
    {
      p.printString("Bases(");
      p.printInteger(bases.size());
      p.printString(":\n");
      FieldElement s(Q);
      for(set<IntegerVector>::const_iterator i=bases.begin();i!=bases.end();i++)
	{
	  p.printVector(*i);
	  FieldElement vol=volume(*i,Atransposed);
	  p.printString("  Vol: ");
	  p.printFieldElement(vol);
	  p.printNewLine();
	  s=s+vol;
	}
      p.printString("Total volume: ");
      p.printFieldElement(s);
      p.printNewLine();
      /*      p.printString("Circuits (inequalities):\n");
      p.printVectorList(inequalities());
      p.printString("Facets:\n");
      p.printVectorList(facets());
      p.printString("Interior point:\n");
      p.printVector(interiorPoint());
      */
    }
    IntegerVectorList facets()const
    {
      return fastNormals(inequalitiesFast());

      IntegerVectorList flipable;
      IntegerVectorList normals=wallRemoveScaledInequalities(inequalitiesFast());
      //      IntegerVectorList normals=wallRemoveScaledInequalities(inequalities());



      for(IntegerVectorList::iterator i=normals.begin();i!=normals.end();i++)
	{
	  //	  if(!termOrder(*i,*i-*i))
	    {
	      if(isFacet(normals,i))
		{
		  //  if(wallContainsPositiveVector(*i))
		    flipable.push_back(*i);
		}
	      else
		{
		  IntegerVectorList::iterator temp=i;
		  temp++;
		  normals.erase(i);
		  temp--;
		  i=temp;
		}
	    }
	}
      /*      AsciiPrinter Q(Stderr);
      Q.printString("Bases:\n");
      for(set<IntegerVector>::const_iterator i=bases.begin();i!=bases.end();i++)
	{
	  Q.printVector(*i);
	  Q.printNewLine();
	}

      AsciiPrinter(Stderr).printVectorList(inequalities());
      AsciiPrinter(Stderr).printVectorList(flipable);
      AsciiPrinter(Stderr).printVectorList(inequalitiesFast());
      log0 fprintf(stderr,"-----------------\n");
      */

      //log0 fprintf(stderr,"%i\n",0);
      //  log0 fprintf(stderr,"Number of facets:%i\n",flipable.size());

      return flipable;
    }
    bool isEmpty()const
    {
      return bases.empty();
    }
/*    void flip(IntegerVector const &normal)
    {
      AsciiPrinter P(Stderr);
      log2 print(P);
      //log0 P.printVector(normal);
      int n=normal.size();
      //      IntegerVectorList l=wallRemoveScaledInequalities(inequalities());// This is not needed - one circuit should be enough
      IntegerVectorList l;l.push_back(normal);// Let's do this instead
      for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
	if(dependent(*i,normal))
	  {
	    log2 AsciiPrinter(Stderr).printVector(*i);
	    for(int k=0;k<normal.size();k++)
	      if((*i)[k]<0)
		{
		  int s1=bases.size();
		  IntegerVector i2=*i;
		  i2[k]=0;
		  bases.insert(i2.supportIndices());
		  int s2=bases.size();
		  assert(s2==s1+1);
		}
	      else if((*i)[k]>0)
		{
		  int s1=bases.size();
		  IntegerVector i2=*i;
		  i2[k]=0;
		  bases.erase(i2.supportIndices());
		  int s2=bases.size();
		  assert(s2==s1-1);
		}
	  }
    }*/
    static set<int> toSet(IntegerVector const &v)
    {
      set<int> ret;
      for(int i=0;i<v.size();i++)ret.insert(v[i]);
      return ret;
    }
    static IntegerVector toVector(set<int> const &s)
    {
      IntegerVector ret(s.size());
      int I=0;
      for(set<int>::const_iterator i=s.begin();i!=s.end();i++,I++)ret[I]=*i;
      return ret;
    }
    bool contains(IntegerVector const &v, set<int> const &s)
    {
      for(set<int>::const_iterator i=s.begin();i!=s.end();i++)
	{
	  bool found=false;
	  for(int j=0;j!=v.size();j++)
	    if(v[j]==*i)
	      {
		found=true;
		break;
	      }
	  if(!found)return false;
	}
      return true;
    }
    set<int> uni(set<int> const &a, set<int> const &b)
    {
      set<int> ret=a;
      for(set<int>::const_iterator i=b.begin();i!=b.end();i++)ret.insert(*i);
      return ret;
    }
    void printSet(set<int> const &s)
    {
      cerr<<"{";
      for(set<int>::const_iterator i=s.begin();i!=s.end();i++)
	{
	  if(i!=s.begin())cerr<<",";
	  cerr<<*i;
	}
      cerr<<"}";
    }
    void printSetSet(set<set<int> > const &s)
    {
      cerr<<"{"<<endl;
      for(set<set<int> >::const_iterator i=s.begin();i!=s.end();i++)
	{
	  if(i!=s.begin())cerr<<endl<<",";
	  printSet(*i);
	}
      cerr<<"}"<<endl;
    }

    void flipNew(IntegerVector const &normal)
    {
      set<set<int> > toBeInserted;
      set<set<int> > toBeErased;
      int n=normal.size();
      set<int> nonZeroIndices;
      for(int i=0;i<n;i++)
	if(normal[i])nonZeroIndices.insert(i);

      for(set<int>::const_iterator i=nonZeroIndices.begin();i!=nonZeroIndices.end();i++)
	{
	  set<int> temp=nonZeroIndices;
	  temp.erase(*i);
	  if(normal[*i]>0)
	    toBeErased.insert(temp);
	  else
	    toBeInserted.insert(temp);
	}

      assert(toBeErased.size());

      /*      cerr<<"To be inserted"<<endl;
      printSetSet(toBeInserted);

      cerr<<"ToBeErased"<<endl;
      printSetSet(toBeErased);
      */

      set<int> A=*toBeErased.begin();
      set<set<int> > link;
      for(set<IntegerVector>::const_iterator i=bases.begin();i!=bases.end();i++)
	if(contains(*i,A))
	  link.insert(difference(*i,A));

      for(set<IntegerVector>::iterator i=bases.begin();i!=bases.end();)
	{
	  bool e=false;
	  for(set<set<int> >::const_iterator j=toBeErased.begin();j!=toBeErased.end();j++)
	    if(contains(*i,*j))
	      {
		e=true;
		break;
	      }
	  if(e)
	    {
	      set<IntegerVector>::iterator I=i;
	      I++;
	      bases.erase(i);
	      i=I;
	    }
	  else
	    i++;
	}

      for(set<set<int> >::iterator i=link.begin();i!=link.end();i++)
	{
	  for(set<set<int> >::iterator j=toBeInserted.begin();j!=toBeInserted.end();j++)
	    {
	      bases.insert(toVector(uni(*i,*j)));
	    }
	}
    }
    list<int> usedRays()const
    {
      list<int> ret;
      for(int i=0;i<A.getWidth();i++)
	{
	  bool contains=false;
	  for(set<IntegerVector>::const_iterator j=bases.begin();j!=bases.end();j++)
	    {
	      for(int k=0;k<j->size();k++)
		{
		  if((*j)[k]==i)
		    {
		      contains=true;
		      goto leave;
		    }
		}
	    }
	leave:
	  if(contains)ret.push_back(i);
	}
      return ret;
    }
    float hirschScore()const
    {
      int timesAttained;
      int nVertices=bases.size();
      int nEdges=coDimensionOneTriangles().size();
      int diameter=edgeGraph().diameter(&timesAttained);
      int dimension=A.getHeight();
      int nFacets=usedRays().size();

      if(diameter+dimension-nFacets>0)
	{
	  cerr<<"Counter example found\n";
	  //	  print();
	  assert(0);
	}

      return diameter+dimension-nFacets+(nFacets>120?(-0.01*nFacets):0);
    }
    void changeToTriangulationInducedBy(TermOrder const &T)
    {
    	//TODO: use ray shooting to find facet instead
//    	log0 debug<<"changing\n";
    	while(1)
    	{
    		IntegerVectorList l=facets();
    		bool found=false;
    		for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    		{
    			if(!T(*i,*i-*i))
    			{
  //  				debug<<*i;
    				flipNew(-*i);

    				found=true;
    				break;
    			}
    		}
    		if(!found)break;
    	}
    //	log0 debug<<"done changing\n";
    }
    };

#endif
