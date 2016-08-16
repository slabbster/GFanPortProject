#include "regularsubdivision.h"
#include "polyhedralcone.h"
#include "log.h"
#include "polyhedralfan.h"

set<set<int> > regularSubdivision(IntegerMatrix const &m, IntegerVector const &w)
{
  IntegerVectorList m2;
  
  int n=m.getWidth();
  
  for(int i=0;i<m.getHeight();i++)
    {
      IntegerVector v=m[i];
      v.grow(v.size()+1);
      v[v.size()-1]=w[i];
      m2.push_back(v);
    }
  
  IntegerVectorList empty;
  PolyhedralCone c(m2,empty,n+1);
  c.canonicalize();
  /*
  log0
    {
      AsciiPrinter P(Stderr);
      c.print(&P);

      PolyhedralFan F(n+1);F.insert(c);
      F.printWithIndices(&P,false,0,false,false);

    }
  */
  
  //  log0 fprintf(Stderr,"A");
  IntegerVectorList lowerNormals=c.extremeRays();
  //  log0 fprintf(Stderr,"B");
  
  //    AsciiPrinter(Stdout).printVectorList(m2);
  //    AsciiPrinter(Stdout).printVectorList(lowerNormals);
  
  set<set<int> > faces;
  for(IntegerVectorList::const_iterator i=lowerNormals.begin();i!=lowerNormals.end();i++)
    {
      set<int> face;
       if((*i)[n]<0) //This is how it was originally
	 //if((*i)[n]>0)
	{
	  int J=0;
	  for(IntegerVectorList::const_iterator j=m2.begin();j!=m2.end();j++,J++)
	    if(dotLong(*j,*i)==0)
	      face.insert(J);
	  faces.insert(face);
	}
    }
  return faces;
}

void printSetSetInt(FILE *f, set<set<int> > const &faces)
  {
    fprintf(f,"{");
    for(set<set<int> >::const_iterator i=faces.begin();i!=faces.end();i++)
      {
	if(i!=faces.begin())fprintf(f,",\n");
	
	fprintf(f,"{");
	for(set<int>::const_iterator j=i->begin();j!=i->end();j++)
	  {
	    if(j!=i->begin())fprintf(f,",");
	    fprintf(f,"%i",*j);
	  }
	fprintf(f,"}");
      }
    fprintf(f,"}\n");
  }


/*
PolyhedralCone secondaryCone(IntegerMatrix const &m, set<set<int> > const &subdivision)
{
  IntegerVectorList equations;
  IntegerVectorList inequalities;

  for(set<set<int> >::const_iterator i=subdivision.begin();i!=subdivision.end();i++)
    {
      for(int j=0;j<matrix.getHeight();j++)
	{
	  IntegerVector inequality(matrix.getHeight());
	  if(!i->contains(j))
	    {
	      for(set<int>::const_iterator k=i->begin();k!=i->end();k++)
		{
		  
		  IntegerMatrix m(m.getWidth(),m.getWidth());
		  int L=0;
		  for(set<int>::const_iterator l=i->begin();l!=i->end();l++)
		    if(l!=k)
		      {
			m[L]=*L;
			L++;
		      }
		  
		}
	    }
	}
    }

}

*/
