#include "reversesearch.h"

#include "buchberger.h"
#include "wallideal.h"
#include "printer.h"
#include "lp.h"
#include "log.h"

bool ReverseSearch::computeSearchEdge(PolynomialSet &groebnerBasis, IntegerVector *edge)
{
  //  fprintf(Stderr,"Computing search edge..");

  IntegerVectorList normals=wallInequalities(groebnerBasis);
  normals.sort(LexicographicTermOrder());// Is this needed to make the interior point computation deterministic?

  //  fprintf(Stderr,"\nBasis with no interior point:\n");
  //  AsciiPrinter(Stderr).printPolynomialSet(groebnerBasis);

  if(normals.empty())
    {
      log3 fprintf(Stderr,"WARNING: reversesearch.cpp - No normals\n");
    }

  IntegerVectorList::const_iterator r=shootRay(normals);

  //  fprintf(Stderr,"..done\n");

  if(r!=normals.end())
    {
      *edge=*r;
      return true;
    }
  return false;

  /*  
  IntegerVectorList normals=wallNormals(groebnerBasis);
  normals.sort(LexicographicTermOrder());

  for(IntegerVectorList::const_iterator i=normals.begin();i!=normals.end();i++)
    if(termOrder(*i,*i-*i))
      if(isFacet(normals,i))
	if(wallContainsPositiveVector(*i))
	  {
            *edge=*i;
            return true;
          }
  return false;
  */
}

/*void ReverseSearch::setProgressPrinting(bool p)
{
  progressPrinting=p;
}
*/
static int depth;

//int ReverseSearch::treeSize(const PolynomialSet &groebnerBasis)
int ReverseSearch::treeSize(PolynomialSet &groebnerBasis)
{
  PolynomialRing theRing=groebnerBasis.getRing();
  depth++;
  //  if(progressPrinting)
  {
    static int n;
    n++;
    if(!(n%10))
      log2 fprintf(Stderr,"%i %i\n",n,depth);
  }

  int s=1;
  if(!targetBasis(groebnerBasis)){broken=true;return s;}

  IntegerVectorList flipable;

  //  fprintf(Stderr,"Number of flipable facets:%i\n",flipable.size());
  if(1)
    {
      //      fprintf(Stderr,"isKnownToBeHomogeneous:%i\n",isKnownToBeHomogeneous);
      //  fprintf(Stderr,"Start finding flipable\n");
 
     flipable=wallFlipableNormals(groebnerBasis,isKnownToBeHomogeneous);
     //  fprintf(Stderr,"done\n");
    }
  else
    {
      // For non-homogeneous ideals the following test does not work since it also findes facets that do not intersect the positive orthant.
      assert(isKnownToBeHomogeneous);

      //Taken from Breadth-first search. Apparently this is faster..
      IntegerVectorList normals=algebraicTest(wallInequalities(groebnerBasis),groebnerBasis);
      //      fprintf(Stderr,"Number of inequalities:%i\n",normals.size());
      
      //      AsciiPrinter(Stderr).printVectorList(normals);
      for(IntegerVectorList::iterator i=normals.begin();i!=normals.end();i++)
	{
	  if(!termOrder(*i,*i-*i))
	    {
	      //	      AsciiPrinter(Stderr).printVector(*i);
	      if(isFacet(normals,i))
		{
		  //		  fprintf(Stderr,"isFACET!n");
		  if(wallContainsPositiveVector(*i))
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
    }
  //  AsciiPrinter(Stderr).printVectorList(flipable);

  //  fprintf(Stderr,"Number of flipable facets:%i\n",flipable.size());
  for(IntegerVectorList::iterator i=flipable.begin();i!=flipable.end();i++)
    {
      if(!termOrder(*i,*i-*i))
	{
	  PolynomialSet neighbour=flip(groebnerBasis,*i);
	  
	  IntegerVector edge;
	  if(computeSearchEdge(neighbour,&edge))
	    if(dependent(edge,*i))
	      {
		groebnerBasis=PolynomialSet(theRing);//forget current
		s+=treeSize(neighbour);
		groebnerBasis=flip(neighbour,edge);//recall
		
		if(broken)return s;
	      }
	}
    }


  depth--;
  return s;
}


PolynomialSet ReverseSearch::findRoot(PolynomialSet groebnerBasis)
{
  log2 fprintf(Stderr,"Computing root\n");
  log2 buchberger(&groebnerBasis,termOrder);

  IntegerVector edge;
  while(computeSearchEdge(groebnerBasis,&edge))
    {
      log2 AsciiPrinter(Stderr).printVector(edge);
      groebnerBasis=flip(groebnerBasis,edge);
    }

  log2 fprintf(Stderr,"Done computing root\n");
  return groebnerBasis;
}


ReverseSearch::ReverseSearch(const TermOrder &termOrder_):
  numberOfVertices(0),
  numberOfEdges(0),
  termOrder(termOrder_),
							 isKnownToBeHomogeneous(false)//,
  //  progressPrinting(false)
{
}


void ReverseSearch::enumerate(const PolynomialSet &groebnerBasis)
{
  broken=false;

  PolynomialSet root=findRoot(groebnerBasis);

  if(!isKnownToBeHomogeneous)isKnownToBeHomogeneous=isIdealHomogeneous(root);

  //  fprintf(Stderr,"HOMOGENEOUS:%i\n",isKnownToBeHomogeneous);

  targetBeginEnumeration(groebnerBasis);
  numberOfVertices=treeSize(root);
  targetEndEnumeration();

  //  fprintf(Stderr,"numberOfVertices:%i\n",numberOfVertices);
}
