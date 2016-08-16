// To do: The facet detection in this file may not work for non-homogeneous ideals since facets pass the test even if they do not contain a positive vector. Check if this is the case and fix the bug.

#include "breadthfirstsearch.h"

#include "buchberger.h"
#include "wallideal.h"
#include "printer.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "lp.h"
#include "log.h"

#include <iostream>


class Boundary2
{
  typedef pair<IntegerVector,IntegerVector> EFirst;
  typedef pair<IntegerVectorList*,IntegerVectorList::iterator> ESecond;
  SymmetryGroup const &sym;
  map<EFirst,ESecond > theSet;
  int theSetSize;
public:
  Boundary2(SymmetryGroup const &sym_):
    sym(sym_),
    theSetSize(0)
  {
  }
  int size()const
  {
    return theSetSize;
  }
  pair<IntegerVector,IntegerVector> normalForm(IntegerVector const &ridge, IntegerVector const &ray)const
  {
    pair<IntegerVector,IntegerVector> ret;
    IntegerVector perm;
    ret.first=sym.orbitRepresentative(ridge,&perm);
    ret.second=sym.orbitRepresentativeFixing(SymmetryGroup::compose(perm,ray),ret.first);
    return ret;
  }
  bool containsFlip(IntegerVector const &ridge, IntegerVector const &ray, IntegerVectorList *storedInList_, IntegerVectorList::iterator listIterator_)
  {
    assert(ridge.size()==ray.size());
    EFirst p=normalForm(ridge,ray);
    if(theSet.count(p)==1)
      {
	theSet[p].first->erase(theSet[p].second);
	theSet.erase(p);
	theSetSize--;
	return true;
      }
    theSet[p]=ESecond(storedInList_,listIterator_);
    theSetSize++;
    return false;
  }
  void removeDuplicates(IntegerVector const &ridge, IntegerVectorList &rays)const
  {
    IntegerVectorList ret;
    set<IntegerVector> representatives;
    for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
      {
	IntegerVector rep=sym.orbitRepresentativeFixing(*i,ridge);
	if(representatives.count(rep)==0)
	  {
	    representatives.insert(rep);
	    ret.push_back(*i);
	  }
      }
    rays=ret;
  }
  void print()const
  {
    cerr<< "Boundary" <<endl;
    for(map<EFirst,ESecond>::const_iterator i=theSet.begin();i!=theSet.end();i++)
      {
	AsciiPrinter P(Stderr);
	P << i->first.first << i->first.second;
	cerr << endl;
      }
    cerr<<endl<<endl;
  }
};


struct pathStepFacet
{
  IntegerVectorList ridges;
};






static void checkSameLeadingTerms(PolynomialSet const &a, PolynomialSet const &b)
{
  assert(a.size()==b.size());
  PolynomialSet::const_iterator A=a.begin();
  for(PolynomialSet::const_iterator B=b.begin();B!=b.end();B++,A++)
    assert(A->getMarked().m.exponent==B->getMarked().m.exponent);
}

static void printMarkedTermIdeal(PolynomialSet const &g, string const &s)
{
  PolynomialSet a=g.markedTermIdeal();
  PolynomialSet b=a;
  minimize(&b);
  cerr << "Printing marked termideal. "<<s<< "size:"<<g.size()<<","<<a.size()<<","<<b.size()<<endl;
  AsciiPrinter P(Stderr);
  cerr << "initial ideal:";
  P<<a;
  cerr << "minimized:";
  P<<b;
  assert(a.size()==b.size());
}




void printStack(list<pathStepFacet> const &facetStack)
{
  list<pathStepFacet>::const_iterator i=facetStack.begin();
  AsciiPrinter P(Stderr);
  cerr<<"STACK:"<<endl;
  do
    {
      cerr<<"FACET:"<<endl;
      P<<i->ridges;
      cerr<<endl;
      i++;
    }
  while(i!=facetStack.end());
}




void coneChangeDebugger(PolynomialSet  coneGroebnerBasis, PolynomialSet idealGroebnerBasis, IntegerVectorList const &ridges, IntegerVectorList const &rays)
{
  AsciiPrinter P(Stderr);
  P<<coneGroebnerBasis<<idealGroebnerBasis;
  IntegerVectorList::const_iterator i=ridges.begin();
  for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++,i++)
    {


      changeCone(coneGroebnerBasis, idealGroebnerBasis,*i,*j);
      P<<"NEW CONEGB:"<<coneGroebnerBasis;
      P<<"NEW IDEALGB:"<<idealGroebnerBasis;

    }
}



BreadthFirstSearch::BreadthFirstSearch(const SymmetryGroup &symmetryGroup_, bool minkowski_):
  numberOfVertices(0),
  numberOfEdges(0),
  symmetryGroup(symmetryGroup_),
  minkowski(minkowski_)
{
}

void BreadthFirstSearch::setSubspace(IntegerVectorList const &subspacePerp_)
{
  subspacePerp=subspacePerp_;
}


void BreadthFirstSearch::enumerate(const PolynomialSet &groebnerBasis)
{
  int n=groebnerBasis.getRing().getNumberOfVariables();

  targetBeginEnumeration(groebnerBasis);

  PolynomialSet current=groebnerBasis;

  SymmetryGroup localSymmetryGroup(n);
  if(!symmetryGroup)symmetryGroup=&localSymmetryGroup;

  Boundary boundary(*symmetryGroup);
  list<pathStepFacet> facetStack;

  int numberOfCompletedFacets=0;
  int numberOfCompletedRidges=0;
  int stackSize=0;

  PolyhedralCone currentFacet(n);
  IntegerVector facetUniqueVector;
  goto entry;
  while(1)
    {
      if(!facetStack.front().ridges.empty())
	{
	  IntegerVector normal=facetStack.front().ridges.front();

	  current=(minkowski)?flipMinkowski(g,normal):flip(g,normal);


	  //P<<"RIDGE"<<facetStack.front().ridges.front()<<"\nRAYS"<<rays;
	  cerr<<"5";

	  boundary.removeDuplicates(top.parentRidge,rays);

	  cerr<<"6";
	  ridgeStack.push_front(top);stackSize++;
	  IntegerVector ridgeRidgeRidge=facetStack.front().ridges.front();
	  for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
	    {
	      ridgeStack.front().rays.push_front(*i);
	      if(boundary.containsFlip(ridgeRidgeRidge,*i,&ridgeStack.front().rays,ridgeStack.front().rays.begin()))
		ridgeStack.front().rays.pop_front();
	    }
	  cerr<<"7";
	  numberOfCompletedRidges++;
	}
      else
	{
	  facetStack.pop_front();stackSize--;
	  if(facetStack.empty())break;
	}
    entry:
      log1 fprintf(Stderr,"\n-------------------------------------\n");
      log1 fprintf(Stderr,"Boundary edges in bipartite graph: %i, Completed ridges: %i, Completed facets: %i, Recursion depth:%i\n",boundary.size(),numberOfCompletedRidges,numberOfCompletedFacets,stackSize);
      log1 fprintf(Stderr,"-------------------------------------\n");
	  
      currentFacet=PolyhedralCone(fastNormals(wallInequalities(current)),n);
      cerr<<"A";
      currentFacet.canonicalize();
      cerr<<"B";
      
      if(!targetBasis(g))break;

      facetUniqueVector=currentFacet.getUniquePoint();
      IntegerVectorList facetNormals=currentFacet.getHalfSpaces();
	  
      pathStepFacet stepFacet;
      IntegerVectorList ridges;
      cerr<<"C";
      for(IntegerVectorList::const_iterator i=facetNormals.begin();i!=facetNormals.end();i++)
	{
	  ridges.push_back(ridgeCone.normalized(*i));
	}
      cerr<<"D";
      IntegerVector temp(n);
      boundary.removeDuplicates(temp,ridges);
      
      cerr<<"E";
      facetStack.push_front(stepFacet);stackSize++;
      for(IntegerVectorList::const_iterator i=ridges.begin();i!=ridges.end();i++)
	{
	  PolyhedralCone rayCone=currentFacet.link(*i);
	  rayCone.canonicalize();
	  IntegerVector rayUniqueVector=rayCone.getUniquePoint();
	  facetStack.front().ridges.push_front(*i);
	  if(boundary.containsFlip(*i,rayUniqueVector,&facetStack.front().ridges,facetStack.front().ridges.begin()))
	    facetStack.front().ridges.pop_front();
	}
      cerr<<"F";
      numberOfCompletedFacets++;
    }

  else
	{
	  changeCone(coneGroebnerBasis, idealGroebnerBasis, ridgeStack.front().parentRidge,ridgeStack.front().parentRay);
	  currentFacet=PolyhedralCone(fastNormals(wallInequalities(idealGroebnerBasis)),wallInequalities(coneGroebnerBasis),n);
	  currentFacet.canonicalize();

	  ridgeStack.pop_front();stackSize--;
	  cerr<<"BACK"<<endl;
	  for(IntegerVectorList::const_iterator i=facetStack.front().ridges.begin();i!=facetStack.front().ridges.end();i++)
	    {
	      assert(idealGroebnerBasis.containsInClosedGroebnerCone(*i));
	      assert(coneGroebnerBasis.isHomogeneous(*i));
	    }


	}





	    }
	}
    }

  targetEndEnumeration();
}
