#include "tropicaltraverse.h"
#include "bergman.h"
#include "enumeration.h"
#include "reversesearch.h"
#include "tropical.h"
#include "buchberger.h"
#include "division.h"
#include "dimension.h"
#include "wallideal.h"
#include "lp.h"
#include "subspace.h"
#include "symmetry.h"
#include "tropical2.h"
#include "tropicalbasis.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "multiplicity.h"
#include "log.h"
#include "restrictedautoreduction.h"
#include "groebnerengine.h"

#include <iostream>

/* Faster version of the code in bergman.cpp. */

/** The hypergraph of ridges and facets can be considered as a usual
    bipartite graph where the right nodes are the ridges and the left
    nodes are the facets.  We wish to make a traversal of this
    bipartite graph keeping track of the boundary edges of the
    traversed set. The ConeOrbit object represents the orbit of a
    ridge. The edges of the ridge are listed but only those which
    belong to the boundary of the set of ridges seen so far. When a
    ridge is discovered the ConeOrbit object will be created with all
    its edges present (except the one it was reached by). As progress
    in the computation is made these edges will be deleted.
 */


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

/**
   Rewrite these comments.

   During traversal the path from the current facet to the starting
   facet is stored on a stack. The elements of the stack are objects
   of the class pathStep. The top element on the stack tells through
   which ridge the current facet was reached. This is done by the
   value parent ridge which is the unique ray of the ridge.  In order
   not to recompute the ridge the path facet contains rays of the link
   of the ridge represented by their unique vector. - or rather only
   the rays that are yet to be processed are stored in ridgeRays. In
   order to trace the path back the unique point of the ray from which
   the ridge was reached is stored in parentRay.
 */
struct pathStepRidge
{
  IntegerVector parentRidge;
  IntegerVectorList rays;
  IntegerVector parentRay;
};

struct pathStepFacet
{
  IntegerVectorList ridges;
};

/**
  We need to simulate two mutually recursive functions. An actual
  implementation of these two functions would propably not work since
  the recursion depth could easily be 10000.

  Here follows a sketch of the simulation

lav kegle
find ridges
skriv ned i objekt
put paa stakken

L1:
if ridges in top element
  compute tropical curve
  construct stak object with rays; set parrentRidge,ridgeRays
  push ridge
else
  pop facet
  if empty break;

goto L2

L2:
if(ridgeRays not empty)
   change CONE
   <---entry point
   push facet
else
  pop ridge
  change CONE

goto L1


The strategy for marking is as follows Before a vertex is pushed the
edges that needs to be taken are written in its data. A edge is only
written if its orbit has not been marked. Each time an edge is written
it is simultaneously marked.

*/

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



static void changeCone(PolynomialSet &coneGroebnerBasis, PolynomialSet &idealGroebnerBasis, IntegerVector const &ridgeVector, IntegerVector const &rayVector)
{
  assert(idealGroebnerBasis.containsInClosedGroebnerCone(ridgeVector));
  log2 cerr<<endl<<"Changing cone"<<endl;

  assert(!containsMonomial(coneGroebnerBasis));

  AsciiPrinter P(Stderr);
  //P<<idealGroebnerBasis;
  //P<<ridgeVector<<rayVector;


  PolynomialSet ridgeIdeal=initialFormsAssumeMarked(idealGroebnerBasis,ridgeVector);
  PolynomialSet ridgeIdealOld=ridgeIdeal;
  WeightReverseLexicographicTermOrder T(rayVector);

  //  P<<ridgeIdeal;
  log2 cerr<<"Computing initial Groebner basis"<<endl;
  //  buchberger(&ridgeIdeal,T);

  ridgeIdeal=GE_groebnerBasis(ridgeIdeal,T,true);
  //printMarkedTermIdeal(ridgeIdeal,"ridgeIdeal");
  coneGroebnerBasis=initialFormsAssumeMarked(ridgeIdeal,rayVector);

  PolynomialSet g2(coneGroebnerBasis.getRing());
  //  WeightTermOrder termOrder(termorderWeight(ridgeIdeal));
  WeightTermOrder termOrder(termorderWeight(ridgeIdealOld));

  log2 cerr<<"Lifting"<<endl;
  PolynomialSet temp=ridgeIdealOld;
  temp.markAndScale(T);
  temp=temp.markedTermIdeal();
  //  P<<temp;
  checkSameLeadingTerms(ridgeIdealOld,idealGroebnerBasis);
  for(PolynomialSet::const_iterator j=ridgeIdeal.begin();j!=ridgeIdeal.end();j++)
    {
      /*      {
	Term m=j->getMarked();
	P.printPolynomial(m);
       	if(division(m,temp,T).isZero()){cerr<<"YES";}
	cerr<<endl;
	}*/

      g2.push_back(divisionLift(*j, ridgeIdealOld, idealGroebnerBasis, termOrder));
      cerr<<"*";
    }
  assert(g2.isMarked());
  //printMarkedTermIdeal(g2,"g2");
  log2 cerr<<"Autoreducing"<<endl;

  //  autoReduce(&g2,LexicographicTermOrder());
  //PolynomialSet g2Old=g2;
  int oldSize=g2.size();
    idealGroebnerBasis=GE_autoReduce(g2);

  //printMarkedTermIdeal(g2,"g2");
  /*  if(g2.size()!=oldSize)
    {
      P<<g2Old;
      P<<g2;
      }*/
  assert(idealGroebnerBasis.size()==oldSize);
  //  idealGroebnerBasis=g2;
  assert(!containsMonomial(coneGroebnerBasis));
  log2 cerr<<"Done changing cone"<<endl<<endl;

  //  P<<coneGroebnerBasis;
  //  P<<idealGroebnerBasis;

  /*  log0{
    WeightReverseLexicographicTermOrder T(ridgeVector);//    WeightTermOrder
    PolynomialSet A=idealGroebnerBasis;
    buchberger(&A,T);
    //cerr<<A;
    PolynomialSet B=initialFormsAssumeMarked(A,ridgeVector);
    WeightReverseLexicographicTermOrder T2(rayVector);
    buchberger(&B,T2);
    cerr<<"RIGHT";
    P<<initialFormsAssumeMarked(B,rayVector);
    cerr<<"RIGHT?";
    P<<coneGroebnerBasis;
    }*/

  log2 cerr << "Number of terms in new basis: "<< g2.totalNumberOfTerms()<<endl;
}

void printStack(list<pathStepFacet> const &facetStack, list<pathStepRidge> const &ridgeStack)
{
  list<pathStepFacet>::const_iterator i=facetStack.begin();
  list<pathStepRidge>::const_iterator j=ridgeStack.begin();
  AsciiPrinter P(Stderr);
  cerr<<"STACK:"<<endl;
  if(facetStack.size()>ridgeStack.size())goto entry;
  do
    {
      cerr<<"RIDGE:"<<endl;
      P<<j->parentRidge<<j->rays<<j->parentRay;
      cerr<<endl;

      j++;
    entry:
      cerr<<"FACET:"<<endl;
      P<<i->ridges;
      cerr<<endl;
      i++;
    }
  while(i!=facetStack.end());

  int a;
  //cin >> a;
}

PolyhedralFan tropicalTraverse(PolynomialSet coneGroebnerBasis, PolynomialSet idealGroebnerBasis, SymmetryGroup const *symmetryGroup)
{
  PolynomialSet ideal=idealGroebnerBasis;

  PolynomialRing theRing=coneGroebnerBasis.getRing();
  assert(coneGroebnerBasis.numberOfVariablesInRing()==idealGroebnerBasis.numberOfVariablesInRing());
  int n=coneGroebnerBasis.numberOfVariablesInRing();

  PolyhedralFan ret(n);

  PolyhedralCone homogeneitySpac=homogeneitySpace(coneGroebnerBasis);
  int d=homogeneitySpac.dimensionOfLinealitySpace();

  SymmetryGroup localSymmetryGroup(n);
  if(!symmetryGroup)symmetryGroup=&localSymmetryGroup;

  Boundary2 boundary(*symmetryGroup);
  list<pathStepFacet> facetStack;
  list<pathStepRidge> ridgeStack;

  int numberOfCompletedFacets=0;
  int numberOfCompletedRidges=0;
  int stackSize=0;

  PolyhedralCone currentFacet(n);
  IntegerVector facetUniqueVector;
  goto entry;
  while(1)
    {
    L1:
      //  boundary.print();
      //printStack(facetStack,ridgeStack);

      if(!facetStack.front().ridges.empty())
	{
	  cerr<<"1";
	  pathStepRidge top;
	  PolyhedralCone link=currentFacet.link(facetStack.front().ridges.front());
	  link.canonicalize();
	  cerr<<"2";
	  top.parentRidge=facetStack.front().ridges.front();
	  top.parentRay=link.getUniquePoint();

	  cerr<<"3";

	  AsciiPrinter P(Stderr);
	  //	  P<<"Cone groebner basis"<<coneGroebnerBasis;
	  // P<<"Ideal groebner basis";
	  //P<<idealGroebnerBasis;
	  assert(idealGroebnerBasis.containsInClosedGroebnerCone(facetStack.front().ridges.front()));
	  PolynomialSet tempIdeal=initialFormsAssumeMarked(idealGroebnerBasis,facetStack.front().ridges.front());

	  //P<<tempIdeal;
	  tempIdeal=saturatedIdeal(tempIdeal);
	  IntegerVectorList rays;

	  PolyhedralCone saturatedHomogeneitySpace=homogeneitySpace(tempIdeal);
	  if(saturatedHomogeneitySpace.dimensionOfLinealitySpace()==d)//if saturation of the initial ideal changed the homogeneity space things are easy
	    {
	      rays.push_back(top.parentRay);
	      rays.push_back(-top.parentRay);
	    }
	  else
	    {
	      //P<<tempIdeal;
	      BergmanFan b=bergmanRayIntersection(tempIdeal);
	      //	  BergmanFan b=bergmanRayIntersection(initialIdeal(idealGroebnerBasis,facetStack.front().ridges.front()));


	      cerr<<"4";
	      bool trouble=false;
	      for(BergmanFan::MaximalConeList::const_iterator i=b.cones.begin();i!=b.cones.end();i++)
		{
		  PolyhedralCone rayCone=i->theCone;
		  rayCone.canonicalize();
		  {
		    if(rayCone.getUniquePoint().isZero())trouble=true;
		  }
		  rays.push_back(rayCone.getUniquePoint());
		}
	      if(trouble)
		{
		  b.print(P);
		  P<<tempIdeal;
		  assert(!trouble);
	      }
	    }


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
    L2:
      //      boundary.print();
      //printStack(facetStack,ridgeStack);
      if(!ridgeStack.front().rays.empty())
	{
	  changeCone(coneGroebnerBasis, idealGroebnerBasis, ridgeStack.front().parentRidge,ridgeStack.front().rays.front());
	entry:
	  //	  ret.cones.push_back(BergmanFan::MaximalCone(coneGroebnerBasis,idealGroebnerBasis,false,0,0));

	  log1 fprintf(Stderr,"\n-------------------------------------\n");
	  log1 fprintf(Stderr,"Boundary edges in bipartite graph: %i, Completed ridges: %i, Completed facets: %i, Recursion depth:%i\n",boundary.size(),numberOfCompletedRidges,numberOfCompletedFacets,stackSize);
	  log1 fprintf(Stderr,"-------------------------------------\n");

	  currentFacet=PolyhedralCone(fastNormals(wallInequalities(idealGroebnerBasis)),wallInequalities(coneGroebnerBasis),n);
	  cerr<<"A";
	  currentFacet.canonicalize();
	  cerr<<"B";

	  ret.insert(currentFacet);

	  IntegerVectorList inequalities=currentFacet.getHalfSpaces();
	  IntegerVectorList equations=currentFacet.getEquations();
	  facetUniqueVector=currentFacet.getUniquePoint();
	  IntegerVectorList facetNormals=currentFacet.getHalfSpaces();

	  pathStepFacet stepFacet;
	  IntegerVectorList ridges;
	  cerr<<"C";
	  for(IntegerVectorList::const_iterator i=facetNormals.begin();i!=facetNormals.end();i++)
	    {
	      equations.push_back(*i);
	      PolyhedralCone ridgeCone(inequalities,equations,n);
	      equations.pop_back();
	      ridgeCone.canonicalize();
	      ridges.push_back(ridgeCone.getUniquePoint());
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
    }//goto L1
  return ret;
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
