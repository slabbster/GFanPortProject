#include "symmetrictraversal.h"

#include <map>
#include <algorithm>
#include <iostream>

#include "log.h"

using namespace std;

/**
 * SymmetricTargetFanBuilder
 */

SymmetricTargetFanBuilder::SymmetricTargetFanBuilder(int n, SymmetryGroup const &sym):
	coneCollection(n)
{
}


bool SymmetricTargetFanBuilder::process(ConeTraverser &traverser)
{
	PolyhedralCone cone2=traverser.refToPolyhedralCone();
	cone2.canonicalize();
	coneCollection.insert(cone2);
	return true;
}




/**
 * Classes
 */

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


class Boundary
{
  typedef pair<IntegerVector,IntegerVector> EFirst;
  typedef pair<IntegerVectorList*,IntegerVectorList::iterator> ESecond;
  SymmetryGroup const &sym;
  map<EFirst,ESecond > theSet;
  int theSetSize;
public:
  Boundary(SymmetryGroup const &sym_):
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
  /**
   * This routine remove rays from rays, such that only one ridge-ray pair is left for each orbit.
   * The routine allows an additional list of vectors with the same number of elements as rays to be passed.
   * The routine will remove those vectors from this set which correspond to rays removed from rays.
   *
   * To do this it must know the symmetry group.
   */
  void removeDuplicates(IntegerVector const &ridge, IntegerVectorList &rays, IntegerVectorList *normals=0)const
  {
    IntegerVectorList ret;
    IntegerVectorList normalsRet;
    set<IntegerVector> representatives;
    IntegerVectorList::const_iterator I;
    if(normals)I=normals->begin();
    for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
      {
	IntegerVector rep=sym.orbitRepresentativeFixing(*i,ridge);
	if(representatives.count(rep)==0)
	  {
	    representatives.insert(rep);
	    ret.push_back(*i);
	    if(normals)normalsRet.push_back(*I);
	  }
	if(normals)I++;
      }
    rays=ret;
    if(normals)*normals=normalsRet;
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
  implementation of these two functions would probably not work since
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

/*static void printMarkedTermIdeal(PolynomialSet const &g, string const &s)
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
}*/




static void printStack(list<pathStepFacet> const &facetStack, list<pathStepRidge> const &ridgeStack)
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



void symmetricTraverse(ConeTraverser &traverser, SymmetricTarget &target, SymmetryGroup const *symmetryGroup)
{
	int lastNumberOfEdges=0;
	float averageEdge=0;
	int n=traverser.refToPolyhedralCone().ambientDimension();//symmetryGroup->sizeOfBaseSet();
	  SymmetryGroup localSymmetryGroup(n);
	  if(!symmetryGroup)symmetryGroup=&localSymmetryGroup;

	IntegerVectorList linealitySpaceGenerators=traverser.refToPolyhedralCone().generatorsOfLinealitySpace();

	int d=traverser.refToPolyhedralCone().dimension();

	  Boundary boundary(*symmetryGroup);
	  list<pathStepFacet> facetStack;
	  list<pathStepRidge> ridgeStack;

	  int numberOfCompletedFacets=0;
	  int numberOfCompletedRidges=0;
	  int stackSize=0;

	  IntegerVector facetUniqueVector;
	  goto entry;
	  while(1)
	    {
	    L1:
//	    printStack(facetStack,ridgeStack);
	    //if we have more ProcessRidge calls to do
	      if(!facetStack.front().ridges.empty())
		{
	    	  //ProcessRidge "called"
		  pathStepRidge top;



		  PolyhedralCone link=traverser.refToPolyhedralCone().link(facetStack.front().ridges.front());
		  link.canonicalize();
		  top.parentRay=link.getUniquePoint();


		  top.parentRidge=facetStack.front().ridges.front();
//		  AsciiPrinter(Stderr)<<top.parentRay<<"--------------------------------++++\n";
		  IntegerVectorList rays=traverser.link(facetStack.front().ridges.front());
		  boundary.removeDuplicates(top.parentRidge,rays);
		  ridgeStack.push_front(top);stackSize++;
		  IntegerVector ridgeRidgeRidge=facetStack.front().ridges.front();
		  for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
		    {
		      ridgeStack.front().rays.push_front(*i);
		      if(boundary.containsFlip(ridgeRidgeRidge,*i,&ridgeStack.front().rays,ridgeStack.front().rays.begin()))
			ridgeStack.front().rays.pop_front();
		    }
		  // "state saved" ready to do calls to ProcessFacet
		  numberOfCompletedRidges++;
		}
	      else
		{
	    	  // No more calls to do - we now return from ProcessFacet
	    	  //THIS IS THE PLACE TO CHANGE CONE BACK
		  facetStack.pop_front();stackSize--;
		  if(facetStack.empty())break;
    	 log1 cerr<<"BACK"<<endl;
	  traverser.changeCone(ridgeStack.front().parentRidge,ridgeStack.front().parentRay);
		}
	    L2:
//	    printStack(facetStack,ridgeStack);
	    //check if ProcessRidge needs to make more ProcessFacet calls
	      if(!ridgeStack.front().rays.empty())
		{
	    	  log1 cerr<<"FORWARD"<<endl;
	    	  traverser.changeCone(ridgeStack.front().parentRidge,ridgeStack.front().rays.front());
		entry:
		//ProcessFacet()
		averageEdge=0.99*averageEdge+0.01*(boundary.size()-lastNumberOfEdges);
		log1 fprintf(Stderr,"\n-------------------------------------\n");
		  log1 fprintf(Stderr,"Boundary edges in bipartite graph: %i, Completed ridges: %i, Completed facets: %i, Recursion depth:%i Average new edge/facet:%0.2f\n",boundary.size(),numberOfCompletedRidges,numberOfCompletedFacets,stackSize,averageEdge);
		  log1 fprintf(Stderr,"-------------------------------------\n");
		  lastNumberOfEdges=boundary.size();

//		  target.process(traverser);//Postponed until extrem rays have been computed
		  IntegerVectorList extremeRays=traverser.refToPolyhedralCone().extremeRays(&linealitySpaceGenerators);
                  target.process(traverser);

//		  IntegerVectorList inequalities=traverser.refToPolyhedralCone().getHalfSpaces();
		  IntegerVectorList equations=traverser.refToPolyhedralCone().getEquations();
//		  facetUniqueVector=traverser.refToPolyhedralCone().getUniquePoint();
		  facetUniqueVector=traverser.refToPolyhedralCone().getUniquePointFromExtremeRays(extremeRays);
		  IntegerVectorList facetNormals=traverser.refToPolyhedralCone().getHalfSpaces();

		  pathStepFacet stepFacet;
		  IntegerVectorList ridges;

		  for(IntegerVectorList::iterator i=facetNormals.begin();i!=facetNormals.end();i++)
		    {
			  if(1)
			  {
				  IntegerVector v(n);
				  for(IntegerVectorList::const_iterator j=extremeRays.begin();j!=extremeRays.end();j++)
					  if(dotLong(*i,*j)==0)v+=*j;
				  ridges.push_back(v);
			  }
			  else
			  {
		      equations.push_back(*i);
//		      PolyhedralCone ridgeCone(inequalities,equations,n);


		    	  IntegerVectorList::iterator I=i;
		    	  i++;
		    	  facetNormals.erase(I);
//		    	  PolyhedralCone ridgeCone(facetNormals,equations,n);
		    	  PolyhedralCone ridgeCone =PolyhedralCone::polyhedralConeWithKnownImpliedEquations(facetNormals,equations,n);
		    	  facetNormals.insert(i,equations.back());
		    	  i--;

		      equations.pop_back();
		      ridgeCone.canonicalize();
//		      ridges.push_back(ridgeCone.getUniquePoint());
		      ridges.push_back(ridgeCone.getUniquePointFromExtremeRays(extremeRays));
//		      debug << *i <<"\n"<<ridgeCone.getUniquePointFromExtremeRays(extremeRays)<<"\n\n";
			  }
		    }

		  IntegerVector temp(n);
//		  boundary.removeDuplicates(temp,ridges);//use facetUniqueVector instead
		  boundary.removeDuplicates(facetUniqueVector,ridges,&facetNormals);//use facetUniqueVector instead

		  facetStack.push_front(stepFacet);stackSize++;
		  IntegerVectorList::const_iterator I=facetNormals.begin();
		  for(IntegerVectorList::const_iterator i=ridges.begin();i!=ridges.end();i++,I++)
		    {
			  IntegerVector rayUniqueVector;

			  if(d==n)
			  {
				rayUniqueVector =normalized(*I);
//				if(dotLong(rayUniqueVector,*I)
			  }
			  else
			  {
				  PolyhedralCone rayCone=traverser.refToPolyhedralCone().link(*i);
				  rayCone.canonicalize();
				  rayUniqueVector=rayCone.getUniquePoint();
//				  debug<<traverser.refToPolyhedralCone();
//				  debug<<rayCone;
			  }
/*if(0)
			  if(!(rayUniqueVector-normalized(*I)).isZero())
				  if(!(rayUniqueVector+normalized(*I)).isZero())

				  {
					  debug << "extreme rays" << extremeRays;
					  debug << "facet normal" << *I;
					  debug << "products "<<rowsToIntegerMatrix(extremeRays).vectormultiply(*I);
					  debug << "ridge pt "<<*i;
				      equations.push_back(*I);
				      PolyhedralCone ridgeCone(inequalities,equations,n);
				      equations.pop_back();
				      ridgeCone.canonicalize();
					  debug<<ridgeCone.getUniquePointFromExtremeRays(extremeRays);
					  debug << "CONE"<<traverser.refToPolyhedralCone();
					  debug << "LINK"<<traverser.refToPolyhedralCone().link(*i);

					  debug << rayUniqueVector<<"\n";
			  debug << normalized(*I)<<"\n\n";
			  assert(0);
				  }
				  */
		      facetStack.front().ridges.push_front(*i);
		      if(boundary.containsFlip(*i,rayUniqueVector,&facetStack.front().ridges,facetStack.front().ridges.begin()))
			facetStack.front().ridges.pop_front();
		    }
		  //"State pushed" ready to call ProcessRidge

		  numberOfCompletedFacets++;
		}
	      else
		{
	    	  //ProcessRidge is done making its calls to ProcessFacet so we can return from ProcessRidge
//	    	  cerr<<"BACK"<<endl;
//		  traverser.changeCone(ridgeStack.front().parentRidge,ridgeStack.front().parentRay);
		  ridgeStack.pop_front();stackSize--;

/*		  for(IntegerVectorList::const_iterator i=facetStack.front().ridges.begin();i!=facetStack.front().ridges.end();i++)
		    {
		      assert(idealGroebnerBasis.containsInClosedGroebnerCone(*i));
		      assert(coneGroebnerBasis.isHomogeneous(*i));
		    }*/
		}
	    }//goto L1
	  log1 fprintf(Stderr,"\n-------------------------------------\n");
	  log1 fprintf(Stderr,"Boundary edges in bipartite graph: %i, Completed ridges: %i, Completed facets: %i, Recursion depth:%i\n",boundary.size(),numberOfCompletedRidges,numberOfCompletedFacets,stackSize);
	  log1 fprintf(Stderr,"-------------------------------------\n");
//	  return ret;
}
