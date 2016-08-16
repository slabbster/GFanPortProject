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

static Polynomial wallPolynomial(Polynomial const &p, Subspace const &subspace)
  // This routine should be deleted since it is inefficient
  // hmm... does this actually work?
{
  Polynomial r(p.getRing());
  IntegerVector markedExponent=p.getMarked().m.exponent;

  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      IntegerVector dif=markedExponent-i->first.exponent;

      if(subspace.contains(dif))
	r+=Polynomial(Term(i->second,i->first));
    }

  r.mark(Monomial(p.getRing(),markedExponent));

  return r;
}

static PolynomialSet wallIdeal(PolynomialSet const &groebnerBasis, Subspace const &subspace)
{
  //  fprintf(Stderr,"wallIdeal %i",perp.size());
  PolynomialSet ret(groebnerBasis.getRing());
  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    ret.push_back(wallPolynomial(*i,subspace));

  //  fprintf(Stderr,"done\n");
  return ret;
}


BergmanFan bergmanRayIntersection(PolynomialSet const &idealGroebnerBasis)
  // Input ideal is assumed to be homogeneous with respect to a positive vector
  // Input ideal is assumed not to contain a monomial
  // Call the Krull dimensionn of the ring/ideal d
  // The ideal is homogenous with respect to a vector space of dimension d-1
  // All d dimensional cones in the Gro\"obner fan with no monomials are computed
{
  BergmanFan bfan;

  int n=idealGroebnerBasis.numberOfVariablesInRing();

  //  log0 AsciiPrinter(Stderr).printPolynomialSet(idealGroebnerBasis);


  int d=krullDimension(idealGroebnerBasis); //There should be a better way of doing this
  PolyhedralFan bergmanFan(n);
  PolynomialSet tropBasis=tropicalBasisOfCurve(n,idealGroebnerBasis,&bergmanFan,d-1);
  //  PolyhedralFan bergmanFan=tropicalPrincipalIntersection(n,tropBasis,d-1);
  IntegerVectorList rays=bergmanFan.getRays(d);

  int maximalConeLabel=0;

  //  fprintf(Stderr,"---------------------------------------------------------\n");
  //  AsciiPrinter(Stderr).printVectorList(rays);
  //  fprintf(Stderr,"---------------------------------------------------------\n");


  log2  cerr<<"BBLABL";
  for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
    {
      log2 cerr<<"1";
      PolynomialSet g=idealGroebnerBasis;
      g=GE_groebnerBasis(g,WeightReverseLexicographicTermOrder(*i),true);
      //      buchberger(&g,WeightReverseLexicographicTermOrder(*i));
      PolynomialSet cg=initialFormsAssumeMarked(g,*i);
      log2 cerr<<"2";
 
      bool inList=false;
      for(BergmanFan::MaximalConeList::const_iterator j=bfan.cones.begin();j!=bfan.cones.end();j++)
	{
	  if(areIdealsEqual(j->coneGroebnerBasis,cg))
	    {
	      inList=true;
	      break;
	    }
	}
      log2 cerr<<"3";
      if(!inList)
	{
	  bfan.cones.push_back(BergmanFan::MaximalCone(cg,g,true,maximalConeLabel++));
	}
      log2 cerr<<"4";
    }
  log2 cerr<<"BBLAB  END";
  //  AsciiPrinter temp(Stderr);
  //  bfan.print(temp);
  return bfan;
}

BergmanFan bergmanRay(PolynomialSet const &idealGroebnerBasis)
  // Input ideal is assumed to be homogeneous with respect to a positive vector
  // Input ideal is assumed not to contain a monomial
  // Call the Krull dimensionn of the ring/ideal d
  // The ideal is homogenous with respect to a vector space of dimension d-1
  // All d dimensional cones in the Gro\"obner fan with no monomials are computed
{
  BergmanFan bfan;

  EnumerationTargetCollector gfan;

  LexicographicTermOrder myTermOrder;
  ReverseSearch rs(myTermOrder);
  rs.setEnumerationTarget(&gfan);
  rs.enumerate(idealGroebnerBasis);

  int n=idealGroebnerBasis.numberOfVariablesInRing();
  fprintf(Stderr,"rankOfMatrix(wallin.idealGroebnerBasis=%i\n",rankOfMatrix(wallInequalities(idealGroebnerBasis)));
  AsciiPrinter(Stderr).printVectorList(wallInequalities(idealGroebnerBasis));
  int d=n-rankOfMatrix(wallInequalities(idealGroebnerBasis))+1;
  int krull=krullDimension(idealGroebnerBasis);

  //AsciiPrinter(Stderr).printVectorList(wallInequalities(idealGroebnerBasis));
  assert(rankOfMatrix(wallInequalities(idealGroebnerBasis))==Subspace(wallInequalities(idealGroebnerBasis)).dimension());

  //  fprintf(Stderr,"d: %i krull: %i\n",d,krull);
  //  assert(d==krull);

  int maximalConeLabel=0;

  for(PolynomialSetList::const_iterator g=gfan.theList.begin();g!=gfan.theList.end();g++)
    {      
      PolynomialSetList s;
      //fprintf(Stderr,"current gbasis:\n");
      //AsciiPrinter(Stderr).printPolynomialSet(*g);
      if(0)
	{
	  s=fullColoredIdeals(*g,false);
	  fprintf(Stderr,"Full colored ideals computed, #=%i\n",s.size());
	}
      else
	{
	  IntegerVectorList inequalities=wallInequalities(*g);
	  inequalities=wallFlipableNormals(*g,true);
	  int isize=inequalities.size();
	  //	  fprintf(Stderr,"cdd facets to rays ");
	  // AsciiPrinter(Stderr).printVectorList(inequalities);
	  IntegerVectorList rays=extremeRaysInequalityIndices(inequalities);
	  //fprintf(Stderr,"done\n");
	  //AsciiPrinter(Stderr).printVectorList(rays);
	  //  AsciiPrinter(Stderr).printVectorList(rays);
	  for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
	    if(i->size()!=isize)
	      {
		IntegerVectorList perp;
		
		int j=0;
		int K=0;
		for(IntegerVectorList::const_iterator k=inequalities.begin();k!=inequalities.end()&&j<i->size();k++,K++)
		  if((*i)[j]==K)
		    {
		      perp.push_back(*k);
		      j++;
		    }
		
		s.push_back(wallIdeal(*g,Subspace(perp)));
	      }
	}

      //fprintf(Stderr,"Number of face ideals to check:%i\n",s.size());

      for(PolynomialSetList::const_iterator i=s.begin();i!=s.end();i++)
	{
	  //fprintf(Stderr,"d:%i\n",d);
	  //fprintf(Stderr,"krull:%i\n",krull);
	  //fprintf(Stderr,"n:%i\n",n);
	  //fprintf(Stderr,"rank: %i\n",rankOfMatrix(wallInequalities(*i)));
	  assert(i->isMarked());
	  //	  fprintf(Stderr,"d: %i, rank %i",d,n-rankOfMatrix(wallInequalities(*i)));
	if(d==n-rankOfMatrix(wallInequalities(*i)))//dimension check
	  {
	    //  fprintf(Stderr,"Checking monomial containment");
	  if(!containsMonomial(*i))
	    {
	      //  fprintf(Stderr,"Done - no\n");
	      PolynomialSet cg=*i;
	      
	      //buchberger(&cg,StandardGradedLexicographicTermOrder()); The cg is already a Groebner basis
	      
	      bool inList=false;
	      for(BergmanFan::MaximalConeList::const_iterator j=bfan.cones.begin();j!=bfan.cones.end();j++)
		{
		  if(areIdealsEqual(j->coneGroebnerBasis,cg))
		    {
		      inList=true;
		      break;
		    }
		}
	      if(!inList)
		{
		  bfan.cones.push_back(BergmanFan::MaximalCone(cg,*g,true,maximalConeLabel++));
		}
	    }
	  //	  else
	    // fprintf(Stderr,"Done - yes\n");
	  }
	else
	  fprintf(Stderr,"dimension check fails\n");
	}
    }
  //fprintf(Stderr,"No duplicates:\n");  
  
  /*  for(PolynomialSetList::const_iterator i=tropical.begin();i!=tropical.end();i++)
    {
      AsciiPrinter(Stdout).printPolynomialList(i->coneGroebnerBasis);
      int coDim=rankOfMatrix(wallInequalities(*i));
      int d=i->numberOfVariablesInRing()-coDim;
      fprintf(Stderr,"%i\n",d);
    }
  */
  return bfan;
}


/*static bool staticInOrbit(PolynomialSet const &groebnerBasis1, PolynomialSet const &groebnerBasis2, SymmetryGroup const &s)
{
  for(SymmetryGroup::ElementContainer::const_iterator j=s.elements.begin();j!=s.elements.end();j++)
    if(areIdealsEqual(SymmetryGroup::permutePolynomialSet(groebnerBasis2,*j),groebnerBasis1))return true;
  return false;
  }*/

/*static bool staticPermutationFixesCone(PolynomialSet const &groebnerBasis, IntegerVector const &v)
{ // Cone is fixed iff the cone ideal is fixed.
  PolynomialSet q(groebnerBasis.getRing());
  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    {
      q.push_back(SymmetryGroup::permutePolynomial(*i,v));
    }
  return areIdealsEqual(q,groebnerBasis);
  }*/


/*static int staticOrbitSize(PolynomialSet const &groebnerBasis, SymmetryGroup const &s)
{
  int groupSize=s.elements.size();
  
  int numFixed=0;
  for(SymmetryGroup::ElementContainer::const_iterator j=s.elements.begin();j!=s.elements.end();j++)
    if(staticPermutationFixesCone(groebnerBasis,*j))numFixed++;
  
  //  fprintf(Stderr,"groupSize = %i, numFixed = %i\n",groupSize,numFixed);
  return groupSize/numFixed;
}
*/



/** Objects of this class represent ConeOrbits of facet of a tropical variety.
 */

class ConeOrbit
{
public:
  const SymmetryGroup &s;
  PolynomialSet coneGroebnerBasis;
  PolynomialSet idealGroebnerBasis;
  int label;
  //  PolynomialSetList markedFacets;
  PolyhedralCone theCone;
  //  IntegerVectorList markedFacets; //marked facets of dim d cone. The facet of a facet is a ridge of the tropical variety.
  IntegerVector stableRay; //stable relative interior point
  IntegerVectorList unMarkedFacetNormals;

  ConeOrbit(const SymmetryGroup &s_, PolynomialSet const &coneGroebnerBasis_, PolynomialSet const &idealGroebnerBasis_, int label_):
    coneGroebnerBasis(coneGroebnerBasis_),
    idealGroebnerBasis(idealGroebnerBasis_),
    label(label_),
    s(s_),
    theCone(fastNormals(wallInequalities(idealGroebnerBasis_)),//wallFlipableNormals(idealGroebnerBasis_,false),
	    wallInequalities(coneGroebnerBasis_),
	    idealGroebnerBasis_.getRing().getNumberOfVariables())
  {
    /*
log0 AsciiPrinter(Stderr).printPolynomialSet(coneGroebnerBasis_);
log0 AsciiPrinter(Stderr).printPolynomialSet(idealGroebnerBasis_);

    log0 fprintf(Stderr,"inequalities:\n");
    log0 AsciiPrinter(Stderr).printVectorList(wallInequalities(coneGroebnerBasis_));
    log0 fprintf(Stderr,"equalities:\n");
    log0 AsciiPrinter(Stderr).printVectorList(wallFlipableNormals(idealGroebnerBasis_,false));

    AsciiPrinter P(Stderr);
    log0 theCone.print(&P);
    */
    theCone.findFacets();
    //    theCone.canonicalize();//REMOVE ME JAN 2009
    stableRay=PolyhedralFan::stableRay(theCone,&s);
    unMarkedFacetNormals=theCone.getHalfSpaces();
  }
  IntegerVector getStableRay()const
  {
    return stableRay;
  }
  bool doesPermutationFixCone(IntegerVector const &v)const
  { // Cone is fixed iff stableRay is fixed.
    return (SymmetryGroup::compose(v,stableRay)==stableRay);
  }
  int orbitSize()const
  {
    int groupSize=s.elements.size();
    
    int numFixed=0;
    for(SymmetryGroup::ElementContainer::const_iterator j=s.elements.begin();j!=s.elements.end();j++)
      if(doesPermutationFixCone(*j))numFixed++;
    
    //  fprintf(Stderr,"groupSize = %i, numFixed = %i\n",groupSize,numFixed);
    return groupSize/numFixed;
  }

  void markFacet(IntegerVector const &relIntRidgeVector)
  {
    for(SymmetryGroup::ElementContainer::const_iterator j=s.elements.begin();j!=s.elements.end();j++)
      if(doesPermutationFixCone(*j))
	{
	  IntegerVector relIntRidgeVector2=SymmetryGroup::compose(*j,relIntRidgeVector);
	  {
	    for(IntegerVectorList::iterator i=unMarkedFacetNormals.begin();i!=unMarkedFacetNormals.end();)
	      if(dotLong(*i,relIntRidgeVector2)==0)
		{
		  IntegerVectorList::iterator temp=i;
		  temp++;
		  unMarkedFacetNormals.erase(i);
		  i=temp;
		}
	      else
		i++;
	  }
	}
  }
  
  /*  bool containsAndMark(PolynomialSet const &coneGroebnerBasis_, PolynomialSet const &facetIdeal, IntegerVector *labelPermutation)
  {
    for(SymmetryGroup::ElementContainer::const_iterator j=s.elements.begin();j!=s.elements.end();j++)
      if(areIdealsEqual(coneGroebnerBasis,SymmetryGroup::permutePolynomialSet(coneGroebnerBasis_,*j)))
	{
	  PolynomialSet facetIdeal2=SymmetryGroup::permutePolynomialSet(facetIdeal,*j);
	  bool found=false;
	  for(PolynomialSetList::const_iterator i=markedFacets.begin();i!=markedFacets.end();i++)
	    if(areIdealsEqual(*i,facetIdeal2))
	      {
		found=true;
		break;
	      }
	  //	  assert(!found); //this is not a mistake is it?
	  markedFacets.push_back(facetIdeal2);
	  if(labelPermutation)*labelPermutation=*j;
	  return true;
	}
    return false;
    }*/
  /* Alternative using only geometric information.
   */
  bool containsAndMark(IntegerVector const &relIntFacetVector, IntegerVector const &relIntRidgeVector, IntegerVector *labelPermutation)
  {
    for(SymmetryGroup::ElementContainer::const_iterator j=s.elements.begin();j!=s.elements.end();j++)
      if(theCone.contains(SymmetryGroup::compose(*j,relIntFacetVector)))
	{
	  markFacet(SymmetryGroup::compose(*j,relIntRidgeVector));
	  if(labelPermutation)*labelPermutation=*j;
	  return true;
	}
    return false;
  }
  bool hasUnmarkedFacet()const
  {
    return(!unMarkedFacetNormals.empty());    
  }
  IntegerVector getUnmarkedFacet()const
  {
    assert(!unMarkedFacetNormals.empty());
    return *unMarkedFacetNormals.begin();
  }
  /*  bool isMarkedFacet(PolynomialSet const &f)
  {
    for(SymmetryGroup::ElementContainer::const_iterator j=s.elements.begin();j!=s.elements.end();j++)
      if(staticPermutationFixesCone(coneGroebnerBasis,*j))
	for(PolynomialSetList::const_iterator i=markedFacets.begin();i!=markedFacets.end();i++)
	  if(areIdealsEqual(SymmetryGroup::permutePolynomialSet(f,*j),*i))return true;
    return false;
    }*/
  void print(AsciiPrinter &p)const
  {
    p.printString("ConeOrbit{\n");
    p.printInteger(label);
    p.printNewLine();
    IntegerVector v=getStableRay();
    /*    if(v.size()==25)
      {
	p.printVector(v.subvector(0,5));p.printNewLine();
	p.printVector(v.subvector(5,10));p.printNewLine();
	p.printVector(v.subvector(10,15));p.printNewLine();
	p.printVector(v.subvector(15,20));p.printNewLine();
	p.printVector(v.subvector(20,25));p.printNewLine();	
	}*/
    p.printVector(v);
    p.printNewLine();
    p.printString("Orbit size:");
    p.printInteger(orbitSize());

    p.printString("\nConeIdeal:\n");
    p.printPolynomialSet(coneGroebnerBasis);
    p.printString("\nFullIdeal:\n");
    p.printPolynomialSet(idealGroebnerBasis);
    //    p.printString("Unmarked facet normals:\n");
    //    p.printVectorList(unMarkedFacetNormals);
    theCone.print(&p);
    p.printString("}ConeOrbit\n");
  }
  IntegerVectorList getEquations()const
  {
    return theCone.getEquations();
  }
  IntegerVectorList getInequalities()const
  {
    return theCone.getHalfSpaces();
  }
};

class ConeOrbitContainer
{
  typedef list<ConeOrbit> ConeOrbitList;
  ConeOrbitList l;
public:
  void push_back(const ConeOrbit &orbit)
  {
    l.push_back(orbit);
  }
  bool empty()
  {
    return l.empty();
  }
  ConeOrbit &front()
  {
    return *l.begin();
  }
  int size()
  {
    return l.size();
  }
  void pop_front()
  {
    l.pop_front();
  }
  void print(AsciiPrinter &p)
  {
    p.printString("OrbitList{\n");
    for(ConeOrbitList::const_iterator i=l.begin();i!=l.end();i++)
      {
	i->print(p);
	//	p.printPolynomialSet(i->coneGroebnerBasis);
	//p.printNewLine();
      }
    p.printString("}OrbitList\n");
  }
  bool containsAndMark(IntegerVector const &relIntFacetVector, IntegerVector const &relIntRidgeVector, int *label, IntegerVector *labelPermutation)
  //  bool containsAndMark(PolynomialSet const &coneGroebnerBasis, PolynomialSet const &facetIdeal, int *label, IntegerVector *labelPermutation)
  {
    //    fprintf(Stderr,"listlength:%i",l.size());
    //  int I=0;
    // for(ConeOrbitList::const_iterator i=l.begin();i!=l.end();i++)
    //  {
    //	fprintf(Stderr,"%i",I++);
    //  }
    for(ConeOrbitList::iterator i=l.begin();i!=l.end();i++)
      {
	/*	fprintf(Stderr,"Comparing:\n");
	AsciiPrinter(Stderr).printPolynomialSet(coneGroebnerBasis);
	AsciiPrinter(Stderr).printPolynomialSet(i->coneGroebnerBasis);
	*/
	//	if(i->containsAndMark(coneGroebnerBasis,facetIdeal,labelPermutation))
	if(i->containsAndMark(relIntFacetVector,relIntRidgeVector,labelPermutation))
	   //	if(areIdealsEqual(coneGroebnerBasis,i->coneGroebnerBasis))// this could be slow!
	  {
	    *label=i->label;
	    return true;
	  }
	//	if(i->coneGroebnerBasis==coneGroebnerBasis)return true;
      }
    /*    fprintf(Stderr,"________________NOT IN LIST:\n");
    AsciiPrinter(Stderr).printPolynomialSet(coneGroebnerBasis);
    */
    return false;
  }
};

static int staticOrbitSize(PolynomialSet const &coneGroebnerBasis, PolynomialSet const &idealGroebnerBasis, SymmetryGroup const &s)
{
  return ConeOrbit(s, coneGroebnerBasis, idealGroebnerBasis, -1).orbitSize();
}

BergmanFan bergman(PolynomialSet const &coneGroebnerBasis1, PolynomialSet const &idealGroebnerBasis1, SymmetryGroup const *symmetryGroup)
{
  PolynomialRing theRing=coneGroebnerBasis1.getRing();
  bool useFanIntersection=true;
  bool isSimplicial=true;

  assert(coneGroebnerBasis1.numberOfVariablesInRing()==idealGroebnerBasis1.numberOfVariablesInRing());

  int n=coneGroebnerBasis1.numberOfVariablesInRing();

  SymmetryGroup localSymmetryGroup(n);
  if(!symmetryGroup)symmetryGroup=&localSymmetryGroup;

  BergmanFan ret;
  ret.setSymmetryGroup(*symmetryGroup);

  ConeOrbitContainer active;

  int maximalConeLabel=0;

  {
    ConeOrbit newConeOrbit(*symmetryGroup,coneGroebnerBasis1,idealGroebnerBasis1,maximalConeLabel++);
    log1 fprintf(Stderr,"Adding orbit of size: %i\n",newConeOrbit.orbitSize());
    active.push_back(newConeOrbit);
  }
  while(!active.empty())
    {
      log1 fprintf(Stderr,"\n-------------------------------------\n");
      log1 fprintf(Stderr,"Size of active set: %i, Completed: %i\n",active.size(),ret.cones.size());
      log1 fprintf(Stderr,"-------------------------------------\n");
      AsciiPrinter p(Stderr);
      
      
      ConeOrbit &current=active.front();
      
      /*      AsciiPrinter p2(stdout);//REMOVE THESE THREE LINES
      current.print(p2);
      fflush(stdout);
      */

      //      IntegerVectorList equalities=wallInequalities(current.coneGroebnerBasis);

      //      IntegerVectorList normals=wallFlipableNormals(current.idealGroebnerBasis,false);
      /*
      {
	PolyhedralCone p(normals,equalities);
	p.findFacets();
	isSimplicial&=p.isSimplicial();
	}HUSK AT INDSÆTTE DETTE IGEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      */
      //      removeRedundantRows(&normals,&equalities,true);//IS THIS RIGHT?

      IntegerVectorList facePerp=current.getEquations();//equalities;
	   /*
      int numberOfNormals=normals.size();
      int numberOfEqualities=equalities.size();
      IntegerVectorList inequalities=equalities;
      inequalities.splice(inequalities.begin(),normals,normals.begin(),normals.end());

      IntegerVector equalitySet(inequalities.size());
      for(int i=0;i<numberOfEqualities;i++)
	equalitySet[i+numberOfNormals]=1;

      IntegerVectorList::const_iterator i=inequalities.begin();
      */
      int numberOfValidFacets=0;
      int numberOfAlreadyMarkedFacets=0;
      while(current.hasUnmarkedFacet())
	//for(int I=0;I<numberOfNormals;I++)
	{

	  //      log1 active.print(p);
	  IntegerVector theUnmarkedFacetNormal=current.getUnmarkedFacet();


	  IntegerVector relIntRidgeVector;
	  {
	    IntegerVectorList equations=current.getEquations();
	    IntegerVectorList inequalities=current.getInequalities();
	    equations.push_back(theUnmarkedFacetNormal);
	    PolyhedralCone C(inequalities,equations,n);
	    relIntRidgeVector=C.getRelativeInteriorPoint();
	  }
	  //log0 fprintf(Stderr,"RIDGE:");
	  //log0 p.printVector(relIntRidgeVector);

	  //	  equalitySet[I]=1;
	  //	  if(hasInteriorPoint(inequalities,true,&equalitySet))
	    {
	      // compute initial ideal
	      PolynomialSet initialIdeal(theRing);
	      int oldDim=rankOfMatrix(facePerp);

	      //facePerp.push_back(*i);
	      facePerp.push_back(theUnmarkedFacetNormal);
	      int newDim=rankOfMatrix(facePerp);
	      initialIdeal=wallIdeal(current.idealGroebnerBasis,Subspace(facePerp));
	      facePerp.pop_back();


	      if(oldDim!=newDim)
		{
		  numberOfValidFacets++;
		  //  if(!current.isMarkedFacet(initialIdeal))
		    {
		      BergmanFan b=useFanIntersection?bergmanRayIntersection(initialIdeal):bergmanRay(initialIdeal);
		      
		      ret.codimensionOneCones.push_back(BergmanFan::CodimensionOneCone(initialIdeal));

		      assert(b.cones.size());


		      for(BergmanFan::MaximalConeList::const_iterator i=b.cones.begin();i!=b.cones.end();i++)
			{
			  /* Unfortunately we have to make a lift to
			     compute the relative interior point. This
			     was not needed in earlier versions where
			     symmetry checks were algebraic and not
			     geometric*/
			  //  assert(i->idealGroebnerBasis.isMarked());
			  log2 cerr<<"8";

			  PolynomialSet g2(theRing);
			  WeightTermOrder termOrder(termorderWeight(initialIdeal));

			  for(PolynomialSet::const_iterator j=i->idealGroebnerBasis.begin();j!=i->idealGroebnerBasis.end();j++)
			    g2.push_back(divisionLift(*j, initialIdeal, current.idealGroebnerBasis, termOrder));
			  assert(g2.isMarked());
			  log2 cerr<<"9";
			
			 
			  if(1)
			    {
			      log2 cerr << "AUTOREDUCTION";
			      autoReduce(&g2,LexicographicTermOrder());
			      log2 cerr << "AUTOREDUCTION END"<<endl;
			    }
			  else
			    {
			      log0 cerr << "RESTRICTED AUTOREDUCTION";
			      autoReduceRestricted(&g2,i->idealGroebnerBasis);
			      log0 cerr << "RESTRICTED AUTOREDUCTION END"<<endl;
			    }
			  ConeOrbit tempConeOrbit(*symmetryGroup,i->coneGroebnerBasis,g2,-1);
			  log2 cerr<<"0";
			  IntegerVector relIntFacetVector=tempConeOrbit.getStableRay();

			  log2 cerr<<"A";

			  int label=-1;
			  IntegerVector labelPermutation;
			  if(!active.containsAndMark(relIntFacetVector,relIntRidgeVector,&label,&labelPermutation))
			    //			  if(!active.containsAndMark(i->coneGroebnerBasis,initialIdeal,&label,&labelPermutation))
			    {
			  log2 cerr<<"B";
			      if(!ret.contains(i->coneGroebnerBasis))
				{
				  label=maximalConeLabel++;
				  labelPermutation=SymmetryGroup::identity(n);
				  ConeOrbit newConeOrbit(*symmetryGroup,i->coneGroebnerBasis,g2,label);
				  log1 fprintf(Stderr,"Adding orbit of size: %i\n",newConeOrbit.orbitSize());
				  //newConeOrbit.markFacet(initialIdeal);
				  newConeOrbit.markFacet(relIntRidgeVector);
				  active.push_back(newConeOrbit);
				}
			      else
				{
				  assert(0);//Can this ever happen? REMOVE JAN 2009
				}
			  log2 cerr<<"C";
			    }
			  log2 cerr<<"D";
			  ret.codimensionOneCones.back().incidenceList.push_back(label);
			  log2 cerr<<"E";
			  ret.codimensionOneCones.back().incidencePermutationList.push_back(labelPermutation);
			  log2 cerr<<"F";
			}
		    }
		    //		  else
		    //  numberOfAlreadyMarkedFacets++;
		}

	      else
		{
		  assert(0);
		}






	    }
	    //	  equalitySet[I]=0;
	    //  i++;

	    //assert(0);//for timing purposes
	}
      log1 fprintf(Stderr,"Done processing this orbit - Number of valid facets: %i   Number of already marked facets: %i\n",numberOfValidFacets,numberOfAlreadyMarkedFacets);
      ret.cones.push_back(BergmanFan::MaximalCone(current.coneGroebnerBasis,current.idealGroebnerBasis,false,current.label,numberOfValidFacets));
      active.pop_front();
    }
  ret.setSimplicial(isSimplicial);
  return ret;
}


//--------------------------------------
// BergmanFan
//--------------------------------------

BergmanFan::MaximalCone::MaximalCone(PolynomialSet const &coneGroebnerBasis_, PolynomialSet const &idealGroebnerBasis_, bool storeIdealBasis, int label_, int numberOfFacets_):
  idealBasisStored(storeIdealBasis),
  coneGroebnerBasis(coneGroebnerBasis_),
  idealGroebnerBasis(idealGroebnerBasis_),
  label(label_),
  numberOfFacets(numberOfFacets_),
  multiplicity(-1),
  theCone(fastNormals(wallInequalities(idealGroebnerBasis_)), 
	  wallInequalities(coneGroebnerBasis_),
	  idealGroebnerBasis_.getRing().getNumberOfVariables())
{
  if(!storeIdealBasis)idealGroebnerBasis=PolynomialSet(coneGroebnerBasis_.getRing());
  theCone.canonicalize();
}

int BergmanFan::numberOfMaximalCones()const
{
  int ret=0;;
  for(MaximalConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      IntegerVector stableRay=PolyhedralFan::stableRay(i->theCone,&symmetryGroup);
      //ConeOrbit C(symmetryGroup, i->coneGroebnerBasis, i->idealGroebnerBasis, -1);
      ret+=symmetryGroup.orbitSize(stableRay);//C.orbitSize();
	//	ret+=staticOrbitSize(i->coneGroebnerBasis, symmetryGroup);
    }
  return ret;
}

void BergmanFan::print(Printer &p)
{
  int numberOfMaximalCones=0;
  p.printString("Printing tropical variety modulo symmetry\n");
  p.printString("-----------------");
  p.printNewLine();
  p.printString("1. Maximal cones:\n");
  p.printString("-----------------");
  p.printNewLine();

  for(MaximalConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      p.printString("Orbit index: ");
      p.printInteger(i->label);
      p.printString("\n");
      p.printString("Initial ideal:\n");
      p.printPolynomialSet(i->coneGroebnerBasis);
      //p.printPolynomialSet(i->idealGroebnerBasis);
      int orbitSize=(symmetryGroup.sizeOfBaseSet())?staticOrbitSize(i->coneGroebnerBasis,i->coneGroebnerBasis,symmetryGroup):1;
      p.printString("OrbitSize:");
      p.printInteger(orbitSize);
      p.printString("\n");
      p.printString("NumberOfFacets:");
      p.printInteger(i->numberOfFacets);
      p.printString("\n\n");
      numberOfMaximalCones+=orbitSize;

      /*      {
	IntegerVectorList normals=wallInequalities(i->idealGroebnerBasis);
	IntegerVectorList equations=wallInequalities(i->coneGroebnerBasis);
	PolyhedralCone c(normals,equations);
	c.canonicalize();
	c.print(&p);
      }
      */

      if(symmetryGroup.sizeOfBaseSet())
      {
	list<int> indices;
	int index=0;


	/*  //This needs to be changed so that torus permutations are considered too.
	for(MaximalConeList::const_iterator j=cones.begin();j!=cones.end();j++,index++)
	  if(staticInOrbit(j->coneGroebnerBasis,i->coneGroebnerBasis,symmetryGroup))indices.push_back(index);
	*/

	if(indices.size()>1)
	  {
	    fprintf(Stderr,"Conflicting orbits!!!!:");
	    for(list<int>::const_iterator j=indices.begin();j!=indices.end();j++)
	      fprintf(Stderr," %i ",*j);
	    fprintf(Stderr,"\n");
	  }
      }
    }
  p.printString("-----------------------");
  p.printNewLine();
  p.printString("2. Codimension 1 cones:\n");
  p.printString("-----------------------");
  p.printNewLine();
  for(CodimensionOneConeList::const_iterator i=codimensionOneCones.begin();i!=codimensionOneCones.end();i++)
    {
      //p.printString("----------------------");
      //      p.printNewLine();
      p.printString("Groebner basis of initial ideal:\n");
      p.printPolynomialSet(i->idealGroebnerBasis);
      //      p.printNewLine();
      p.printString("Adjacent maximal cones (orbit index, permutation):\n");
      p.printString("(");
      IntegerVectorList::const_iterator J=i->incidencePermutationList.begin();
      for(list<int>::const_iterator j=i->incidenceList.begin();j!=i->incidenceList.end() && J!=i->incidencePermutationList.end();j++,J++)
	{
	  if(j!=i->incidenceList.begin())p.printString(",\n");
	  p.printString("(");
	  p.printInteger(*j);
	  p.printString(", ");
	  p.printVector(*J);
	  p.printString(")");
	}
      p.printString(")\n");
      {
	int index=*(i->incidenceList.begin());
	MaximalConeList::const_iterator j=cones.begin();
	while(index>0){index--;j++;}
	//	IntegerVectorList normals=wallInequalities(SymmetryGroup::permutePolynomialSet(j->idealGroebnerBasis,*(i->incidencePermutationList.begin())));

	/*
	  IntegerVectorList normals=wallInequalities(SymmetryGroup::permutePolynomialSet(j->idealGroebnerBasis,SymmetryGroup::inverse(*(i->incidencePermutationList.begin()))));
	IntegerVectorList equations=wallInequalities(i->idealGroebnerBasis);
	PolyhedralCone c(normals,equations);
	c.canonicalize();
	c.print(&p);
	*/

	/*	p.printString("-----------");
	p.printNewLine();
	PolyhedralFan F=PolyhedralFan::facetsOfCone(c);
	F.print(&p);
	p.printString("-----------");
	p.printNewLine();
	*/
	//p.printString("----------------------");
	p.printNewLine();
      }
    }
  p.printString("Done printing tropical variety - #maxcones=");
  p.printInteger(numberOfMaximalCones);
  p.printString(" (");
  p.printInteger(cones.size());
  p.printString(")  #codim1cones= ? (");
  p.printInteger(codimensionOneCones.size());
  p.printString(")");
  p.printNewLine();
}


bool BergmanFan::contains(PolynomialSet const &g)
{
  for(MaximalConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      if(areIdealsEqual(g,i->coneGroebnerBasis))
	{

	  return true;
	}
    }
  return false;
}


void BergmanFan::setSymmetryGroup(SymmetryGroup const &s)
{
  symmetryGroup=s;
}


PolyhedralFan BergmanFan::toPolyhedralFan()const
{
  assert(!cones.empty());
  int n=cones.begin()->idealGroebnerBasis.numberOfVariablesInRing();
  PolyhedralFan ret(n);

  for(MaximalConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      //PolyhedralCone c1(wallInequalities(i->idealGroebnerBasis),wallInequalities(i->coneGroebnerBasis));
      //      PolyhedralCone c1(wallFlipableNormals(i->idealGroebnerBasis,false),wallInequalities(i->coneGroebnerBasis),n);
      PolyhedralCone c1=i->theCone;
      log2 fprintf(Stderr,"Cononicalising...\n");
      c1.canonicalize();
      log2 fprintf(Stderr,"... done canonicalising...\n");
      //      fprintf(Stderr,"a\n");
      //      for(SymmetryGroup::ElementContainer::const_iterator j=symmetryGroup.elements.begin();j!=symmetryGroup.elements.end();j++)
      //	{
	  /*	  IntegerVectorList normals=wallInequalities(SymmetryGroup::permutePolynomialSet(i->idealGroebnerBasis,*j));
	  IntegerVectorList equations=wallInequalities(SymmetryGroup::permutePolynomialSet(i->coneGroebnerBasis,*j));
	  PolyhedralCone c(normals,equations);
	  */
      //	  PolyhedralCone c=c1.permuted(*j);
      //	  c.canonicalize();
      c1.setMultiplicity(i->multiplicity);
	  ret.insert(c1);
	  //	}
    }
  return ret;
}


void BergmanFan::setSimplicial(bool b)
{
  simplicial=b;
}


bool BergmanFan::isSimplicial()const
{
  return simplicial;
}


void BergmanFan::computeMultiplicities()
{
  for(MaximalConeList::iterator i=cones.begin();i!=cones.end();i++)
    i->multiplicity=multiplicity(i->coneGroebnerBasis);
      //MULTIPLICITY TEST
  //      AsciiPrinter(Stderr).printPolynomialSet(current.coneGroebnerBasis);
  //      fprintf(Stderr,"MULTIPLICITY :%i\n",multiplicity(current.coneGroebnerBasis));
}
