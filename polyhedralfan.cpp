#include <sstream>
#include "polyhedralfan.h"
#include "reversesearch.h"
#include "wallideal.h"
#include "buchberger.h"
#include "printer.h"
#include "timer.h"
#include "symmetry.h"
#include "polymakefile.h"
#include "symmetriccomplex.h"
#include "linalg.h"
#include "lp.h"
#include "codimoneconnectedness.h"
#include "symmetrictraversal.h"
#include "traverser_groebnerfan.h"
#include "log.h"

static Timer polyhedralFanRefinementTimer("Polyhedral fan refinement",1);

PolyhedralFan::PolyhedralFan(int ambientDimension):
  n(ambientDimension)
{
}


PolyhedralFan PolyhedralFan::fullSpace(int n)
{
  PolyhedralFan ret(n);

  PolyhedralCone temp(n);
  temp.canonicalize();
  ret.cones.insert(temp);

  return ret;
}


PolyhedralFan PolyhedralFan::halfSpace(int n, int i)
{
  assert(0<=i);
  assert(i<n);
  PolyhedralFan ret(n);

  IntegerVector v(n);
  v[i]=-1;
  IntegerVectorList l;
  IntegerVectorList empty;
  l.push_back(v);
  ret.cones.insert(PolyhedralCone(l,empty,n));

  return ret;
}


PolyhedralFan PolyhedralFan::facetsOfCone(PolyhedralCone const &c)
{
  PolyhedralCone C(c);
  C.canonicalize();
  PolyhedralFan ret(C.ambientDimension());

  IntegerVectorList halfSpaces=C.getHalfSpaces();

  for(IntegerVectorList::const_iterator i=halfSpaces.begin();i!=halfSpaces.end();i++)
    {
      IntegerVectorList a;
      IntegerVectorList b;
      b.push_back(*i);
      PolyhedralCone n=intersection(PolyhedralCone(a,b),c);
      n.canonicalize();
      ret.cones.insert(n);
    }
  return ret;
}

PolyhedralFan PolyhedralFan::complementOfCone(PolyhedralCone const &c, bool includec)
{
  PolyhedralCone C=c;
  C.canonicalize();
  IntegerVectorList inequalities=C.getHalfSpaces();
  IntegerVectorList equations=C.getEquations();

  for(IntegerVectorList::const_iterator i=equations.begin();i!=equations.end();i++)
    inequalities.push_back(*i);

  int n=C.ambientDimension();

  PolyhedralFan ret=PolyhedralFan::fullSpace(n);
  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
    {
      PolyhedralFan temp(C.ambientDimension());
      for(int j=-1;j<2;j+=2)
	{
	  IntegerVectorList inequality;
	  inequality.push_back(j* *i);
	  IntegerVectorList empty;
	  PolyhedralCone tempC(inequality,empty,n);
	  tempC.canonicalize();
	  temp.insert(tempC);
	}
      ret=refinement(ret,temp);
    }
  if(!includec)ret.remove(C);
  return ret;
}

PolyhedralFan PolyhedralFan::bergmanOfPrincipalIdeal(Polynomial const &p1)
{
  PolynomialRing theRing=p1.getRing();
  PolynomialRing theSecondRing=theRing.withVariablesAppended("H");
  Polynomial p=p1.homogenization(theSecondRing);
  PolyhedralFan ret(p1.getNumberOfVariables());
  PolynomialSet g(theSecondRing);
  g.push_back(p);
  buchberger(&g,LexicographicTermOrder());

  EnumerationTargetCollector gfan;
  LexicographicTermOrder myTermOrder;
  ReverseSearch rs(myTermOrder);

  rs.setEnumerationTarget(&gfan);

  rs.enumerate(g);

  PolynomialSetList theList=gfan.getList();
  for(PolynomialSetList::const_iterator i=theList.begin();i!=theList.end();i++)
    {
      //      AsciiPrinter(Stderr).printPolynomialSet(*i);
      IntegerVectorList inequalities=wallInequalities(*i);
      IntegerVectorList f=wallFlipableNormals(*i,true);
      for(IntegerVectorList::const_iterator j=f.begin();j!=f.end();j++)
	{
	  //	  AsciiPrinter(Stderr).printVector(*j);
	  if(myTermOrder(*j,*j-*j))
	    {
	      IntegerVectorList equalities;
	      equalities.push_back(*j);
	      PolyhedralCone c=PolyhedralCone(inequalities,equalities).withLastCoordinateRemoved();
	      c.canonicalize();
	      c.setLinearForm(i->begin()->getMarked().m.exponent.subvector(0,p1.getNumberOfVariables()));
	      c.setMultiplicity(gcdOfVector(j->subvector(0,j->size()-1)));
	      ret.cones.insert(c);
	    }
	}
    }

  return ret;
}


PolyhedralFan PolyhedralFan::normalFanOfNewtonPolytope(Polynomial const &p1)
{
  PolynomialRing theRing=p1.getRing();
  PolynomialRing theSecondRing=theRing.withVariablesAppended("H");
  Polynomial p=p1.homogenization(theSecondRing);
  PolyhedralFan ret(p1.getNumberOfVariables());
  PolynomialSet g(theSecondRing);
  //  PolynomialSet g(theRing);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  g.push_back(p);
  buchberger(&g,LexicographicTermOrder());

  EnumerationTargetCollector gfan;
/*  {//old enumeration strategy
  LexicographicTermOrder myTermOrder;
  ReverseSearch rs(myTermOrder);

  rs.setEnumerationTarget(&gfan);

  rs.enumerate(g);
  }
*/
  {//new enumeration strategy
    GroebnerFanTraverser traverser(g);
    traverser.setIsKnownToBeComplete(true);
    TargetGlue target(gfan);
    symmetricTraverse(traverser,target);
  }

  PolynomialSetList theList=gfan.getList();
  for(PolynomialSetList::const_iterator i=theList.begin();i!=theList.end();i++)
    {
      IntegerVectorList inequalities=wallInequalities(*i);
      IntegerVectorList equalities;

      PolyhedralCone c=PolyhedralCone(inequalities,equalities).withLastCoordinateRemoved();
      c.canonicalize();
      c.setLinearForm(i->begin()->getMarked().m.exponent.subvector(0,c.ambientDimension()));
      ret.cones.insert(c);
    }

  return ret;
}


void PolyhedralFan::print(class Printer *p)const
{
  p->printString("Printing PolyhedralFan");
  p->printNewLine();
  p->printString("Ambient dimension: ");
  p->printInteger(n);
  p->printNewLine();
  p->printString("Number of cones: ");
  p->printInteger(cones.size());
  p->printNewLine();
  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      p->printNewLine();
      p->printPolyhedralCone(*i);
    }
  p->printString("Done printing PolyhedralFan.");
  p->printNewLine();
}

int PolyhedralFan::getAmbientDimension()const
{
  return n;
}

bool PolyhedralFan::isEmpty()const
{
  return cones.empty();
}

int PolyhedralFan::getMaxDimension()const
{
  assert(!cones.empty());

  return cones.begin()->dimension();
}

int PolyhedralFan::getMinDimension()const
{
  assert(!cones.empty());

  return cones.rbegin()->dimension();
}

PolyhedralFan refinement(const PolyhedralFan &a, const PolyhedralFan &b, int cutOffDimension, bool allowASingleConeOfCutOffDimension)
{
  TimerScope ts(&polyhedralFanRefinementTimer);
  //  fprintf(Stderr,"PolyhedralFan refinement: #A=%i #B=%i\n",a.cones.size(),b.cones.size());
  int conesSkipped=0;
  int numberOfComputedCones=0;
  bool linealityConeFound=!allowASingleConeOfCutOffDimension;
  assert(a.getAmbientDimension()==b.getAmbientDimension());

  PolyhedralFan ret(a.getAmbientDimension());

  for(PolyhedralConeList::const_iterator A=a.cones.begin();A!=a.cones.end();A++)
    for(PolyhedralConeList::const_iterator B=b.cones.begin();B!=b.cones.end();B++)
      {
	PolyhedralCone c=intersection(*A,*B);
	int cdim=c.dimension();
	//	assert(cdim>=linealitySpaceDimension);
	bool thisIsLinealityCone=(cutOffDimension>=cdim);
	if((!thisIsLinealityCone)||(!linealityConeFound))
	  {
	    numberOfComputedCones++;
	    c.canonicalize();
	    ret.cones.insert(c);
	    linealityConeFound=linealityConeFound || thisIsLinealityCone;
	  }
	else
	  {
	    conesSkipped++;
	  }
      }
  //  fprintf(Stderr,"Number of skipped cones: %i, lineality cone found: %i\n",conesSkipped,linealityConeFound);
  //  fprintf(Stderr,"Number of computed cones: %i, number of unique cones: %i\n",numberOfComputedCones,ret.cones.size());

  return ret;
}


PolyhedralFan product(const PolyhedralFan &a, const PolyhedralFan &b)
{
  PolyhedralFan ret(a.getAmbientDimension()+b.getAmbientDimension());

  for(PolyhedralConeList::const_iterator A=a.cones.begin();A!=a.cones.end();A++)
    for(PolyhedralConeList::const_iterator B=b.cones.begin();B!=b.cones.end();B++)
      ret.insert(product(*A,*B));

  return ret;
}


IntegerVectorList PolyhedralFan::getRays(int dim)
{
  IntegerVectorList ret;
  for(PolyhedralConeList::iterator i=cones.begin();i!=cones.end();i++)
    {
      if(i->dimension()==dim)
	ret.push_back(i->getRelativeInteriorPoint());
    }
  return ret;
}


IntegerVectorList PolyhedralFan::getRelativeInteriorPoints()
{
  IntegerVectorList ret;
  for(PolyhedralConeList::iterator i=cones.begin();i!=cones.end();i++)
    {
      ret.push_back(i->getRelativeInteriorPoint());
    }
  return ret;
}


PolyhedralCone const& PolyhedralFan::highestDimensionalCone()const
{
  return *cones.rbegin();
}

void PolyhedralFan::insert(PolyhedralCone const &c)
{
  cones.insert(c);
}


void PolyhedralFan::insertFacetsOfCone(PolyhedralCone const &c)
{
  PolyhedralFan facets=facetsOfCone(c);
  for(PolyhedralConeList::const_iterator i=facets.cones.begin();i!=facets.cones.end();i++)insert(*i);
}


void PolyhedralFan::remove(PolyhedralCone const &c)
{
  cones.erase(c);
}

void PolyhedralFan::removeAllExcept(int a)
{
  PolyhedralConeList::iterator i=cones.begin();
  while(a>0)
    {
      if(i==cones.end())return;
      i++;
      a--;
    }
  cones.erase(i,cones.end());
}

void PolyhedralFan::removeAllLowerDimensional()
{
  if(!cones.empty())
    {
      int d=getMaxDimension();
      PolyhedralConeList::iterator i=cones.begin();
      while(i!=cones.end() && i->dimension()==d)i++;
      cones.erase(i,cones.end());
    }
}


PolyhedralFan PolyhedralFan::facetComplex()const
{
  //  fprintf(Stderr,"Computing facet complex...\n");
  PolyhedralFan ret(n);

  for(PolyhedralConeList::iterator i=cones.begin();i!=cones.end();i++)
    {
      PolyhedralFan a=facetsOfCone(*i);
      for(PolyhedralConeList::const_iterator j=a.cones.begin();j!=a.cones.end();j++)
	ret.insert(*j);
    }
  //  fprintf(Stderr,"Done computing facet complex.\n");
  return ret;
}


PolyhedralFan PolyhedralFan::fullComplex()const
{
  PolyhedralFan ret=*this;

  while(1)
    {
      log2 debug<<"looping";
      bool doLoop=false;
      PolyhedralFan facets=ret.facetComplex();
      log2 debug<<"number of facets"<<facets.size()<<"\n";
      for(PolyhedralConeList::const_iterator i=facets.cones.begin();i!=facets.cones.end();i++)
	if(!ret.contains(*i))
	  {
	    ret.insert(*i);
	    doLoop=true;
	  }
      if(!doLoop)break;
    }
  return ret;
}


/*
PolyhedralFan PolyhedralFan::facetComplexSymmetry(SymmetryGroup const &sym, bool keepRays, bool dropLinealitySpace)const
{
  log1 fprintf(Stderr,"Computing facet complex...\n");
  PolyhedralFan ret(n);

  if(keepRays)
    for(PolyhedralConeList::iterator i=cones.begin();i!=cones.end();i++)
      if(i->dimension()==i->dimensionOfLinealitySpace()+1)ret.insert(*i);

  for(PolyhedralConeList::iterator i=cones.begin();i!=cones.end();i++)
    {
      PolyhedralFan a=facetsOfCone(*i);
      for(PolyhedralConeList::const_iterator j=a.cones.begin();j!=a.cones.end();j++)
	{
	  if((!dropLinealitySpace) || j->dimension()!=j->dimensionOfLinealitySpace())
	    {
	      IntegerVector v=j->getRelativeInteriorPoint();
	      bool alreadyInRet=false;
	      for(SymmetryGroup::ElementContainer::const_iterator k=sym.elements.begin();k!=sym.elements.end();k++)
		{
		  IntegerVector u=SymmetryGroup::compose(*k,v);
		  if(!j->containsRelatively(u))
		    {
		      for(PolyhedralConeList::const_iterator l=ret.cones.begin();l!=ret.cones.end();l++)
			if(l->containsRelatively(u))alreadyInRet=true;
		    }
		}
	      if(!alreadyInRet)ret.insert(*j);
	    }
	}
    }
  log1 fprintf(Stderr,"Done computing facet complex.\n");
  return ret;
}
*/

PolyhedralFan PolyhedralFan::facetComplexSymmetry(SymmetryGroup const &sym, bool keepRays, bool dropLinealitySpace)const
{
  log1 fprintf(Stderr,"Computing facet complex...\n");
  PolyhedralFan ret(n);

  vector<IntegerVector> relIntTable;
  vector<int> dimensionTable;

  if(keepRays)
    for(PolyhedralConeList::iterator i=cones.begin();i!=cones.end();i++)
      if(i->dimension()==i->dimensionOfLinealitySpace()+1)
	{
	  relIntTable.push_back(i->getRelativeInteriorPoint());
	  dimensionTable.push_back(i->dimension());
	  ret.insert(*i);
	}

  for(PolyhedralConeList::iterator i=cones.begin();i!=cones.end();i++)
    {
      int iDim=i->dimension();
      if(dropLinealitySpace && (i->dimension()==i->dimensionOfLinealitySpace()+1))continue;

      //      i->findFacets();
      IntegerVectorList normals=i->getHalfSpaces();
      for(IntegerVectorList::const_iterator j=normals.begin();j!=normals.end();j++)
	{
	  bool alreadyInRet=false;
	  for(int l=0;l<relIntTable.size();l++)
	    {
	      if(dimensionTable[l]==iDim-1)
		for(SymmetryGroup::ElementContainer::const_iterator k=sym.elements.begin();k!=sym.elements.end();k++)
		  {
		    IntegerVector u=SymmetryGroup::compose(*k,relIntTable[l]);
		    if((dotLong(*j,u)==0) && (i->contains(u)))
		      {
			alreadyInRet=true;
			goto exitLoop;
		      }
		  }
	    }
	exitLoop:
	  if(!alreadyInRet)
	    {
	      IntegerVectorList equations=i->getEquations();
	      IntegerVectorList inequalities=i->getHalfSpaces();
	      equations.push_back(*j);
	      PolyhedralCone c(inequalities,equations,n);
	      c.canonicalize();
	      ret.insert(c);
	      relIntTable.push_back(c.getRelativeInteriorPoint());
	      dimensionTable.push_back(c.dimension());
	    }
	}
    }
  log1 fprintf(Stderr,"Done computing facet complex.\n");
  return ret;
}

/*
IntegerVectorList PolyhedralFan::getRaysInPrintingOrder(SymmetryGroup *sym)const
{
  assert(!cones.empty());
  int h=cones.begin()->dimensionOfLargestContainedSubspace();
  PolyhedralFan f=*this;
  while(f.getMaxDimension()!=h+1)
    {
      f=f.facetComplex();
    }

  IntegerVectorList rays;

  PolyhedralFan done(n);
  for(PolyhedralConeList::const_iterator i=f.cones.begin();i!=f.cones.end();i++)
    if(!done.contains(*i))
      for(SymmetryGroup::ElementContainer::const_iterator k=sym->elements.begin();k!=sym->elements.end();k++)
	{
	  PolyhedralCone cone=i->permuted(*k);
	  if(!done.contains(cone))
	    {
	      rays.push_back(cone.getRelativeInteriorPoint());
	      done.insert(cone);
	    }
	}
  return rays;
}
*/


IntegerVector PolyhedralFan::stableRay(PolyhedralCone const &c, SymmetryGroup const *sym)
{
  PolyhedralCone C=c;//cast away const instead?

  IntegerVector v=C.getRelativeInteriorPoint();

  IntegerVector ret(v.size());
  for(SymmetryGroup::ElementContainer::const_iterator k=sym->elements.begin();k!=sym->elements.end();k++)
    {
      IntegerVector v2=SymmetryGroup::compose(*k,v);
      if(c.contains(v2))ret+=v2;
    }
  return normalized(ret);
}


IntegerVector PolyhedralFan::stableRay(IntegerVector const &v, IntegerVectorList const &equations, IntegerVectorList const &inequalities, SymmetryGroup const *sym)
{
  IntegerVector ret(v.size());
  for(SymmetryGroup::ElementContainer::const_iterator k=sym->elements.begin();k!=sym->elements.end();k++)
    {
      IntegerVector v2=SymmetryGroup::compose(*k,v);
      bool containsV2=true;

      for(IntegerVectorList::const_iterator l=equations.begin();l!=equations.end();l++)
	if(dotLong(*l,v2)!=0)
	  {
	    containsV2=false;
	    goto leave;
	  }
      for(IntegerVectorList::const_iterator l=inequalities.begin();l!=inequalities.end();l++)
	if(dotLong(*l,v2)<0)
	  {
	    containsV2=false;
	    goto leave;
	  }
    leave:
      if(containsV2)ret+=v2;
    }
  return normalized(ret);
}

/* Slow version using facetComplexSymmetry()
IntegerVectorList PolyhedralFan::getRaysInPrintingOrder(SymmetryGroup *sym)const
{
  assert(!cones.empty());
  int h=cones.begin()->dimensionOfLinealitySpace();
  PolyhedralFan f=*this;
  if(f.getMaxDimension()==h)return IntegerVectorList();
  while(f.getMaxDimension()>h+1)
    {
      f=f.facetComplexSymmetry(*sym,true,true);
    }
  IntegerVectorList rays;

  for(PolyhedralConeList::const_iterator i=f.cones.begin();i!=f.cones.end();i++)
    {
      if(i->dimension()!=i->dimensionOfLinealitySpace())//This check is needed since the above while loop may not be run and therefore the lineality space may not have been removed.
	{
	  bool found=false;
	  for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++)
	    if(i->contains(*j))found=true;

	  if(!found)
	    {
	      //IntegerVector interiorPointForOrbit=i->getRelativeInteriorPoint();
	      IntegerVector interiorPointForOrbit=stableRay(*i,sym);
	      PolyhedralFan done(n);
	      for(SymmetryGroup::ElementContainer::const_iterator k=sym->elements.begin();k!=sym->elements.end();k++)
		{
		  PolyhedralCone cone=i->permuted(*k);

		  if(!done.contains(cone))
		    {
		      rays.push_back(SymmetryGroup::compose(*k,interiorPointForOrbit));
		      done.insert(cone);
		    }
		}
	    }
	}
    }
  return rays;
  }*/

PolyhedralFan PolyhedralFan::rayComplexSymmetry(SymmetryGroup const &sym)const
{
  //  log0 fprintf(Stderr,"rayComplexSymmetry - begin\n");
  PolyhedralFan ret(n);
  log1 fprintf(Stderr,"Computing rays of %i cones\n",cones.size());
  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      {
	static int t;
	if(!((t++)%10))log1 fprintf(Stderr,"%i\n",t);
      }
      //  log0 fprintf(Stderr,"calling\n");
      IntegerVectorList rays=i->extremeRays();
      //log0 fprintf(Stderr,"returning\n");
      for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++)
	{
	  bool alreadyInRet=false;
	  for(PolyhedralConeList::const_iterator I=ret.cones.begin();I!=ret.cones.end();I++)
	    for(SymmetryGroup::ElementContainer::const_iterator k=sym.elements.begin();k!=sym.elements.end();k++)
	      {
		IntegerVector u=SymmetryGroup::compose(*k,*j);
		if((I->contains(u)))
		  {
		    alreadyInRet=true;
		    goto exitLoop;
		  }
	      }
	exitLoop:
	  IntegerVectorList equations=i->getEquations();
	  IntegerVectorList inequalities1=i->getHalfSpaces();
	  IntegerVectorList inequalities2;
	  for(IntegerVectorList::const_iterator k=inequalities1.begin();k!=inequalities1.end();k++)
	    {
	      if(dotLong(*j,*k))
		inequalities2.push_back(*k);
	      else
		equations.push_back(*k);
	    }
	  PolyhedralCone C(inequalities2,equations,n);
	  C.canonicalize();
	  ret.insert(C);
	}
    }
  //  log0 fprintf(Stderr,"rayComplexSymmetry - end\n");
  return ret;
}


#if 0
IntegerVectorList PolyhedralFan::getRaysInPrintingOrder(SymmetryGroup const *sym)const
{
  assert(!cones.empty());
  int h=cones.begin()->dimensionOfLinealitySpace();

  /*
  PolyhedralFan f=*this;
  if(f.getMaxDimension()==h)return IntegerVectorList();
  while(f.getMaxDimension()>h+1)
    {
      f=f.facetComplexSymmetry(*sym,true,true);
    }
  */
  PolyhedralFan f=rayComplexSymmetry(*sym);
  IntegerVectorList rays;

  log1 fprintf(Stderr,"Number of cones in RayComplex: %i\n",f.cones.size());

  for(PolyhedralConeList::const_iterator i=f.cones.begin();i!=f.cones.end();i++)
    {
      static int t;
      log1 fprintf(Stderr,"%i\n",t++);
      if(i->dimension()!=i->dimensionOfLinealitySpace())//This check is needed since the above while loop may not be run and therefore the lineality space may not have been removed.
	{
	  bool found=false;
	  for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++)
	    if(i->contains(*j))
	      {
		found=true;
		break;
	      }
	  if(!found)
	    {
	      //IntegerVector interiorPointForOrbit=i->getRelativeInteriorPoint();
	      IntegerVector interiorPointForOrbit=stableRay(*i,sym);
	      //    PolyhedralFan done(n);

	      //Check that this works:
	      set<IntegerVector> thisOrbitsRays;

	      for(SymmetryGroup::ElementContainer::const_iterator k=sym->elements.begin();k!=sym->elements.end();k++)
		{
		  IntegerVector temp=SymmetryGroup::compose(*k,interiorPointForOrbit);
		  thisOrbitsRays.insert(temp);
		}
	      for(set<IntegerVector>::const_iterator i=thisOrbitsRays.begin();i!=thisOrbitsRays.end();i++)rays.push_back(*i);
	      //Instead of this:
	      /*	      for(SymmetryGroup::ElementContainer::const_iterator k=sym->elements.begin();k!=sym->elements.end();k++)
		{
		  bool found=false;
		  IntegerVector temp=SymmetryGroup::compose(*k,interiorPointForOrbit);
		  for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++)//REWRITE WITH LOGARITHMIC SEARCH
		    if(*j==temp)
		      {
			found=true;
			break;
		      }
		  if(!found)
		    {
		      PolyhedralCone cone=i->permuted(*k);
		      rays.push_back(SymmetryGroup::compose(*k,interiorPointForOrbit));
		      //      done.insert(cone);
		    }
		}
	      */
	    }
	}
    }
  return rays;
}

#elseif 0
//version used until Sep 2010
IntegerVectorList PolyhedralFan::getRaysInPrintingOrder(SymmetryGroup const *sym, bool upToSymmetry)const
{
	SymmetryGroup localsym(n);
	if(!sym)sym=&localsym;
  IntegerVectorList rays;
  log1 fprintf(Stderr,"Computing rays of %i cones\n",cones.size());
  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      {
	static int t;
	if(!((t++)%10))log1 fprintf(Stderr,"%i\n",t);
      }
      IntegerVectorList temp=i->extremeRays();
      //      AsciiPrinter(Stderr).printVectorList(temp);

      for(IntegerVectorList::const_iterator j=temp.begin();j!=temp.end();j++)
	{
	  IntegerVectorList equations=i->getEquations();
	  IntegerVectorList inequalities1=i->getHalfSpaces();
	  IntegerVectorList inequalities2;
	  for(IntegerVectorList::const_iterator k=inequalities1.begin();k!=inequalities1.end();k++)
	    {
	      if(dotLong(*j,*k))
		inequalities2.push_back(*k);
	      else
		equations.push_back(*k);
	    }
	  bool isFound=false;
	  for(IntegerVectorList::const_iterator j2=rays.begin();j2!=rays.end();j2++)
	    for(SymmetryGroup::ElementContainer::const_iterator k=sym->elements.begin();k!=sym->elements.end();k++)
	      {
		bool isInCone=true;
		IntegerVector v=SymmetryGroup::compose(*k,*j2);
		for(IntegerVectorList::const_iterator l=equations.begin();l!=equations.end();l++)
		  if(dotLong(*l,v)!=0)
		    {
		      isInCone=false;
		      goto leave;
		    }
		for(IntegerVectorList::const_iterator l=inequalities2.begin();l!=inequalities2.end();l++)
		  if(dotLong(*l,v)<0)
		    {
		      isInCone=false;
		      goto leave;
		    }
	      leave:
		if(isInCone)
		  {
		    isFound=true;
		    goto leave2;
		  }
	      }
	leave2:
	  if(!isFound)
	    {
	      IntegerVector ray=stableRay(*j,equations,inequalities2,sym);
	      rays.push_back(ray);
	    }
	}
    }
  rays.sort();
  if(upToSymmetry)return rays;
  IntegerVectorList ret;
  for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
    {
      set<IntegerVector> thisOrbitsRays;
      for(SymmetryGroup::ElementContainer::const_iterator k=sym->elements.begin();k!=sym->elements.end();k++)
	{
	  IntegerVector temp=SymmetryGroup::compose(*k,*i);
	  thisOrbitsRays.insert(temp);
	}
      for(set<IntegerVector>::const_iterator i=thisOrbitsRays.begin();i!=thisOrbitsRays.end();i++)ret.push_back(*i);
    }
  return ret;
}
#else
IntegerVectorList PolyhedralFan::getRaysInPrintingOrder(SymmetryGroup const *sym, bool upToSymmetry)const
{
  /*
   * The ordering changed in this version. Previously the orbit representatives stored in "rays" were
   * just the first extreme ray from the orbit that appeared. Now it is gotten using "orbitRepresentative"
   * which causes the ordering in which the different orbits appear to change.
   */

  if(cones.empty())return IntegerVectorList();
  IntegerVectorList generatorsOfLinealitySpace=cones.begin()->generatorsOfLinealitySpace();//all cones have the same lineality space

        SymmetryGroup localsym(n);
        if(!sym)sym=&localsym;
  set<IntegerVector> rays;
  log1 fprintf(Stderr,"Computing rays of %i cones\n",cones.size());
  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      {
        static int t;
        if(!((t++)%10))log1 fprintf(Stderr,"%i\n",t);
      }
      IntegerVectorList temp=i->extremeRays(&generatorsOfLinealitySpace);
      for(IntegerVectorList::const_iterator j=temp.begin();j!=temp.end();j++)
        rays.insert(sym->orbitRepresentative(*j));
    }
  IntegerVectorList ret;
  if(upToSymmetry)
    for(set<IntegerVector>::const_iterator i=rays.begin();i!=rays.end();i++)ret.push_back(*i);
  else
    for(set<IntegerVector>::const_iterator i=rays.begin();i!=rays.end();i++)
      {
        set<IntegerVector> thisOrbitsRays;
        for(SymmetryGroup::ElementContainer::const_iterator k=sym->elements.begin();k!=sym->elements.end();k++)
          {
            IntegerVector temp=SymmetryGroup::compose(*k,*i);
            thisOrbitsRays.insert(temp);
          }
        for(set<IntegerVector>::const_iterator i=thisOrbitsRays.begin();i!=thisOrbitsRays.end();i++)ret.push_back(*i);
      }
  return ret;
}


#endif


/*void PolyhedralFan::printWithIndices(class Printer *p, SymmetryGroup *sym)const //fan must be pure
{
  //  print(p);
  SymmetryGroup symm(n);
  if(!sym)sym=&symm;
  assert(!cones.empty());
  int h=cones.begin()->dimensionOfLargestContainedSubspace();
  fprintf(Stdout,"Rays:\n");
  //  IntegerVectorList rays;//=f.getRelativeInteriorPoints();
  fprintf(Stderr,"Computing rays...\n");
  IntegerVectorList rays=getRaysInPrintingOrder(sym);
  fprintf(Stderr,"Done computing rays.\n");
  p->printVectorList(rays,true);
  PolyhedralFan f=*this;

  //  while(f.getMaxDimension()>=h)
  IntegerVector fvector(f.getMaxDimension()-h);
  while(f.getMaxDimension()!=h)
    {
      int currentDimension=f.getMaxDimension()-h;
      fprintf(Stderr,"Processing dimension %i cones...\n",currentDimension);
      PolyhedralFan done(n);
      IntegerVector rayIncidenceCounter(rays.size());
      p->printString("Printing index list for dimension ");
      p->printInteger(currentDimension);
      p->printString(" cones:\n");
      p->printString("{\n");
      int faceIndex =0;
      bool split=false;
      bool addComma=false;
      for(PolyhedralConeList::const_iterator i=f.cones.begin();i!=f.cones.end();i++)
	{
	  if(!done.contains(*i))
	    {
	      for(SymmetryGroup::ElementContainer::const_iterator k=sym->elements.begin();k!=sym->elements.end();k++)
		{
		  PolyhedralCone cone=i->permuted(*k);
		  if(!done.contains(cone))
		    {
		      // p->printString("Face ");
		      // p->printInteger(faceIndex);
		      // p->printString(": ");
		      int rayIndex=0;
		      IntegerVector indices(0);
		      for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++)
			{
			  if(cone.contains(*j))
			    {
			      indices.grow(indices.size()+1);
			      indices[indices.size()-1]=rayIndex;
			      rayIncidenceCounter[rayIndex]++;
			    }
			  rayIndex++;
			}
		      if(addComma)
			{
			  p->printString(",");
			  p->printNewLine();
			  if(split)
			    p->printNewLine();
			  split=false;
			}
		      p->printVector(indices,true);
		      addComma=true;
		      faceIndex++;
		      done.insert(cone);
		    }
		}
	      split=true;//p->printNewLine();
	    }
	}
      p->printString("}\n");
      p->printString("Number of dimension ");
      p->printInteger(f.getMaxDimension()-h);
      p->printString(" cones incident to each ray:\n");
      p->printVector(rayIncidenceCounter);
      p->printNewLine();
      fvector[f.getMaxDimension()-h-1]=faceIndex;
      f=f.facetComplex();
      fprintf(Stderr,"Done processing dimension %i cones.\n",currentDimension);
      //      fvector.grow(fvector.size()+1);fvector[fvector.size()-1]=faceIndex;
    }
  p->printString("F-vector:\n");
  p->printVector(fvector);
  p->printNewLine();
}
*/

 /*
void PolyhedralFan::printWithIndices(class Printer *p, SymmetryGroup *sym, const char *polymakeFileName)const //fan must be pure
{
  IntegerVector multiplicities;
  SymmetryGroup symm(n);
  PolymakeFile polymakeFile;






  if(polymakeFileName)
    {
      polymakeFile.create(polymakeFileName,"PolyhedralFan","PolyhedralFan");
    }

  if(!sym)sym=&symm;
  assert(!cones.empty());
  int h=cones.begin()->dimensionOfLinealitySpace();
  fprintf(Stdout,"Rays:\n");
  //  IntegerVectorList rays;//=f.getRelativeInteriorPoints();
  fprintf(Stderr,"Computing rays...\n");
  IntegerVectorList rays=getRaysInPrintingOrder(sym);
  fprintf(Stderr,"Done computing rays.\n");


  SymmetricComplex symCom(n,rays,*sym);


  if(p)
    p->printVectorList(rays,true);
  if(polymakeFileName)
    {
      polymakeFile.writeCardinalProperty("AMBIENT_DIM",n);
      polymakeFile.writeCardinalProperty("DIM",getMaxDimension());
      polymakeFile.writeCardinalProperty("LINEALITY_DIM",h);
      polymakeFile.writeMatrixProperty("RAYS",rowsToIntegerMatrix(rays,n));
      polymakeFile.writeCardinalProperty("N_RAYS",rays.size());
      polymakeFile.writeMatrixProperty("LINEALITY_SPACE",rowsToIntegerMatrix(highestDimensionalCone().linealitySpace().dualCone().getEquations()));
      polymakeFile.writeMatrixProperty("ORTH_LINEALITY_SPACE",rowsToIntegerMatrix(highestDimensionalCone().linealitySpace().getEquations()));
      polymakeFile.writeCardinalProperty("PURE",1);
    }


  PolyhedralFan f=*this;
  stringstream s;

  //  while(f.getMaxDimension()>=h)
  IntegerVector fvector(f.getMaxDimension()-h);
  bool isHighestDimension=true;
  while(f.getMaxDimension()!=h)
    {
      int currentDimension=f.getMaxDimension()-h;
      fprintf(Stderr,"Processing dimension %i cones...\n",currentDimension);
      IntegerVector rayIncidenceCounter(rays.size());
      if(p)
	{
	  p->printString("Printing index list for dimension ");
	  p->printInteger(currentDimension);
	  p->printString(" cones:\n");
	  p->printString("{\n");
	}
      int faceIndex =0;
      bool split=false;
      bool addComma=false;
      for(PolyhedralConeList::const_iterator i=f.cones.begin();i!=f.cones.end();i++)
	{
	  {
	    SymmetricComplex::Cone c;
	    c.dimension=i->dimension();

	    int rayIndex=0;
	    for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++)
	      {
		if(i->contains(*j))c.insert(rayIndex);
		rayIndex++;
	      }
	    symCom.insert(c);
	  }


	  PolyhedralFan done(n);
	  for(SymmetryGroup::ElementContainer::const_iterator k=sym->elements.begin();k!=sym->elements.end();k++)
	    {
	      PolyhedralCone cone=i->permuted(*k);
	      if(!done.contains(cone))
		{
		  if(isHighestDimension)
		    {
		      multiplicities.grow(multiplicities.size()+1);
		      multiplicities[multiplicities.size()-1]=i->getMultiplicity();
		    }
		  // p->printString("Face ");
		  // p->printInteger(faceIndex);
		  // p->printString(": ");
		  int rayIndex=0;
		  IntegerVector indices(0);
		  for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++)
		    {
		      if(cone.contains(*j))
			{
			  indices.grow(indices.size()+1);
			  indices[indices.size()-1]=rayIndex;
			  rayIncidenceCounter[rayIndex]++;
			}
		      rayIndex++;
		    }
		  if(addComma)
		    {
		      if(p)
			{
			  p->printString(",");
			  p->printNewLine();
			}
		      if(isHighestDimension)s<<endl;
		      if(split)
			{
			  if(p)p->printNewLine();
			  //	  s<<endl;
			}
		      split=false;
		    }
		  if(p)
		    p->printVector(indices,true);
		  if(isHighestDimension)
		    {
		      s << '{';
		      for(int i=0;i<indices.size();i++)
			{
			  if(i)s<<" ";
			  s<<indices[i];
			}
		      s << '}';
		    }
		  addComma=true;
		  faceIndex++;
		  done.insert(cone);
		}
	    }
	  split=true;//p->printNewLine();
	}
      if(p)
	{
	  p->printString("}\n");
	  p->printString("Number of dimension ");
	  p->printInteger(f.getMaxDimension()-h);
	  p->printString(" cones incident to each ray:\n");
	  p->printVector(rayIncidenceCounter);
	  p->printNewLine();
	}
      fvector[f.getMaxDimension()-h-1]=faceIndex;
      fprintf(Stderr,"TESTESTSETST\n");
      f=f.facetComplexSymmetry(*sym);
      fprintf(Stderr,"TESTESTSETST\n");
      fprintf(Stderr,"Done processing dimension %i cones.\n",currentDimension);
      //      fvector.grow(fvector.size()+1);fvector[fvector.size()-1]=faceIndex;
      isHighestDimension=false;
    }
  if(p)
    {
      p->printString("Multiplicities:\n");
      p->printVector(multiplicities);
      p->printNewLine();
      p->printString("F-vector:\n");
      p->printVector(fvector);
      p->printNewLine();
    }

  if(polymakeFileName)
    {
      polymakeFile.writeCardinalVectorProperty("F_VECTOR",fvector);
      s<<endl;
      polymakeFile.writeStringProperty("MAXIMAL_CONES",s.str());

      IntegerVectorList m;
      m.push_back(multiplicities);
      polymakeFile.writeMatrixProperty("MULTIPLICITIES",rowsToIntegerMatrix(m).transposed());
    }
  if(polymakeFileName)
    polymakeFile.close();


  AsciiPrinter(Stdout).printString(symCom.toString(symCom.getMinDim(),symCom.getMaxDim(),true,true));
}
*/



  /*MARKS CONES AS NONMAXIMAL IN THE SYMMETRIC COMPLEX IN WHICH THEY WILL BE INSERTED -not to be confused with the facet testing in the code
   */
 static list<SymmetricComplex::Cone> computeFacets(SymmetricComplex::Cone const &theCone, IntegerMatrix const &rays, IntegerVectorList const &facetCandidates, SymmetricComplex const &theComplex/*, int linealityDim*/)
{
  set<SymmetricComplex::Cone> ret;

  for(IntegerVectorList::const_iterator i=facetCandidates.begin();i!=facetCandidates.end();i++)
    {
      set<int> indices;

      bool notAll=false;
      for(vector<int>::const_iterator j=theCone.indices.begin();j!=theCone.indices.end();j++)
	if(dotLong(rays[*j],*i)==0)
	  indices.insert(*j);
	else
	  notAll=true;

      SymmetricComplex::Cone temp(indices,theCone.dimension-1,0,false,theComplex);
      /*      temp.multiplicity=0;
      temp.dimension=theCone.dimension-1;
      temp.setIgnoreSymmetry(true);
      */
      if(notAll)ret.insert(temp);

    }
  //  fprintf(Stderr,"HEJ!!!!\n");

  list<SymmetricComplex::Cone> ret2;
  for(set<SymmetricComplex::Cone>::const_iterator i=ret.begin();i!=ret.end();i++)
    {
      bool isMaximal=true;

      /*      if(i->indices.size()+linealityDim<i->dimension)//#3
	isMaximal=false;
	else*/
	for(set<SymmetricComplex::Cone>::const_iterator j=ret.begin();j!=ret.end();j++)
	  if(i!=j && i->isSubsetOf(*j))
	    {
	      isMaximal=false;
	      break;
	    }
      if(isMaximal)
	{
	  SymmetricComplex::Cone temp(i->indexSet(),i->dimension,i->multiplicity,true,theComplex);
	  temp.setKnownToBeNonMaximal(); // THIS IS WHERE WE SET THE CONES NON-MAXIMAL FLAG
	  //	  temp.setIgnoreSymmetry(false);
	  ret2.push_back(temp);
	}
    }
  return ret2;
}


void addFacesToSymmetricComplex(SymmetricComplex &c, PolyhedralCone const &cone, IntegerVectorList const &facetCandidates, IntegerVectorList const &generatorsOfLinealitySpace)
{
  IntegerMatrix const &rays=c.getVertices();
  set<int> indices;

//  for(int j=0;j<rays.getHeight();j++)if(cone.contains(rays[j]))indices.insert(j);

  IntegerVectorList l=cone.extremeRays(&generatorsOfLinealitySpace);
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)indices.insert(c.indexOfVertex(*i));

  addFacesToSymmetricComplex(c,indices,facetCandidates,cone.dimension(),cone.getMultiplicity());
}

void addFacesToSymmetricComplex(SymmetricComplex &c, set<int> const &indices, IntegerVectorList const &facetCandidates, int dimension, int multiplicity)
{
  IntegerMatrix const &rays=c.getVertices();
  list<SymmetricComplex::Cone> clist;
  {

    SymmetricComplex::Cone temp(indices,dimension,multiplicity,true,c);
    //    temp.dimension=cone.dimension();
    //   temp.multiplicity=cone.getMultiplicity();
    clist.push_back(temp);
  }

  //  int linealityDim=cone.dimensionOfLinealitySpace();

  while(!clist.empty())
    {
      SymmetricComplex::Cone &theCone=clist.front();

      if(!c.contains(theCone))
	{
	  log2
	  {
	    static int t;
	    if((t&1023)==0)
	      {
		fprintf(Stderr,"clist size:%i\n",clist.size());
	      }
	    t++;
	  }

	  c.insert(theCone);
	  //	  log0 fprintf(Stderr,"ADDING\n");
	  list<SymmetricComplex::Cone> facets=computeFacets(theCone,rays,facetCandidates,c/*,linealityDim*/);
	  clist.splice(clist.end(),facets);
	}
      clist.pop_front();
    }

}


/**
   Produce strings that express the vectors in terms of rays of the fan modulo the lineality space. Symmetry is ignored??
 */
vector<string> PolyhedralFan::renamingStrings(IntegerVectorList const &theVectors, IntegerVectorList const &originalRays, IntegerVectorList const &linealitySpace, SymmetryGroup *sym)const
{
  vector<string> ret;
  for(IntegerVectorList::const_iterator i=theVectors.begin();i!=theVectors.end();i++)
    {
      for(PolyhedralConeList::const_iterator j=cones.begin();j!=cones.end();j++)
	{
	  if(j->contains(*i))
	    {
	      vector<int> relevantIndices;
	      IntegerVectorList relevantRays=linealitySpace;
	      int K=0;
	      for(IntegerVectorList::const_iterator k=originalRays.begin();k!=originalRays.end();k++,K++)
		if(j->contains(*k))
		  {
		    relevantIndices.push_back(K);
		    relevantRays.push_back(*k);
		  }

	      FieldMatrix LFA(Q,relevantRays.size(),n);
	      int J=0;
	      for(IntegerVectorList::const_iterator j=relevantRays.begin();j!=relevantRays.end();j++,J++)
		LFA[J]=integerVectorToFieldVector(*j,Q);
	      FieldVector LFB=concatenation(integerVectorToFieldVector(*i,Q),FieldVector(Q,relevantRays.size()));
	      LFA=LFA.transposed();
	      FieldVector LFX=LFA.solver().canonicalize(LFB);
	      stringstream s;
	      if(LFX.subvector(0,n).isZero())
	        {
		  s<<"Was:";
	          FieldVector S=LFX.subvector(n+linealitySpace.size(),LFX.size());
		  for(int k=0;k<S.size();k++)
		    if(!S[k].isZero())
		      s<<"+"<<S[k].toString()<<"*["<<relevantIndices[k]<<"] ";
	        }
	      ret.push_back(s.str());
	      break;
	    }
	}
    }
  return ret;
}

SymmetricComplex PolyhedralFan::toSymmetricComplex(SymmetryGroup *sym)
{
  SymmetryGroup symm(n);
	  if(!sym)sym=&symm;

	  IntegerVectorList rays=getRaysInPrintingOrder(sym);

	  SymmetricComplex symCom(n,rays,*sym);

	  if(cones.empty())return symCom;
	  IntegerVectorList generatorsOfLinealitySpace=cones.begin()->generatorsOfLinealitySpace();

	  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
	    {
	      {
		static int t;
		log1 fprintf(Stderr,"Adding faces of cone %i\n",t++);
	      }
	      log2 fprintf(Stderr,"Dim: %i\n",i->dimension());

	      addFacesToSymmetricComplex(symCom,*i,i->getHalfSpaces(),generatorsOfLinealitySpace);
	    }

	  log1 cerr<<"Remapping";
	  symCom.remap();
	  log1 cerr<<"Done remapping";
	  return symCom;
}

void PolyhedralFan::printWithIndices(class Printer *p, int flags, SymmetryGroup *sym, vector<string> const *comments)const
//void PolyhedralFan::printWithIndices(class Printer *p, bool printMultiplicities, SymmetryGroup *sym, bool group, bool ignoreCones, bool xml, bool tPlaneSort, vector<string> const *comments)const
{
  assert(p);
  //  IntegerVector multiplicities;
  SymmetryGroup symm(n);

  PolymakeFile polymakeFile;
//  polymakeFile.create("NONAME","fan","PolyhedralFan",flags&FPF_xml);
  polymakeFile.create("NONAME","fan","SymmetricFan",flags&FPF_xml);

  bool produceXml=polymakeFile.isXmlFormat();


  if(!sym)sym=&symm;

  if(cones.empty())
    {
      p->printString("Polyhedral fan is empty. Printing not supported.\n");
      return;
    }

  int h=cones.begin()->dimensionOfLinealitySpace();

  log1 fprintf(Stderr,"Computing rays.\n");
  IntegerVectorList rays=getRaysInPrintingOrder(sym);

  SymmetricComplex symCom(n,rays,*sym);

  polymakeFile.writeCardinalProperty("AMBIENT_DIM",n);
  polymakeFile.writeCardinalProperty("DIM",getMaxDimension());
  polymakeFile.writeCardinalProperty("LINEALITY_DIM",h);
  polymakeFile.writeMatrixProperty("RAYS",rowsToIntegerMatrix(rays,n),true,comments);
  polymakeFile.writeCardinalProperty("N_RAYS",rays.size());
  IntegerVectorList linealitySpaceGenerators=highestDimensionalCone().linealitySpace().dualCone().getEquations();
  polymakeFile.writeMatrixProperty("LINEALITY_SPACE",rowsToIntegerMatrix(linealitySpaceGenerators,n));
  polymakeFile.writeMatrixProperty("ORTH_LINEALITY_SPACE",rowsToIntegerMatrix(highestDimensionalCone().linealitySpace().getEquations(),n));

  if(flags & FPF_primitiveRays)
  {
	 IntegerVectorList primitiveRays;
	 for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
		 for(PolyhedralConeList::const_iterator j=cones.begin();j!=cones.end();j++)
			 if(j->contains(*i)&&(j->dimensionOfLinealitySpace()+1==j->dimension()))
					 primitiveRays.push_back(j->semiGroupGeneratorOfRay());

	  polymakeFile.writeMatrixProperty("PRIMITIVE_RAYS",rowsToIntegerMatrix(primitiveRays,n));
  }


  IntegerVectorList generatorsOfLinealitySpace=cones.begin()->generatorsOfLinealitySpace();

  log1 fprintf(Stderr,"Building symmetric complex.\n");
  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      {
	static int t;
	log1 fprintf(Stderr,"Adding faces of cone %i\n",t++);
      }
      log2 fprintf(Stderr,"Dim: %i\n",i->dimension());

      addFacesToSymmetricComplex(symCom,*i,i->getHalfSpaces(),generatorsOfLinealitySpace);
    }

  log1 cerr<<"Remapping";
  symCom.remap();
  log1 cerr<<"Done remapping";


  PolyhedralFan f=*this;

  //  IntegerVector fvector(f.getMaxDimension()-h);



  //fprintf(Stderr,"maxdim %i h %i\n",f.getMaxDimension(),h);
  /*  while(!f.cones.empty())
    {
      int currentDimension=f.getMaxDimension()-h;
      IntegerVector rayIncidenceCounter(rays.size());
      int faceIndex =0;
      bool split=false;
      bool addComma=false;
      for(PolyhedralConeList::const_iterator i=f.cones.begin();i!=f.cones.end();i++)
	{
	  {
	    SymmetricComplex::Cone c;
	    c.dimension=i->dimension();
	    c.multiplicity=i->getMultiplicity();

	    int rayIndex=0;
	    for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++)
	      {
		if(i->contains(*j))c.insert(rayIndex);
		rayIndex++;
	      }
	    symCom.insert(c);
	  }
	}
      //      fvector[f.getMaxDimension()-h-1]=faceIndex;
      f=f.facetComplexSymmetry(*sym);
    }
  */
  log1 fprintf(Stderr,"Computing f-vector.\n");
  IntegerVector fvector=symCom.fvector();
  polymakeFile.writeCardinalVectorProperty("F_VECTOR",fvector);
  log1 fprintf(Stderr,"Done computing f-vector.\n");

  if(flags&FPF_boundedInfo)
    {
      log1 fprintf(Stderr,"Computing bounded f-vector.\n");
      IntegerVector fvectorBounded=symCom.fvector(true);
      polymakeFile.writeCardinalVectorProperty("F_VECTOR_BOUNDED",fvectorBounded);
      log1 fprintf(Stderr,"Done computing bounded f-vector.\n");
    }
/*
 * Removed to make the Polymake people happy.
 *    {
    int euler=0;
    int mul=-1;
    for(int i=0;i<fvector.size();i++,mul*=-1)euler+=mul*fvector[i];
    polymakeFile.writeCardinalProperty("MY_EULER",euler);
  }
*/
  log1 fprintf(Stderr,"Checking if complex is simplicial and pure.\n");
  polymakeFile.writeBooleanProperty("SIMPLICIAL",symCom.isSimplicial());
  polymakeFile.writeBooleanProperty("PURE",symCom.isPure());
  log1 fprintf(Stderr,"Done checking.\n");


  if(flags&FPF_conesCompressed)
  {
    polymakeFile.writeArrayArrayIntProperty("SYMMETRY_GENERATORS",rowsToIntegerMatrix(sym->getUniqueGenerators(),n));
    log1 fprintf(Stderr,"Producing list of cones up to symmetry.\n");
    polymakeFile.writeStringProperty("CONES_ORBITS",symCom.toString(symCom.getMinDim(),symCom.getMaxDim(),false,flags&FPF_group,0,true,flags&FPF_tPlaneSort,produceXml));
    log1 fprintf(Stderr,"Done producing list of cones up to symmetry.\n");
    log1 fprintf(Stderr,"Producing list of maximal cones up to symmetry.\n");
    stringstream multiplicities;
    polymakeFile.writeStringProperty("MAXIMAL_CONES_ORBITS",symCom.toString(symCom.getMinDim(),symCom.getMaxDim(),true,flags&FPF_group, &multiplicities,true,flags&FPF_tPlaneSort,produceXml));
    if(flags&FPF_multiplicities)polymakeFile.writeStringProperty("MULTIPLICITIES_ORBITS",multiplicities.str());
    log1 fprintf(Stderr,"Done producing list of maximal cones up to symmetry.\n");
  }

  if(flags&FPF_conesExpanded)
    {
      if(flags&FPF_cones)
	{
	  log1 fprintf(Stderr,"Producing list of cones.\n");
	  polymakeFile.writeStringProperty("CONES",symCom.toString(symCom.getMinDim(),symCom.getMaxDim(),false,flags&FPF_group,0,false,flags&FPF_tPlaneSort,produceXml));
	  log1 fprintf(Stderr,"Done producing list of cones.\n");
	}
      if(flags&FPF_maximalCones)
	{
	  log1 fprintf(Stderr,"Producing list of maximal cones.\n");
	  stringstream multiplicities;
	  polymakeFile.writeStringProperty("MAXIMAL_CONES",symCom.toString(symCom.getMinDim(),symCom.getMaxDim(),true,flags&FPF_group, &multiplicities,false,flags&FPF_tPlaneSort,produceXml));
	  if(flags&FPF_multiplicities)polymakeFile.writeStringProperty("MULTIPLICITIES",multiplicities.str());
	  log1 fprintf(Stderr,"Done producing list of maximal cones.\n");
	}
    }

  if(flags&FPF_values)
    {
      {
	IntegerVectorList values;
	for(IntegerVectorList::const_iterator i=linealitySpaceGenerators.begin();i!=linealitySpaceGenerators.end();i++)
	  {
	    IntegerVector v(1);
	    v[0]=evaluatePiecewiseLinearFunction(*i);
	    values.push_back(v);
	  }
	polymakeFile.writeMatrixProperty("LINEALITY_VALUES",rowsToIntegerMatrix(values,1));
      }
      {
	IntegerVectorList values;
	for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
	  {
	    IntegerVector v(1);
	    v[0]=evaluatePiecewiseLinearFunction(*i);
	    values.push_back(v);
	  }
	polymakeFile.writeMatrixProperty("RAY_VALUES",rowsToIntegerMatrix(values,1));
      }
    }


  log1 fprintf(Stderr,"Producing final string for output.\n");
  stringstream s;
  polymakeFile.writeStream(s);
  string S=s.str();
  log1 fprintf(Stderr,"Printing string.\n");
  p->printString(S.c_str());
  log1 fprintf(Stderr,"Done printing string.\n");
}


PolyhedralFan PolyhedralFan::readFan(string const &filename, bool onlyMaximal, IntegerVector *w, set<int> const *coneIndices, SymmetryGroup const *sym, bool readCompressedIfNotSym)
{
    PolymakeFile inFile;
    inFile.open(filename.c_str());

    int n=inFile.readCardinalProperty("AMBIENT_DIM");
    int nRays=inFile.readCardinalProperty("N_RAYS");
    IntegerMatrix rays=inFile.readMatrixProperty("RAYS",nRays,n);
    int linealityDim=inFile.readCardinalProperty("LINEALITY_DIM");
    IntegerMatrix linealitySpace=inFile.readMatrixProperty("LINEALITY_SPACE",linealityDim,n);


    const char *sectionName=0;
    const char *sectionNameMultiplicities=0;
    if(sym || readCompressedIfNotSym)
    {
      sectionName=(onlyMaximal)?"MAXIMAL_CONES_ORBITS":"CONES_ORBITS";
      sectionNameMultiplicities=(onlyMaximal)?"MULTIPLICITIES_ORBITS":"DUMMY123";
    }
      else
      {  sectionName=(onlyMaximal)?"MAXIMAL_CONES":"CONES";
      sectionNameMultiplicities=(onlyMaximal)?"MULTIPLICITIES":"DUMMY123";
      }


    IntegerVector w2(n);
    if(w==0)w=&w2;

    SymmetryGroup sym2(n);
    if(sym==0)sym=&sym2;

    vector<list<int> > cones=inFile.readMatrixIncidenceProperty(sectionName);
    IntegerVectorList r;

    bool hasMultiplicities=inFile.hasProperty(sectionNameMultiplicities);
    IntegerMatrix multiplicities(0,0);
    if(hasMultiplicities)multiplicities=inFile.readMatrixProperty(sectionNameMultiplicities,cones.size(),1);


    PolyhedralFan ret(n);

    log2 cerr<< "Number of orbits to expand "<<cones.size()<<endl;
    for(int i=0;i<cones.size();i++)
      if(coneIndices==0 || coneIndices->count(i))
	{
	  log2 cerr<<"Expanding symmetries of cone"<<endl;
	  /*	  for(SymmetryGroup::ElementContainer::const_iterator perm=sym->elements.begin();perm!=sym->elements.end();perm++)
	    {
	      IntegerVectorList coneRays;
	      for(list<int>::const_iterator j=cones[i].begin();j!=cones[i].end();j++)
		coneRays.push_back(SymmetryGroup::compose(*perm,rays[*j]));
	      if(isInNonNegativeSpan(*w,coneRays,linealitySpace.getRows()))
		{
		  PolyhedralCone C=PolyhedralCone::givenByRays(coneRays,linealitySpace.getRows(),n);
		  C.canonicalize();
		  ret.insert(C);
		}
	    }
	  */
	  {
	    IntegerVectorList coneRays;
	    for(list<int>::const_iterator j=cones[i].begin();j!=cones[i].end();j++)
	      coneRays.push_back((rays[*j]));
	    PolyhedralCone C=PolyhedralCone::givenByRays(coneRays,linealitySpace.getRows(),n);
	    if(hasMultiplicities)C.setMultiplicity(multiplicities[i][0]);
	    for(SymmetryGroup::ElementContainer::const_iterator perm=sym->elements.begin();perm!=sym->elements.end();perm++)
	      {
		if(C.contains(SymmetryGroup::composeInverse(*perm,*w)))
		  {
		    PolyhedralCone C2=C.permuted(*perm);
		    C2.canonicalize();
		    ret.insert(C2);
		  }
	      }
	  }
	}
    return ret;
}


IncidenceList PolyhedralFan::getIncidenceList(SymmetryGroup *sym)const //fan must be pure
{
  IncidenceList ret;
  SymmetryGroup symm(n);
  if(!sym)sym=&symm;
  assert(!cones.empty());
  int h=cones.begin()->dimensionOfLinealitySpace();
  IntegerVectorList rays=getRaysInPrintingOrder(sym);
  PolyhedralFan f=*this;

  while(f.getMaxDimension()!=h)
    {
      IntegerVectorList l;
      PolyhedralFan done(n);
      IntegerVector rayIncidenceCounter(rays.size());
      int faceIndex =0;
      for(PolyhedralConeList::const_iterator i=f.cones.begin();i!=f.cones.end();i++)
	{
	  if(!done.contains(*i))
	    {
	      for(SymmetryGroup::ElementContainer::const_iterator k=sym->elements.begin();k!=sym->elements.end();k++)
		{
		  PolyhedralCone cone=i->permuted(*k);
		  if(!done.contains(cone))
		    {
		      int rayIndex=0;
		      IntegerVector indices(0);
		      for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++)
			{
			  if(cone.contains(*j))
			    {
			      indices.grow(indices.size()+1);
			      indices[indices.size()-1]=rayIndex;
			      rayIncidenceCounter[rayIndex]++;
			    }
			  rayIndex++;
			}
		      l.push_back(indices);
		      faceIndex++;
		      done.insert(cone);
		    }
		}
	    }
	}
      ret[f.getMaxDimension()]=l;
      f=f.facetComplex();
    }
  return ret;
}


void PolyhedralFan::makePure()
{
  if(getMaxDimension()!=getMinDimension())removeAllLowerDimensional();
}

bool PolyhedralFan::contains(PolyhedralCone const &c)const
{
  return cones.count(c);
}


PolyhedralCone PolyhedralFan::coneContaining(IntegerVector const &v)const
{
  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    if(i->contains(v))return i->faceContaining(v);
  debug<<"Vector "<<v<<" not contained in support of fan\n";
  assert(0);
}


PolyhedralFan::coneIterator PolyhedralFan::conesBegin()const
{
  return cones.begin();
}


PolyhedralFan::coneIterator PolyhedralFan::conesEnd()const
{
  return cones.end();
}


bool PolyhedralFan::isRefinementOf(PolyhedralFan const &f)const
{
  /*  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      static int t;
      fprintf(Stderr,"%i\n",t++);
      for(PolyhedralConeList::const_iterator j=f.cones.begin();j!=f.cones.end();j++)
	{
	  PolyhedralCone c=intersection(*i,*j);
	  if(c.dimension()==n)
	    {

	      if(!j->contains(*i))
		return false;
	    }
	}
	}*/

  int t=0;
  for(PolyhedralConeList::const_iterator j=f.cones.begin();j!=f.cones.end();j++)
    {
      int t1=0;
      for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
	{
	  PolyhedralCone c=intersection(*i,*j);
	  if(c.dimension()==n)
	    {
	      t1++;
	      if(!j->contains(*i))
		return false;
	    }
	}
      fprintf(Stderr,"%i\n",t1);
      t+=t1;
    }
  fprintf(Stderr,"%i\n",t);

      /*      bool found=false;
      for(PolyhedralConeList::const_iterator j=f.cones.begin();j!=f.cones.end();j++)
	  if(j->contains(*i)){
	    found=true;
	    break;
	  }
      if(!found)
	{
	  for(PolyhedralConeList::const_iterator j=f.cones.begin();j!=f.cones.end();j++)
	    {
	      PolyhedralCone c=intersection(*i,*j);
	      if(c.dimension()==n){
		AsciiPrinter(Stderr).printVector(c.getRelativeInteriorPoint());
	      }
	    }

	  return false;
	  }*/

  return true;
}

/*	    for(SymmetryGroup::ElementContainer::const_iterator perm=sym->elements.begin();perm!=sym->elements.end();perm++)
	      {
		if(C.contains(SymmetryGroup::composeInverse(*perm,*w)))
		  {
		    PolyhedralCone C2=C.permuted(*perm);
		    C2.canonicalize();
		    ret.insert(C2);
		  }
	      }
*/
PolyhedralFan PolyhedralFan::link(IntegerVector const &w, SymmetryGroup *sym)const
{
  SymmetryGroup symL(n);
  if(!sym)sym=&symL;

  PolyhedralFan ret(n);

  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      for(SymmetryGroup::ElementContainer::const_iterator perm=sym->elements.begin();perm!=sym->elements.end();perm++)
	{
	  IntegerVector w2=SymmetryGroup::composeInverse(*perm,w);
	  if(i->contains(w2))
	    {
	      IntegerVectorList equations=i->getEquations();
	      IntegerVectorList inequalities1=i->getHalfSpaces();
	      IntegerVectorList inequalities2;
	      for(IntegerVectorList::const_iterator j=inequalities1.begin();j!=inequalities1.end();j++)
		if(dotLong(w2,*j)==0)inequalities2.push_back(*j);
	      PolyhedralCone C(inequalities2,equations,n);
	      C.canonicalize();
	      C.setLinearForm(i->getLinearForm());
	      PolyhedralCone C2=C.permuted(*perm);
	      C2.canonicalize();
	      C2.setMultiplicity(i->getMultiplicity());
	      ret.insert(C2);
	    }
	}
    }
  return ret;
}

PolyhedralFan PolyhedralFan::link(IntegerVector const &w)const
{
  PolyhedralFan ret(n);

  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      if(i->contains(w))
	{
	  IntegerVectorList equations=i->getEquations();
	  IntegerVectorList inequalities1=i->getHalfSpaces();
	  IntegerVectorList inequalities2;
	  for(IntegerVectorList::const_iterator j=inequalities1.begin();j!=inequalities1.end();j++)
	    if(dotLong(w,*j)==0)inequalities2.push_back(*j);
	  PolyhedralCone C(inequalities2,equations,n);
	  C.canonicalize();
	  C.setLinearForm(i->getLinearForm());
      C.setMultiplicity(i->getMultiplicity());
	  ret.insert(C);
	}
    }
  return ret;
}


int64 PolyhedralFan::evaluatePiecewiseLinearFunction(IntegerVector const &x)const
{
  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      if(i->contains(x))return dotLong(i->getLinearForm(),x);
    }
  assert(0);
  return 0;
}


FieldElement PolyhedralFan::volume(int d, SymmetryGroup *sym)const
{
  SymmetryGroup symL(n);
  if(!sym)sym=&symL;

  FieldElement ret(Q);

  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      if(i->dimension()==d)
	{
	  IntegerVector w=stableRay(*i,sym);
	  ret=ret+Q.zHomomorphism(sym->orbitSize(w))*i->volume();
	}
    }
  return ret;
}


bool PolyhedralFan::isConnected(SymmetryGroup *sym)const
{
  SymmetryGroup symL(n);
  if(!sym)sym=&symL;

  CodimOneConnectednessTester ct;

  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      log2 cerr<<"Computing ridges of facet." << endl;
      PolyhedralFan ridges=facetsOfCone(*i);
      IntegerVectorList interiorPoints;
      for(PolyhedralConeList::const_iterator j=ridges.cones.begin();j!=ridges.cones.end();j++)
	interiorPoints.push_back(sym->orbitRepresentative(j->getUniquePoint()));
      ct.insertFacetOrbit(interiorPoints);
    }
  return ct.isConnected();
}


int PolyhedralFan::size()const
{
  return cones.size();
}

int PolyhedralFan::dimensionOfLinealitySpace()const
{
  assert(cones.size());//slow!
  return cones.begin()->dimensionOfLinealitySpace();
}

PolyhedralFan PolyhedralFan::negated()const
{
  PolyhedralFan ret(n);

  for(coneIterator i=conesBegin();i!=conesEnd();i++)
    ret.insert(i->negated());
  return ret;
}


bool PolyhedralFan::isCompatible(PolyhedralCone const &c)const
{
  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      PolyhedralCone C=intersection(c,*i);
      C.canonicalize();
      if(!c.hasFace(C))return false;
      if(!i->hasFace(C))return false;
    }
  return true;
}

void PolyhedralFan::merge(PolyhedralCone const &c)
{
  AsciiPrinter P(Stderr);

  if(isCompatible(c))insert(c);
  else
    {
      assert(0);//Does not work in general.
    }
  /*
  //  P<<"BEFORE MERGE-------------------------" <<*this;

  PolyhedralFan ret=complementOfCone(c,false);
  //  P<<"COMPLEMENT OF CONE-------------------------" <<ret;
  ret=refinement(*this,ret);

  PolyhedralFan C(c.ambientDimension());
  C.insert(c);

  for(PolyhedralConeList::const_iterator i=cones.begin();i!=cones.end();i++)
    C=refinement(complementOfCone(*i,true),C);

  for(PolyhedralConeList::const_iterator i=C.cones.begin();i!=C.cones.end();i++)
    ret.insert(*i);
  *this=ret;

  // P<<"AFTER MERGE" <<*this;
  */
}



void PolyhedralFan::removeNonMaximal()
{
  for(PolyhedralConeList::iterator i=cones.begin();i!=cones.end();)
    {
      IntegerVector w=i->getRelativeInteriorPoint();
      bool containedInOther=false;
      for(PolyhedralConeList::iterator j=cones.begin();j!=cones.end();j++)
	if(j!=i)
	  {
	    if(j->contains(w)){containedInOther=true;break;}
	  }
      if(containedInOther)
	{
	  PolyhedralConeList::iterator k=i;
	  i++;
	  cones.erase(k);
	}
      else i++;
    }
}


