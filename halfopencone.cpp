#include "halfopencone.h"

#include <iostream>
#include "buchberger.h"
#include "enumeration.h"
#include "reversesearch.h"
#include "wallideal.h"

#include "printer.h"
#include "parser.h"

#include "newtonpolytope.h"

#include "lp.h"
#include "log.h"

static void saveList(IntegerVectorList const &l, char const *name)
{
  FILE *f=fopen(name,"w");
  if(f)
    {
      AsciiPrinter P(f);
      P<<l;
      fclose(f);
    }
}

class HalfOpenConeProcessorVoid : public HalfOpenConeProcessor
{
public:
  void process(HalfOpenCone const &c)
  {
  }
};

/* This Cone processor saves memory when doing up to symmetry computations where orbits can be repeated */
class HalfOpenConeProcessorPointCollector : public HalfOpenConeProcessor
{
  int counter;
public:
  IntegerVectorList interiorPoints;
  HalfOpenConeProcessorPointCollector():
    counter(0)
  {
  }
  void process(HalfOpenCone const &c)
  {
    HalfOpenCone c1=c;
    PolyhedralCone c2=c1.closure();
    c2.canonicalize();
    interiorPoints.push_back(c2.getRelativeInteriorPoint());

    if(savePartialResult)
      {
        assert(0);
	counter++;
	if(counter>=1000)
	  saveList(interiorPoints,"partialresult");
	counter=0;
      }
  }
};

/* This Cone processor saves memory when doing up to symmetry computations where orbits can be repeated */
class HalfOpenConeProcessorConeCollector : public HalfOpenConeProcessor
{
public:
  HalfOpenConeList theList;
  void process(HalfOpenCone const &c)
  {
    theList.push_back(c);
  }
};


static void printHalfOpenCone(Printer &P, HalfOpenCone c)
{
  P.printPolyhedralCone(c.closure());
}


static void printHalfOpenConeList(Printer &P, HalfOpenConeList const &l)
{
  P.printString("Begin HalfOpenConeList\n");
  for(HalfOpenConeList::const_iterator i=l.begin();i!=l.end();i++)
    printHalfOpenCone(P,*i);
  P.printString("End HalfOpenConeList\n");
}


bool HalfOpenCone::contains(IntegerVector const &v)const
{
  IntegerVectorList inequalityList=lifted.getHalfSpaces();
  IntegerVectorList equationList=lifted.getEquations();

  IntegerVectorList strict,nonstrict;


  for(IntegerVectorList::const_iterator i=inequalityList.begin();i!=inequalityList.end();i++)
    if((*i)[i->size()-1]<0)
      strict.push_back(*i);
    else if((*i)[i->size()-1]==0)
      nonstrict.push_back(*i);
    else
      {//CHANGED
	assert(i->subvector(0,i->size()-1).isZero());
	strict.push_back(*i);
      }

  for(IntegerVectorList::const_iterator i=equationList.begin();i!=equationList.end();i++)
    if(dotLong(i->subvector(0,i->size()-1),v)!=0)return false;

  for(IntegerVectorList::const_iterator i=nonstrict.begin();i!=nonstrict.end();i++)
    if(dotLong(i->subvector(0,i->size()-1),v)<0)return false;

  for(IntegerVectorList::const_iterator i=strict.begin();i!=strict.end();i++)
    if(dotLong(i->subvector(0,i->size()-1),v)<=0)return false;

  return true;
}

void HalfOpenCone::appendList(IntegerVectorList &to, IntegerVectorList const &from, int appendValue)
{
  for(IntegerVectorList::const_iterator i=from.begin();i!=from.end();i++)
    {
      IntegerVector v=*i;
      v.resize(v.size()+1);
      v[v.size()-1]=appendValue;
      to.push_back(v);
    }
}


HalfOpenCone::HalfOpenCone(int dimension_, PolyhedralCone const &lifted_):
  dimension(dimension_),
  liftedDimension(dimension_+1),
  lifted(lifted_)
{
  //  lifted.findFacets();
}


HalfOpenCone::HalfOpenCone(int dimension_, IntegerVectorList const &equations, IntegerVectorList const &nonstrict, IntegerVectorList const &strict, bool findFacets, bool canonicalize):
  dimension(dimension_),
  liftedDimension(dimension_+1),
  lifted(dimension_+1)
{
  IntegerVectorList equationList,inequalityList;

  appendList(equationList,equations,0);
  appendList(inequalityList,nonstrict,0);
  appendList(inequalityList,strict,-1);
  inequalityList.push_back(IntegerVector::standardVector(liftedDimension,dimension)); //CHANGED
  //      AsciiPrinter(Stderr).printVectorList(inequalityList);
  //      AsciiPrinter(Stderr).printVectorList(equationList);
  //      AsciiPrinter(Stderr).printInteger(liftedDimension);
  lifted=PolyhedralCone(inequalityList,equationList,liftedDimension);
  if(findFacets)lifted.findFacets();
  if(canonicalize)lifted.canonicalize();
}


HalfOpenCone::HalfOpenCone(PolyhedralCone C, TermOrder const &t):
  dimension(C.ambientDimension()),
  liftedDimension(C.ambientDimension()+1),
  lifted(C.ambientDimension()+1)
{
  HalfOpenConeList ret;
  C.findFacets();
  assert(C.dimension()==C.ambientDimension());

  IntegerVectorList facets=C.getHalfSpaces();

  IntegerVectorList strictList,nonStrictList;
  for(IntegerVectorList::const_iterator i=facets.begin();i!=facets.end();i++)
    {
      if(t(*i,*i-*i))
	strictList.push_back(*i);
      else
	nonStrictList.push_back(*i);
    }

  IntegerVectorList inequalityList;
  appendList(inequalityList,nonStrictList,0);
  appendList(inequalityList,strictList,-1);
  inequalityList.push_back(IntegerVector::standardVector(liftedDimension,dimension)); //CHANGED
  IntegerVectorList empty;
  lifted=PolyhedralCone(inequalityList,empty,liftedDimension);
  //  if(findFacets)lifted.findFacets();
}


HalfOpenCone::HalfOpenCone(int ambientDimension):
  dimension(ambientDimension),
  liftedDimension(ambientDimension+1),
  lifted(ambientDimension+1)
{
  IntegerVectorList inequalityList;
  inequalityList.push_back(IntegerVector::standardVector(liftedDimension,dimension));
  IntegerVectorList empty;
  lifted=PolyhedralCone(inequalityList,empty,liftedDimension);
}


static IntegerVectorList swapFirstLast(const IntegerVectorList &l)
{
  IntegerVectorList ret;

  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      IntegerVector v=*i;
      int t=v[0];
      v[0]=v[v.size()-1];
      v[v.size()-1]=t;
      ret.push_back(v);
    }

  return ret;
}


bool HalfOpenCone::isEmpty()
{
    bool ret1=!hasHomogeneousSolution(liftedDimension,
				   swapFirstLast(lifted.getHalfSpaces()),
				   swapFirstLast(lifted.getEquations())
				   );
    /*
  IntegerVectorList inequalityList;
  inequalityList.push_back(IntegerVector::standardVector(liftedDimension,dimension));
  PolyhedralCone temp=intersection(lifted,PolyhedralCone(inequalityList,IntegerVectorList(),liftedDimension));
  IntegerVector v=temp.getRelativeInteriorPoint();
  //  AsciiPrinter(Stderr).printVector(v);
  bool ret2=(v[dimension]==0);
    */

  /*  fprintf(Stderr,"Inequalities:\n");
  AsciiPrinter(Stderr).printVectorList(lifted.getHalfSpaces());
  fprintf(Stderr,"Equations:\n");
  AsciiPrinter(Stderr).printVectorList(lifted.getEquations());
  fprintf(Stderr,"hasSolution=%i\n",ret1);
  */
  //  assert(ret1==ret2);

  return ret1;
}


bool haveEmptyIntersection(const HalfOpenCone &a, const HalfOpenCone &b)
{
  assert(a.dimension==b.dimension);
  IntegerVectorList inequalityList=a.lifted.getHalfSpaces();
  IntegerVectorList equationList=a.lifted.getEquations();

  IntegerVectorList inequalityList2=b.lifted.getHalfSpaces();
  IntegerVectorList equationList2=b.lifted.getEquations();

  inequalityList.splice(inequalityList.begin(),inequalityList2);
  equationList.splice(equationList.begin(),equationList2);

  bool ret1=!hasHomogeneousSolution(a.liftedDimension,swapFirstLast(inequalityList),swapFirstLast(equationList));

  /*  HalfOpenCone c=intersection(a,b);
  if(c.isEmpty()!=ret1)
    {
      AsciiPrinter(Stderr).printVectorList(inequalityList);
      AsciiPrinter(Stderr).printVectorList(equationList);

      AsciiPrinter(Stderr).printVectorList(c.lifted.getHalfSpaces());
      AsciiPrinter(Stderr).printVectorList(c.lifted.getEquations());

      fprintf(Stderr,"hasHomogeneousSolution siger %i\n",!ret1);

      assert(0);
    }
  */

  return ret1;
}


/*bool HalfOpenCone::isEmpty()
{
  IntegerVector v(liftedDimension);
  v[dimension]=-1;

  IntegerVectorList equationList,inequalityList;
  inequalityList.push_back(v);
  PolyhedralCone c(inequalityList,equationList,liftedDimension);
  PolyhedralCone c2=intersection(c,lifted);
  lifted.canonicalize();
  c2.canonicalize();

  return !(c2!=lifted);
  }*/

HalfOpenCone intersection(const HalfOpenCone &a, const HalfOpenCone &b, bool findFacets)
{
  assert(a.dimension==b.dimension);

  /*  fprintf(Stderr,"-----------------------------------------------------------\n");
  fprintf(Stderr,"Intersecting:\n");
  AsciiPrinter P(Stderr);
  printHalfOpenCone(P,a);
  printHalfOpenCone(P,b);
  */



  HalfOpenCone ret=HalfOpenCone(a.dimension,intersection(a.lifted,b.lifted));

  {
    static int t;
    t++;
    if((!(t&7))||findFacets)ret.lifted.findFacets();

    //1 4:53
    //3 3:38
    //7
  }


  /*  fprintf(Stderr,"Result:\n");
  printHalfOpenCone(P,ret);
  fprintf(Stderr,"Is empty:%i\n",ret.isEmpty());
  fprintf(Stderr,"-----------------------------------------------------------\n");
  fprintf(Stderr,"States: %i,%i,%i\n",a.lifted.getState(),b.lifted.getState(),ret.lifted.getState());
  fprintf(Stderr,"-----------------------------------------------------------\n");
  */

  return ret;
}


IntegerVectorList HalfOpenCone::shrink(const IntegerVectorList &l)
{
  IntegerVectorList ret;
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    ret.push_back(i->subvector(0,i->size()-1));
  return ret;
}


PolyhedralCone HalfOpenCone::closure()
{
  lifted.findFacets();
  return PolyhedralCone(shrink(lifted.getHalfSpaces()),shrink(lifted.getEquations()),dimension);
}


/*
 */
HalfOpenConeList orientedBoundary(PolyhedralCone C, TermOrder const &t, HalfOpenCone *restrictingCone)
{
  int dimension=C.ambientDimension();
  HalfOpenConeList ret;
  C.findFacets();
  assert(C.dimension()==C.ambientDimension());

  IntegerVectorList facets=C.getHalfSpaces();

  IntegerVectorList strictList,nonStrictList;
  for(IntegerVectorList::const_iterator i=facets.begin();i!=facets.end();i++)
    {
      if(t(*i,*i-*i))
	strictList.push_back(*i);
      else
	nonStrictList.push_back(*i);
    }
  //  log0 AsciiPrinter(Stderr).printVectorList(strictList);
  //  log0 AsciiPrinter(Stderr).printVectorList(nonStrictList);
  // Let's make the non-strict inequalities strict one at a time and add a cone for each iteration
  while(!nonStrictList.empty())
    {
      //      fprintf(Stderr,"NEWWALL\n");
      IntegerVector v=nonStrictList.front();
      nonStrictList.pop_front();


      IntegerVectorList equationList;
      equationList.push_back(v);
      {
	HalfOpenCone c(dimension,equationList,nonStrictList,strictList,true);
	if(restrictingCone)
	  c=intersection(*restrictingCone,c,true);
	if(!c.isEmpty())ret.push_back(c);
      }
      strictList.push_back(v);
    }
  return ret;
}


HalfOpenConeList tropicalHyperSurface(Polynomial const &p1)
{
  PolynomialRing theRing=p1.getRing();
  PolynomialRing theSecondRing=theRing.withVariablesAppended("H");
  Polynomial p=p1.homogenization(theSecondRing);
  HalfOpenConeList ret;
  PolynomialSet g(theRing);
  g.push_back(p);
  buchberger(&g,LexicographicTermOrder());

  EnumerationTargetCollector gfan;
  LexicographicTermOrder myTermOrder;
  ReverseSearch rs(myTermOrder);

  rs.setEnumerationTarget(&gfan);

  fprintf(Stderr,"Starting enumeratioin\n");
  rs.enumerate(g);
  fprintf(Stderr,"Done\n");

  PolynomialSetList theList=gfan.getList();
  for(PolynomialSetList::const_iterator i=theList.begin();i!=theList.end();i++)
    {
      HalfOpenConeList temp=orientedBoundary(groebnerCone(i->deHomogenization(),false),myTermOrder);
      ret.splice(ret.begin(),temp);
    }

  //  AsciiPrinter P(Stderr);
  //  printHalfOpenConeList(P,ret);

  return ret;
}


static void buildFanFromPolynomial(Polynomial const &p1, HalfOpenConeList *fullDimCones, HalfOpenConeList *coDimOneCones, IntegerVector *parentList, HalfOpenCone *restrictingCone=0)
{

  int n(p1.getRing().getNumberOfVariables());

  IntegerVectorList empty;
  HalfOpenCone dummy(n,empty,empty,empty,true);
  if(restrictingCone==0)restrictingCone=&dummy;

  if(p1.isZero()) //If the polynomial is zero it makes most sense to return a complete fan consisting of one cone in both cases
    {
      IntegerVectorList empty;
      if(fullDimCones)fullDimCones->push_back(intersection(HalfOpenCone(n,empty,empty,empty),*restrictingCone));
      if(coDimOneCones)coDimOneCones->push_back(intersection(HalfOpenCone(n,empty,empty,empty),*restrictingCone));
      return;
    }

  LexicographicTermOrder myTermOrder;
  IntegerVectorList l=newtonPolytope(p1);

  //  log0 AsciiPrinter(Stderr).printVectorList(l);
  assert(!l.empty());
  removeNonExtremeVerticesOfPolytope(l);
  //  log0 AsciiPrinter(Stderr).printVectorList(l);
  IntegerVector parents;
  int numberOfCoDimOneCones=0;
  int numberOfFullDimCones=0;
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      numberOfFullDimCones++;
      IntegerVectorList inequalities;
      for(IntegerVectorList::const_iterator j=l.begin();j!=l.end();j++)
	if(j!=i)
	  inequalities.push_back(*i-*j);

      //      log0 fprintf(Stderr,"ineq\n");
      //  log0 AsciiPrinter(Stderr).printVectorList(inequalities);
      inequalities=normalizedWithSumsAndDuplicatesRemoved(inequalities);
      //log0 fprintf(Stderr,"ineq\n");
      //log0 AsciiPrinter(Stderr).printVectorList(inequalities);

      IntegerVectorList empty;
      PolyhedralCone C(inequalities,empty,n,PCP_impliedEquationsKnown);


      C.findFacets();

      if(coDimOneCones || parentList)
	{
	  HalfOpenConeList temp=orientedBoundary(C,myTermOrder,restrictingCone);
	  //	  fprintf(Stderr,"TESTSET:%i\n",temp.size());
	  for(int i=0;i<temp.size();i++)
	    {
	      numberOfCoDimOneCones++;
	      parents.resize(numberOfCoDimOneCones);
	      parents[numberOfCoDimOneCones-1]=numberOfFullDimCones-1;
	    }
	  if(coDimOneCones)
	    {
	      /*	  fprintf(Stderr,"-------------------------------------------------\n");
	      for(HalfOpenConeList::const_iterator k=temp.begin();k!=temp.end();k++)
		{
		  AsciiPrinter P(Stderr);
		  k->print(P);
		}
		fprintf(Stderr,"-------------------------------------------------\n");*/
	      HalfOpenConeList::iterator k=coDimOneCones->end();
	      coDimOneCones->splice(k,temp);
	    }
	}
      if(fullDimCones)
	{
	  fullDimCones->push_back(intersection(HalfOpenCone(C,myTermOrder),*restrictingCone));
	  /*	  fprintf(Stderr,"Came from-------------------------------------------------\n");
	  AsciiPrinter P(Stderr);
	  fullDimCones->back().print(P);
	  fprintf(Stderr,"-------------------------------------------------\n");
	  */
	}
    }

  if(parentList)
    *parentList=parents;
}


HalfOpenConeList tropicalHyperSurfaceFast(Polynomial const &p1)
{
  HalfOpenConeList ret;
  buildFanFromPolynomial(p1,0,&ret,0);
  return ret;
}


HalfOpenConeList normalFanOfNewtonPolytope(Polynomial const &p1)
{
  HalfOpenConeList ret;
  buildFanFromPolynomial(p1,&ret,0,0);
  return ret;
}


HalfOpenConeList refinement(HalfOpenConeList const &a, HalfOpenConeList const &b)
{
  HalfOpenConeList ret;
  for(HalfOpenConeList::const_iterator i=a.begin();i!=a.end();i++)
    for(HalfOpenConeList::const_iterator j=b.begin();j!=b.end();j++)
      if(!haveEmptyIntersection(*i,*j))
      {
	HalfOpenCone c=intersection(*i,*j);
	//	c.isEmpty();
	//	c.isEmpty();
	//	if(!c.isEmpty())
	  ret.push_back(c);
      }

  return ret;
}


HalfOpenConeList tropicalHyperSurfaceIntersection2(int dimension, PolynomialSet const &g)
{
  HalfOpenConeList intersection;
  IntegerVectorList empty;
  intersection.push_back(HalfOpenCone(dimension,empty,empty,empty));

  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    {
      HalfOpenConeList surface=tropicalHyperSurface(*i);

      fprintf(Stderr,"Number of cones in current intersection:%i\n",intersection.size());
      fprintf(Stderr,"Number of cones in next surface:%i\n",surface.size());

      fprintf(Stderr,"A\n");
      intersection=refinement(intersection,surface);
      fprintf(Stderr,"B\n");
    }
  fprintf(Stderr,"%i",intersection.size());

  return intersection;
}


void tropicalHyperSurfaceIntersectionWithProcessor(int dimension, PolynomialSet const &g, HalfOpenConeProcessor &myProcessor, PolyhedralCone *restrictingCone, bool expand, int intervalLow, int intervalHigh)
{
  assert(expand || !restrictingCone);//Need to fix this if we don't want to expand. Add in a different HalfOpenConeProcessor.
  //  HalfOpenConeList intersection;
  if(restrictingCone)
    {
      IntegerVectorList equations=restrictingCone->getEquations();
      IntegerVectorList nonStrict=restrictingCone->getHalfSpaces();
      IntegerVectorList strict;

      /****************************************************************
       * REMOVE THE LINE BELOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       */
//      strict.push_back(IntegerVector::standardVector(g.getRing().getNumberOfVariables(),0));cerr<<"REMOVE THIS LINE\n";

      HalfOpenCone restrictingCone2(dimension,equations,nonStrict,strict,true);

      //intersection=
      tropicalHyperSurfaceIntersectionInSubspace(dimension,g,&restrictingCone2,myProcessor,intervalLow,intervalHigh);//<<---here
      //intersection=tropicalHyperSurfaceIntersection(dimension,g,&restrictingCone2);
    }
    else
      {
	HalfOpenCone rCone(dimension);
	//intersection=
	tropicalHyperSurfaceIntersectionInSubspace(dimension,g,&rCone,myProcessor,intervalLow,intervalHigh);//<<---here
	/*    PolyhedralFan intersectionOld=tropicalHyperSurfaceIntersection(dimension,g);
      AsciiPrinter P(Stderr);
      intersection.print(&P);
      intersectionOld.print(&P);
	*/
      }
  log2 fprintf(Stderr,"Halfopen intersection completed.\n");
}

PolyhedralFan tropicalHyperSurfaceIntersectionClosed(int dimension, PolynomialSet const &g, PolyhedralCone *restrictingCone, bool expand, bool saveResult, int intervalLow, int intervalHigh)
{
	  HalfOpenConeProcessorPointCollector myProcessor;
	  if(saveResult)myProcessor.setSave();

	  tropicalHyperSurfaceIntersectionWithProcessor(dimension,g,myProcessor,restrictingCone,expand,intervalLow,intervalHigh);

	  if(saveResult)saveList(myProcessor.interiorPoints,"finalresult");


	  //  AsciiPrinter P(Stderr);
	  //  printHalfOpenConeList(intersection, P);

	  PolyhedralFan ret(dimension);
	  //  for(HalfOpenConeList::iterator i=intersection.begin();i!=intersection.end();i++)
	  for(IntegerVectorList::const_iterator i=myProcessor.interiorPoints.begin();i!=myProcessor.interiorPoints.end();i++)
	    {
	      /*PolyhedralCone c=i->closure();
	      c.canonicalize();

	      if(expand)
		c=normalConeOfMinkowskiSum(g,c.getRelativeInteriorPoint());

	      ret.insert(c);
	      */
	      ret.insert(normalConeOfMinkowskiSum(g,*i));
	    }
	  /*      AsciiPrinter P(Stderr);
	      ret.print(&P);
	      cerr<<"ETSTSETST------------------------------------------------------------------------------";
	  */

	  return ret;

}
void HalfOpenCone::splitIntoRelativelyOpenCones(list<HalfOpenCone> &l)
{
  //  fprintf(Stderr,"BEGIN\n");
  //  AsciiPrinter P(Stderr);
  //  print(P);
  lifted.findFacets();
  //  print(P);

  /*	{
	  IntegerVector v=StringParser("(3,0,3,2,0,3)").parseIntegerVector();
	  if(contains(v))fprintf(Stderr,"??????????????????????????????????????????\n");
	  }*/

  IntegerVectorList inequalityList=lifted.getHalfSpaces();
  IntegerVectorList equationList=lifted.getEquations();

  IntegerVectorList strict,nonstrict;


  for(IntegerVectorList::const_iterator i=inequalityList.begin();i!=inequalityList.end();i++)
    if((*i)[i->size()-1]<0)
      strict.push_back(*i);
    else if((*i)[i->size()-1]==0)
      nonstrict.push_back(*i);
    else
      {//CHANGED
	assert(i->subvector(0,i->size()-1).isZero());
	strict.push_back(*i);
      }

  //  AsciiPrinter(Stderr).printVectorList(nonstrict);
  //  AsciiPrinter(Stderr).printVectorList(strict);
  //  AsciiPrinter(Stderr).printVectorList(equationList);

  if(nonstrict.size()==0)
    {
      l.push_back(*this);
    }
  else
    {
      IntegerVector chosen=*nonstrict.begin();
      nonstrict.pop_front();

      strict.push_front(chosen);
      (*strict.begin())[strict.begin()->size()-1]=-1;
      IntegerVectorList a=nonstrict;
      IntegerVectorList tempa=strict;
      a.splice(a.begin(),tempa);

      //  fprintf(Stderr,"New inequalities:\n");
      //  AsciiPrinter(Stderr).printVectorList(a);

      HalfOpenCone A(dimension,PolyhedralCone(a,equationList,liftedDimension));
      A.splitIntoRelativelyOpenCones(l);


      equationList.push_front(chosen);
      strict.pop_front();

      IntegerVectorList b=nonstrict;
      IntegerVectorList tempb=strict;
      b.splice(b.begin(),tempb);
      // fprintf(Stderr,"New inequalities:\n");
      // AsciiPrinter(Stderr).printVectorList(b);
      // fprintf(Stderr,"New equationList:\n");
      // AsciiPrinter(Stderr).printVectorList(equationList);

      HalfOpenCone B(dimension,PolyhedralCone(b,equationList,liftedDimension));
      B.splitIntoRelativelyOpenCones(l);
    }
      //  AsciiPrinter(Stderr).print
  //  fprintf(Stderr,"END\n");
}


void HalfOpenCone::print(class Printer &p)const
{
  p.printString("Printing HalfOpenCone\n");
  lifted.print(&p);
  p.printString("Done printing HalfOpenCone\n");
}

HalfOpenConeList splitIntoRelativelyOpenCones(HalfOpenConeList const &l)
{
  AsciiPrinter P(Stderr);
  HalfOpenConeList ret;
  for(HalfOpenConeList::const_iterator i=l.begin();i!=l.end();i++)
    {
      fprintf(Stderr,"A");
      HalfOpenCone temp=*i;
      HalfOpenConeList tempSplit;
      //  fprintf(Stderr,"---------------------------------------------------------------\n");
      //  temp.print(P);
      //  fprintf(Stderr,"---------------------------------------------------------------\n");
      temp.splitIntoRelativelyOpenCones(tempSplit);

      //  fprintf(Stderr,"Splits into:");
      //  for(HalfOpenConeList::const_iterator i=tempSplit.begin();i!=tempSplit.end();i++)
      //    i->print(P);
      //  fprintf(Stderr,"Splits into End.");

      ret.splice(ret.begin(),tempSplit);

      fprintf(Stderr,"B\n");
    }

  return ret;
}

static bool isSubsetOf(IntegerVector const &v, IntegerVector const &u)
{
  for(int i=0;i<v.size();i++)
    {
      bool found=false;
      for(int j=0;j<u.size();j++)
	if(v[i]==u[j])
	  {
	    found=true;
	    break;
	  }
      if(!found)return false;
    }
  return true;
}

static bool isSubsetOf(IntegerVector const &v, IntegerVectorList const &l)
{
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    if(isSubsetOf(v,*i))return true;

  return false;
}


void printHalfOpenConeList(HalfOpenConeList const &l, class Printer & p)
{
  HalfOpenConeList L=splitIntoRelativelyOpenCones(l);

  list<PolyhedralCone> cones;
  for(HalfOpenConeList::iterator i=L.begin();i!=L.end();i++)
    cones.push_back(i->closure());

  int homog=1000000;
  int largest=0;
  int ambientDimension=-1;
  for(list<PolyhedralCone>::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      if(i->dimension()<homog)homog=i->dimension();
      if(i->dimension()>largest)largest=i->dimension();
      ambientDimension=i->ambientDimension();
    }
  assert(homog!=1000000);

  for(list<PolyhedralCone>::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      assert(i->dimensionOfLinealitySpace()==homog);
    }
  fprintf(Stderr,"Ambient dimension: %i, maximal dimension: %i, dimension of lineality space: %i\n",ambientDimension,largest,homog);

  IntegerVectorList rays;
  for(list<PolyhedralCone>::iterator i=cones.begin();i!=cones.end();i++)
    if(i->dimension()==homog+1)rays.push_back(i->getRelativeInteriorPoint());

  p.printString("Rays:\n");

  p.printVectorList(rays,true);

  list<IntegerVectorList> subsets;
  for(int d=homog;d<=largest;d++)
    {
      IntegerVectorList thisDimension;
      list<PolyhedralCone> cones2;
      for(list<PolyhedralCone>::iterator i=cones.begin();i!=cones.end();i++)
	if(i->dimension()==d)
	  cones2.push_back(*i);

      for(list<PolyhedralCone>::const_iterator i=cones2.begin();i!=cones2.end();i++)
	{
	  IntegerVector v(0);
	  int J=0;
	  for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++)
	    {
	      if(i->contains(*j))
		{
		  v.grow(v.size()+1);
		  v[v.size()-1]=J;
		}
	      J++;
	    }
	  thisDimension.push_back(v);
	}
      subsets.push_back(thisDimension);
    }


  list<IntegerVectorList>::const_iterator subsetIterator=subsets.begin();
  list<IntegerVectorList>::const_iterator subsetIteratorNext=subsets.begin();
  for(int d=homog;d<=largest;d++)
    {
      subsetIteratorNext++;
      IntegerVectorList maximal,nonmaximal;

      if(subsetIteratorNext!=subsets.end())
	{
	  for(IntegerVectorList::const_iterator i=subsetIterator->begin();i!=subsetIterator->end();i++)
	    if(isSubsetOf(*i,*subsetIteratorNext))
	      nonmaximal.push_back(*i);
	    else
	      maximal.push_back(*i);
	}
      else
	maximal=*subsetIterator;

      p.printString("Printing ");p.printInteger(subsetIterator->size());p.printString(" ");p.printInteger(d);p.printString("-dimensional cones (");p.printInteger(maximal.size());p.printString(" maximal cones):\n");
      p.printString("{");
      {
	bool first=true;
	for(IntegerVectorList::const_iterator i=maximal.begin();i!=maximal.end();i++)
	  {
	    if(!first)p.printString(",\n");
	    p.printVector(*i);
	    first=false;
	  }
	if(!first && nonmaximal.size()!=0)p.printString(",\n");
	p.printString("\n");
	first=true;
	for(IntegerVectorList::const_iterator i=nonmaximal.begin();i!=nonmaximal.end();i++)
	  {
	    if(!first)
	      p.printString(",\n");
	    p.printVector(*i);
	    first=false;
	  }
      }
      p.printString("}\n");
      subsetIterator++;
   }
  /*  for(int d=homog;d<=largest;d++)
    {
      list<PolyhedralCone> cones2;
      for(list<PolyhedralCone>::iterator i=cones.begin();i!=cones.end();i++)
	if(i->dimension()==d)
	  cones2.push_back(*i);

      p.printString("Printing ");p.printInteger(cones2.size());p.printString(" ");p.printInteger(d);p.printString("-dimensional cones:\n");
      p.printString("{");
      for(list<PolyhedralCone>::const_iterator i=cones2.begin();i!=cones2.end();i++)
	{
	  IntegerVector v(0);
	  int J=0;
	  for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++)
	    {
	      if(i->contains(*j))
		{
		  v.grow(v.size()+1);
		  v[v.size()-1]=J;
		}
	      J++;
	    }
	  if(i!=cones2.begin())p.printString(",\n");
	  p.printVector(v);
	}
      p.printString("}\n");
    }
*/
}


class BitSet
{
  vector<int> v;
public:
  BitSet()
  {
  }
  BitSet(int n):
    v(n)
  {
    for(int i=0;i<n;i++)v[i]=0;
  }
  int& operator[](int n){assert(n>=0 && n<v.size());return (v[n]);}
  const int& operator[](int n)const{assert(n>=0 && n<v.size());return (v[n]);}
  void add(BitSet const &b)
  {
    assert(b.size()==v.size());
    for(int i=0;i<b.size();i++)
      {
	//	fprintf(Stderr,"%i\n",i);
	v[i]=v[i]||b[i];
      }
  }
  int size()const
  {
    return v.size();
  }
  void print(Printer &P)const
  {
    P.printString("(");
    for(int i=0;i<v.size();i++)
      {
	if(i!=0)P.printString(", ");
	P.printInteger(v[i]);
      }
    P.printString(")\n");
  }
  BitSet negated()const
  {
    BitSet ret(size());
    for(int i=0;i<size();i++)ret[i]=1-v[i];
    return ret;
  }
  int sizeOfSubset()const
  {
    int ret=0;
    for(int i=0;i<size();i++)if(v[i])ret++;
    return ret;
  }
  };
/*
class BitSet
{
#define BITSET_BITS_PER_WORD 32
  vector<int> v;
  int ambientSetSize;
  static int bitsSetInInt(int i)
  {
    if(BITSET_BITS_PER_WORD==32)
      {
	i=(i&0x55555555)+((i>>1)&0x55555555);
	i=(i&0x33333333)+((i>>2)&0x33333333);
	i=(i&0x0f0f0f0f)+((i>>4)&0x0f0f0f0f);
	i=(i)+((i>>8));
	i=(i)+((i>>16));
	return i&255;
      }


    int ret=0;
    while(i&(i-1)){i&=i-1;ret++;}
    return ret;
  }
public:
  BitSet()
  {
  }
  BitSet(int n):
    ambientSetSize(n),
    v(((unsigned int)(ambientSetSize-1))/BITSET_BITS_PER_WORD+1)
  {
    for(int i=0;i<v.size();i++)v[i]=0;
  }
  //  int& operator[](unsigned int n){assert(n>=0 && n<ambientSize);return (v[n/BITSET_BITS_PER_WORD]&&);}
  inline bool operator[](unsigned int n)const{assert(n>=0 && n<ambientSetSize);return (v[n/BITSET_BITS_PER_WORD]&(((int)1)<<(n&(BITSET_BITS_PER_WORD-1))));}
  inline void setValue(unsigned int index, bool value)
  {
    if(value)
      v[index/BITSET_BITS_PER_WORD]|=((int)1)<<(index&(BITSET_BITS_PER_WORD-1));
    else
      v[index/BITSET_BITS_PER_WORD]&=-1-(((int)1)<<(index&(BITSET_BITS_PER_WORD-1)));
  }
  void add(BitSet const &b)
  {
    assert(b.size()==ambientSetSize);
    for(int i=0;i<ambientSetSize;i++)v[i]|=b.v[i];
  }
  int size()const
  {
    return ambientSetSize;
  }
  void print(Printer &P)const
  {
    P.printString("(");
    for(int i=0;i<ambientSetSize;i++)
      {
	if(i!=0)P.printString(", ");
	P.printInteger((*this)[i]);
      }
    P.printString(")\n");
  }
  BitSet negated()const
  {
    BitSet ret(ambientSetSize);
    for(int i=0;i<v.size();i++)ret.v[i]=-1-v[i];
    return ret;
  }
  int sizeOfSubset()const
  {
    int ret=0;
    for(int i=0;i<v.size();i++)ret+=bitsSetInInt(v[i]);
    return ret;
  }
  int sizeOfComplement()const
  {
    return ambientSetSize-sizeOfSubset();
  }
};
*/
class Table
{
  vector<vector<vector<BitSet> > > table;
public:
  Table(vector<vector<HalfOpenCone> > const &l):
    table(l.size())
  {
    int N=l.size();
    for(int i=0;i<N;i++)
      {
	vector<vector<BitSet> > v(N);
	for(int j=0;j<N;j++)
	  {
	    vector<BitSet> w(l[i].size());
	    for(int k=0;k<l[i].size();k++)
	      w[k]=BitSet(l[j].size());
	    v[j]=w;
	  }
	table[i]=v;
      }
  }
  bool lookUp(int fan1, int cone1, int fan2, int cone2)
  {
    assert(fan1<table.size());
    assert(fan2<table[fan1].size());
    assert(cone1<table[fan1][fan2].size());
    assert(cone2<table[fan1][fan2][cone1].size());

    return table[fan1][fan2][cone1][cone2];
  }
  void set(int fan1, int cone1, int fan2, int cone2)
  {
    assert(fan1<table.size());
    assert(fan2<table[fan1].size());
    assert(cone1<table[fan1][fan2].size());
    assert(cone2<table[fan1][fan2][cone1].size());

    table[fan1][fan2][cone1][cone2]=true;
    table[fan2][fan1][cone2][cone1]=true;
    //table[fan1][fan2][cone1].setValue(cone2,true);
    //table[fan2][fan1][cone2].setValue(cone1,true);
  }
  BitSet const& nonCandidates(int fan1, int cone1, int fan2)const
  {
    assert(fan1<table.size());
    assert(fan2<table[fan1].size());
    assert(cone1<table[fan1][fan2].size());

    return table[fan1][fan2][cone1];
  }
  void print()const
  {
    AsciiPrinter P(Stderr);
    for(int i=0;i<table.size();i++)
      for(int j=0;j<table[i].size();j++)
	{
	  fprintf(Stderr,"Entry (%i,%i)\n",i,j);
	  for(int k=0;k<table[i][j].size();k++)
	    table[i][j][k].print(P);
	}
  }
};

class RelationTable
{
  vector<vector<HalfOpenCone> > fanList;
  Table knownEmptyIntersectionInIntersection;
  Table knownNonEmptyIntersection;
public:
  int numberOfSolvedLPs;
  RelationTable(vector<vector<HalfOpenCone> > const &l):
    fanList(l),
    knownEmptyIntersectionInIntersection(l),
    knownNonEmptyIntersection(l),
    numberOfSolvedLPs(0)
  {

  }
  bool knownToIntersectTriviallyInIntersection(int fan1, int cone1, int fan2, int cone2)
  {
    assert(fan1<fanList.size());
    assert(fan2<fanList.size());
    assert(cone1<fanList[fan1].size());
    assert(cone2<fanList[fan2].size());

    return knownEmptyIntersectionInIntersection.lookUp(fan1,cone1,fan2,cone2);
  }
  bool intersectTriviallyInIntersection(int fan1, int cone1, int fan2, int cone2)
  {
    assert(fan1<fanList.size());
    assert(fan2<fanList.size());
    assert(cone1<fanList[fan1].size());
    assert(cone2<fanList[fan2].size());

    if(knownEmptyIntersectionInIntersection.lookUp(fan1,cone1,fan2,cone2))
      return true;
    if(knownNonEmptyIntersection.lookUp(fan1,cone1,fan2,cone2))
      return false;

    //    fprintf(Stderr,"UPDATING:f1:%i,c1:%i,f2:%i,c2:%i\n",fan1,cone1,fan2,cone2);
    bool ret=haveEmptyIntersection(fanList[fan1][cone1],fanList[fan2][cone2]);
    numberOfSolvedLPs++;
    if(ret)
      knownEmptyIntersectionInIntersection.set(fan1,cone1,fan2,cone2);
    else
      knownNonEmptyIntersection.set(fan1,cone1,fan2,cone2);
    return ret;
  }
  const BitSet &getNonCandidates(int fan1, int cone1, int fan2)
  {
    //    for(int c2=0;c2<fanList[fan2].size();c2++)
    //      intersectTriviallyInIntersection(fan1,cone1,fan2,c2);

    return knownEmptyIntersectionInIntersection.nonCandidates(fan1,cone1,fan2);
  }
  void markNoIntersectionInIntersection(int fan1, int cone1, int fan2, int cone2)
  {
    knownEmptyIntersectionInIntersection.set(fan1,cone1,fan2,cone2);
  }
  void markKnownNonEmptyIntersection(int fan1, int cone1, int fan2, int cone2)
  {
	  knownNonEmptyIntersection.set(fan1,cone1,fan2,cone2);
  }
  void print()const
  {
    fprintf(Stderr,"knownEmptyIntersectionInIntersection:");
    knownEmptyIntersectionInIntersection.print();
    fprintf(Stderr,"knownNonEmptyIntersection:");
    knownNonEmptyIntersection.print();
  }
};



struct RecursionData
{
  vector<vector<HalfOpenCone> > fans;
  IntegerVector chosen;
  IntegerVector chosenFans;
  IntegerVector iterators; //just used for printing
  IntegerVector nCandidates; //just used for printing
  BitSet usedFans;
  int numberOfUsefulCalls;
  int totalNumberOfCalls;
  HalfOpenConeProcessor &processor;
  bool intervalSet;
  int intervalLow,intervalHigh;
public:
  RelationTable table;
  RecursionData(vector<vector<HalfOpenCone> > const &fans_, HalfOpenConeProcessor &processor_):
    table(fans_),
    fans(fans_),
    chosen(fans_.size()),
    chosenFans(fans_.size()),
    usedFans(fans_.size()),
    iterators(fans_.size()),
    nCandidates(fans_.size()),
    numberOfUsefulCalls(0),
    totalNumberOfCalls(0),
    processor(processor_),
    intervalSet(false)
  {
  }
  void setInterval(int lo, int hi)
  {
    intervalLow=lo;
    intervalHigh=hi;
    intervalSet=true;
  }
  //  HalfOpenConeList ret;

  bool randBool()const
  {
    return 0;
    static int i;
    i++;
    return (i&3)==0;
  }

  BitSet computeCandidates(int index, int fanNumber)
  {
    /*    BitSet nonCandidates(fans[fanNumber].size());
    for(int i=0;i<index;i++)
      {
	nonCandidates.add(table.getNonCandidates(chosenFans[i],chosen[i],fanNumber));
      }
    return nonCandidates.negated();
    */

    BitSet nonCandidates(fans[fanNumber].size());
    for(int i=0;i<index;i++)
      {
	nonCandidates.add(table.getNonCandidates(chosenFans[i],chosen[i],fanNumber));
      }
    for(int j=0;j<nonCandidates.size();j++)
      if((!nonCandidates[j]) ||randBool())
	for(int i=0;i<index;i++)
	  if(table.intersectTriviallyInIntersection(chosenFans[i], chosen[i], fanNumber, j))
	    {
	      nonCandidates[j]=true;
	      //	      nonCandidates.setValue(j,true);
	      break;
	    }

    return nonCandidates.negated();
  }
  bool rek(int index, HalfOpenCone const &current)
  {
    totalNumberOfCalls++;

    //    if((totalNumberOfCalls&15) == 0)closure();
    //    if((totalNumberOfCalls&255) == 255)closure();

    bool success=false;

    if(index == fans.size())
      {
	log2 fprintf(Stderr,"ADDING CONE\n");
	//ret.push_back(current);
	processor.process(current);
	numberOfUsefulCalls++;
	return true;
      }
    else
      {
	AsciiPrinter P(Stderr);

	int bestIndex=-1;
	int bestNumberOfCandidates=1000000;
	for(int i=0;i<fans.size();i++)
	  {
	    if(!usedFans[i])
	      {
		int n=computeCandidates(index,i).sizeOfSubset();
		//		fprintf(Stderr,"Number of candidates for fan %i: %i\n",i,n);
		if(n<=bestNumberOfCandidates)  //we could choose a strict inequality
		  {
		    bestNumberOfCandidates=n;
		    bestIndex=i;
		  }
	      }
	  }
	assert(bestIndex!=-1);
	BitSet candidates=computeCandidates(index,bestIndex);


	chosenFans[index]=bestIndex;
	usedFans[chosenFans[index]]=true;
	//usedFans.setValue(chosenFans[index],true);


	nCandidates[index]=bestNumberOfCandidates;//just for printing

	static int iterationNumber;
	if(!(iterationNumber++ & 31))
	  log2
	{
	  fprintf(Stderr,"Iteration level:%i, Chosen fan:%i, Number of candidates:%i, Iteration Number:%i, Useful (%i/%i)=%f\n",index,bestIndex,bestNumberOfCandidates,iterationNumber,numberOfUsefulCalls,totalNumberOfCalls,float(numberOfUsefulCalls)/totalNumberOfCalls);
	  fprintf(Stderr,"Chosen fans vector: ");
	  P.printVector(chosenFans,false,2);
	  fprintf(Stderr,"\nChosen cone vector: ");
	  P.printVector(chosen,false,2);
	  fprintf(Stderr,"\nNcandidates vector: ");
	  P.printVector(nCandidates,false,2);
	  fprintf(Stderr,"\nIterator vector:    ");
	  P.printVector(iterators,false,2);
	  fprintf(Stderr,"\n\n");
	}

	//
	/*	if(index>1)
	  {
	    smallest=1000000;
	    bx=-1;
	    by=-1;
	    for(int x=0;x<index;x++)
	      for(int y=x;y<index;y++)
		{
		  BitSet c1=table.getNonCandidates(chosenFans[x],chosen[x],bestIndex);
		  c1.add(table.getNonCandidates(chosenFans[y],chosen[y],bestIndex));
		  int n=c1.sizeOfSubset();
		  if(n<smallest)
		    {
		      bx=x;
		      by=y;
		    }
		}
	    bool skipThisOne=true;
	    for(int i=0;i<fans[chosenFans[index]].size();i++)
	      {
		if(!table.intersectTriviallyInIntersection(chosenFans[x],chosen[x],chosenFans[index],i))skipThisOne=false;
		if(!table.intersectTriviallyInIntersection(chosenFans[y],chosen[y],chosenFans[index],i))skipThisOne=false;
	      }
	      }*/


	// P.printInteger(fans[index].size());
	for(int i=0;i<fans[chosenFans[index]].size();i++)
	  if(candidates[i])
	    {
	      if(intervalSet && (index==0))
		{
		  if(i<intervalLow || i>=intervalHigh)
		    {
		      cerr << "SKIPPING CONE "<<i<<endl;
		      continue;
		    }
		}
	      bool ok=true;
	      for(int j=0;j<index;j++)
		{
		  if(table.intersectTriviallyInIntersection(chosenFans[j],chosen[j],chosenFans[index],i))
		    {
		      ok=false;
		      break;

		    }
		}
	      if(ok && !haveEmptyIntersection(current,fans[chosenFans[index]][i]))
		{
		  chosen[index]=i;

		  //log0 fprintf(Stderr,"A\n");
		  HalfOpenCone next=intersection(current,fans[chosenFans[index]][i],true);
		  //log0 fprintf(Stderr,"B\n");
		  //if((fans.size()!=10))//&&(fans.size()!=9))
		    {
		      bool s=rek(index+1,next);
		      success|=s;
		    }
		  /*		  else //Some old optimization code if we have 10 fans????
		    {
		      vector<vector<HalfOpenCone> > L2;
		      vector<int> indicesOfNonUsedFans;
		      vector<vector<int> > indicesOfCones;
		      {
			for(int i=0;i<fans.size();i++)
			  if(!usedFans[i])
			    {
			      vector<HalfOpenCone> L;
			      indicesOfCones.push_back(vector<int>());
			      for(int j=0;j<fans[i].size();j++)
				{
				  if(!haveEmptyIntersection(next,fans[i][j]))
				    {
				      L.push_back(intersection(next,fans[i][j],false));
				      indicesOfCones.back().push_back(j);
				    }
				}
			      fprintf(Stderr,"New fan size:%i\n",L.size());
			      L2.push_back(L);
			      indicesOfNonUsedFans.push_back(i);
			    }
		      }
		      RecursionData data(L2);
		      for(int f1=0;f1<L2.size();f1++)
			for(int f2=f1+1;f2<L2.size();f2++)
			  for(int c1=0;c1<L2[f1].size();c1++)
			    for(int c2=0;c2<L2[f2].size();c2++)
			      if(table.knownToIntersectTriviallyInIntersection
				 (indicesOfNonUsedFans[f1],indicesOfCones[f1][c1],
				  indicesOfNonUsedFans[f2],indicesOfCones[f2][c2]))
				data.table.markNoIntersectionInIntersection(f1, c1, f2, c2);

		      //data.completeTable();
		      data.transitiveClosure();
		      IntegerVectorList empty;
		      success|=data.rek(0, HalfOpenCone(next.dimension,empty,empty,empty));
		      ret.splice(ret.begin(),data.ret);
		    }*/
		  chosen[index]=-1;//just for printing
		}
	      iterators[index]++;//just for printing
	    }
	/*	if(!success)
	  {
	    for(int x=0;x<index;x++)
	      for(int y=x;y<index;y++)
		{
		  BitSet c1=table.getNonCandidates(chosenFans[x],chosen[x],bestIndex);
		  c1.add(table.getNonCandidates(chosenFans[y],chosen[y],bestIndex));
		  int n=c1.negated().sizeOfSubset();
		  if(n==0)
		    {
		      table.markNoIntersectionInIntersection(chosenFans[x],chosen[x],chosenFans[y],chosen[y]);
		      fprintf(Stderr,"ONE FOR FREE\n");
		    }
		}
		}*/

	nCandidates[index]=-1;//just for printing
	iterators[index]=0;//just for printing

	usedFans[chosenFans[index]]=false;
	//	usedFans.setValue(chosenFans[index],false);
	chosenFans[index]=-1;
      }
    if(success)numberOfUsefulCalls++;
    return success;
  }
  bool closure()
  {
    log2 cerr<<"computing closure"<<endl;
    bool ret=false;
    int a=0;
    for(int f1=0;f1<fans.size();f1++)
      {
	//	fprintf(Stderr,"%i\n",f1);
	for(int f2=f1+1;f2<fans.size();f2++)
	  for(int c1=0;c1<fans[f1].size();c1++)
	    for(int c2=0;c2<fans[f2].size();c2++)
	      {
		//		if(!table.intersectTriviallyInIntersection(f1,c1,f2,c2))
		if(!table.knownToIntersectTriviallyInIntersection(f1,c1,f2,c2))
		  {
		    bool dontintersect=false;
		    for(int f3=0;f3<fans.size();f3++)
		      {
			BitSet c=table.getNonCandidates(f1,c1,f3);
			c.add(table.getNonCandidates(f2,c2,f3));
			if(c.negated().sizeOfSubset()==0)
			  {
			    dontintersect=true;
			    a++;
			    break;
			  }
			if(c.negated().sizeOfSubset()<4 && ((f3&7) ==0))//just an experiment
			  {
			    for(int k=0;k<c.size();k++)
			      if(!c[k])
				{
				  table.intersectTriviallyInIntersection(f1,c1,f3,k);
				  table.intersectTriviallyInIntersection(f2,c2,f3,k);
				}
			  }
		      }
		    if(dontintersect)
		      {
			table.markNoIntersectionInIntersection(f1,c1,f2,c2);
			ret=true;
		      }
		  }
	      }
      }
    log2 fprintf(Stderr,"%i FOR FREE\n",a);
    log2 cerr<<"done computing closure"<<endl;
    return ret;
  }

  void transitiveClosure()
  {
    while(closure());
  }

  void completeTable()
  {
	    for(int f1=0;f1<fans.size();f1++)
	    	for(int c1=0;c1<fans[f1].size();c1++)
	    	  for(int c2=0;c2<fans[f1].size();c2++)
	    		  if(c1!=c2)
	    			  table.markNoIntersectionInIntersection(f1,c1,f1,c2);
//	    		  else
//	    			  table.markKnownNonEmptyIntersection(f1,c1,f1,c2);


	    	for(int f1=0;f1<fans.size();f1++)
      {
	log3 fprintf(Stderr,"%i\n",f1);
	for(int f2=f1+1;f2<fans.size();f2++)
	for(int c1=0;c1<fans[f1].size();c1++)
	  for(int c2=0;c2<fans[f2].size();c2++)
	    table.intersectTriviallyInIntersection(f1,c1,f2,c2);
      }
    transitiveClosure();
//    table.print();//HERE
  }
  void extractInformationFromFullFan(vector<vector<HalfOpenCone> > const &fans2, vector<IntegerVector> parentList)
  {
    int freeLPs=0;
    HalfOpenConeProcessorVoid dummy;
    RecursionData data2(fans2,dummy);
    data2.completeTable();
    data2.transitiveClosure();

    assert(fans.size()==fans2.size());
    assert(fans.size()==parentList.size());
    //    data2.table.print();

    for(int f1=0;f1<fans.size();f1++)
      for(int f2=f1+1;f2<fans.size();f2++)
	for(int c1=0;c1<fans[f1].size();c1++)
	  for(int c2=0;c2<fans[f2].size();c2++)
	    if(data2.table.knownToIntersectTriviallyInIntersection
	       (f1,parentList[f1][c1],
		f2,parentList[f2][c2]))
	      {
		AsciiPrinter P(Stderr);

		/*		if(!table.intersectTriviallyInIntersection(f1,c1,f2,c2))
		  {
		    fprintf(Stderr,"f1: %i c1: %i f2: %i c2: %i p[f1][c1]: %i p[f2][c2]: %i\n",f1,c1,f2,c2,parentList[f1][c1],parentList[f2][c2]);
		    HalfOpenCone A2=fans2[f1][parentList[f1][c1]];
		    HalfOpenCone B2=fans2[f2][parentList[f2][c2]];
		    HalfOpenCone A1=fans[f1][c1];
		    HalfOpenCone B1=fans[f2][c2];

		    fprintf(Stderr,"%i %i\n",haveEmptyIntersection(A1,B1),haveEmptyIntersection(A2,B2));
		    A1.print(P);
		    B1.print(P);
		    A2.print(P);
		    B2.print(P);
		    }*/


		//		assert(table.intersectTriviallyInIntersection(f1,c1,f2,c2));
		table.markNoIntersectionInIntersection(f1, c1, f2, c2);
		freeLPs++;
	      }
    //    table.print();
    log2 fprintf(Stderr,"Number of infeassible LPs discovered from higherdimensional cones: %i\n",freeLPs);
  }
};



void tropicalHyperSurfaceIntersection(int dimension, PolynomialSet const &g, HalfOpenCone *restrictingCone, HalfOpenConeProcessor &processor)
{
  vector<vector<HalfOpenCone> > fullDimFanList;
  vector<vector<HalfOpenCone> > coDimOneFanList;
  vector<IntegerVector> parents(g.size());

  //  AsciiPrinter(Stderr) << g;

  //  fprintf(Stderr,"GSIZE:%i\n",g.size());
  {
    int I=0;
    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
      {
	HalfOpenConeList l,lf;
	buildFanFromPolynomial(*i, &lf, &l, &(parents[I]),restrictingCone);

	log2 AsciiPrinter(Stderr).printVector(parents[I]);

	log2 fprintf(Stderr,"\n");

	vector<HalfOpenCone> L;
	for(HalfOpenConeList::const_iterator i=l.begin();i!=l.end();i++)
	  {
	    L.push_back(*i);
	  }
	vector<HalfOpenCone> F;
	for(HalfOpenConeList::const_iterator i=lf.begin();i!=lf.end();i++)
	  {
	    F.push_back(*i);
	  }
	fullDimFanList.push_back(F);
	coDimOneFanList.push_back(L);
	//if(b==4)break;
	I++;
      }
  }

  RecursionData data(coDimOneFanList,processor);
  data.extractInformationFromFullFan(fullDimFanList,parents);

//  data.completeTable();//HERE
//  data.table.print();//HERE
  IntegerVectorList empty;
  data.rek(0, HalfOpenCone(dimension,empty,empty,empty));
  log2 fprintf(Stderr,"LPs solved:%i for relation table\n",data.table.numberOfSolvedLPs);
}


HalfOpenConeList tropicalHyperSurfaceIntersection(int dimension, PolynomialSet const &g, HalfOpenCone *restrictingCone)
{
  HalfOpenConeProcessorConeCollector collector;
  tropicalHyperSurfaceIntersection(dimension, g, restrictingCone, collector);
  return collector.theList;
}


#include "binomial.h"
#include "linalg.h"

//ignore last coordinate. Get rid of pivots
static IntegerVectorList rewriteVectorList(IntegerVectorList const &l, list<int> nonPivots, FieldMatrix const &A)
{
  IntegerVectorList ret;
  int M=nonPivots.size();

  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      IntegerVector v(M+1);
      v[M]=(*i)[i->size()-1];
      IntegerVector u=A.canonicalize(integerVectorToFieldVector(i->subvector(0,i->size()-1),Q)).primitive();
      int J=0;
      for(list<int>::const_iterator j=nonPivots.begin();j!=nonPivots.end();j++,J++)v[J]=u[*j];
      ret.push_back(v);
    }

  return ret;
}

//add pivot columns
static IntegerVectorList expandVectorList(IntegerVectorList const &l, list<int> pivots)
{
  IntegerVectorList ret;
  int pivotSize=pivots.size();

  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      int newSize=pivotSize+i->size();

      IntegerVector v(newSize);

      list<int>::const_iterator J1=pivots.begin();
      int J2=0;
      for(int j=0;j<newSize;j++)
	{
	  if((J1!=pivots.end())&&(((*J1))==j))
	    {
	      v[j]=0;
	      J1++;
	    }
	  else
	    {
	      v[j]=(*i)[J2];
	      J2++;
	    }
	}

      ret.push_back(v);
    }

  return ret;
}


HalfOpenCone HalfOpenCone::withChosenCoordinates(list<int> chosen)const
{
  chosen.push_back(dimension);
  PolyhedralCone newLifted=PolyhedralCone(
					  subvectorsOfIntegerVectorList(lifted.getHalfSpaces(),chosen),
					  subvectorsOfIntegerVectorList(lifted.getEquations(),chosen),
					  chosen.size()
					  );
  return HalfOpenCone(chosen.size()-1,newLifted);
}


HalfOpenCone HalfOpenCone::rewrite(FieldMatrix const &A, list<int> nonPivots)const
{
  PolyhedralCone newLifted=PolyhedralCone(
					  rewriteVectorList(lifted.getHalfSpaces(),nonPivots,A),
					  rewriteVectorList(lifted.getEquations(),nonPivots,A),
					  nonPivots.size()+1
					  );
  return HalfOpenCone(nonPivots.size(),newLifted);
}


HalfOpenCone HalfOpenCone::rewriteExpand(list<int> pivots, IntegerVectorList const &newEquations)const
{
  IntegerVectorList equations=expandVectorList(lifted.getEquations(),pivots);

  for(IntegerVectorList::const_iterator i=newEquations.begin();i!=newEquations.end();i++)
    {
      IntegerVector v(i->size()+1);
      for(int j=0;j<v.size()-1;j++)v[j]=(*i)[j];
      equations.push_back(v);
    }

  PolyhedralCone newLifted=PolyhedralCone(
					  expandVectorList(lifted.getHalfSpaces(),pivots),
					  equations,
					  dimension+pivots.size()+1
					  );
  return HalfOpenCone(dimension+pivots.size(),newLifted);
}

static void print(HalfOpenConeList const &l, Printer &P)
{
  P << "Printing half Open cone list-----------------------------------------------------------------\n";
  for(HalfOpenConeList::const_iterator i=l.begin();i!=l.end();i++)
    {
      i->lifted.print(&P);
    }
  P << "Done printing half Open cone list-----------------------------------------------------------------\n";
}
static void print2(vector<HalfOpenCone> const &l, Printer &P)
{
  P << "Printing half Open cone list---------------------------------------------------------------\n";
  for(vector<HalfOpenCone>::const_iterator i=l.begin();i!=l.end();i++)
    {
      i->lifted.print(&P);
    }
  P << "Done printing half Open cone list----------------------------------------------------------------!\n";
}

static IntegerVectorList liftEquations(IntegerVectorList const &l, FieldMatrix const &B)
{
  IntegerVectorList ret;
  list<int> nonPivots=B.nonPivotColumns();

  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      FieldVector I=integerVectorToFieldVector(i->subvector(0,i->size()-1),Q);
      IntegerVector v(B.getWidth()+1);
      FieldVector v2(Q,B.getWidth());
      for(int j=0;j<B.getWidth();j++)
	v2[j]=dot(B.canonicalize(integerVectorToFieldVector(IntegerVector::standardVector(B.getWidth(),j),Q)).subvector(nonPivots),I);
      IntegerVector v3=v2.primitive();
      for(int j=0;j<v3.size();j++)v[j]=v3[j];
      v[B.getWidth()]=(*i)[i->size()-1];
      ret.push_back(v);
    }

  return ret;
}

HalfOpenCone liftEquations(HalfOpenCone const &c, FieldMatrix const &B)
{
  PolyhedralCone newCone=PolyhedralCone(
					  liftEquations(c.lifted.getHalfSpaces(),B),
					  liftEquations(c.lifted.getEquations(),B),
					  B.getWidth()+1
					  );
  return HalfOpenCone(B.getWidth(),newCone);
}

/* Use full for transforming the cone back into original coordinates*/
class HalfOpenConeProcessorAffineChange : public HalfOpenConeProcessor
{
  IntegerVectorList interiorPoints;
  HalfOpenConeProcessor &parent;
  list<int> pivots;
  IntegerVectorList newEquations;
  FieldMatrix B;
public:
  HalfOpenConeProcessorAffineChange(HalfOpenConeProcessor &parent_, list<int> const &pivots_, IntegerVectorList const &newEquations_, FieldMatrix const &B_):
    parent(parent_),
    pivots(pivots_),
    newEquations(newEquations_),
    B(B_)
  {
    if(parent_.savePartialResult)setSave();
  }
  void process(HalfOpenCone const &c)
  {
    HalfOpenCone c2=liftEquations(c.rewriteExpand(pivots,newEquations),B);
    parent.process(c2);
  }
};


//HalfOpenConeList
void tropicalHyperSurfaceIntersectionInSubspace(int dimension, PolynomialSet const &G, HalfOpenCone *restrictingCone, HalfOpenConeProcessor &processor, int intervalLow, int intervalHigh)
{
  /*
    Removing the lineality space. (This could be carried out on the
    level of polyhedra, which may or may not be more efficient.)

    We first find generators for the linalityspace, write them as rows
    of a matrix B.  We reduce B and wish to ignore the coordinates
    whose columns contains a pivot.

    Afterwards we need to expand the inequalities/equations to the
    entire space. The matrix B tels us how to project to
    lowerdimensional space (by taking normals forms). The equations
    are gotten as the pull-back of the lower dimensional equations by
    this normal form map.
   */
  log2 cerr<<"Projecting Newton polytopes modulo the homogeneity space.";

  int N=G.getRing().getNumberOfVariables();
  IntegerVectorList w=wallInequalities(G);
  FieldMatrix W=integerMatrixToFieldMatrix(rowsToIntegerMatrix(w,N),Q);
  if(restrictingCone)
    {
      PolyhedralCone temp=restrictingCone->closure();
      W=combineOnTop(combineOnTop(W,
					  integerMatrixToFieldMatrix(rowsToIntegerMatrix(temp.getEquations(),N),Q)),
			                  integerMatrixToFieldMatrix(rowsToIntegerMatrix(temp.getHalfSpaces(),N),Q));
    }
  FieldMatrix B=W.reduceAndComputeKernel();
  B.reduce();
  list<int> BNonPivots=B.nonPivotColumns();
  PolynomialRing R(G.getRing().getField(),BNonPivots.size());

  PolynomialSet g=G.embeddedInto(R,&BNonPivots);

  HalfOpenCone restrictedConeNew(0,PolyhedralCone(1));
  if(restrictingCone)
    {
      restrictedConeNew=restrictingCone->withChosenCoordinates(BNonPivots);
      restrictingCone=&restrictedConeNew;
    }

  log2 cerr<<"Done projecting Newton polytopes modulo the homogeneity space.";
  /*
    Now do the computation with the new set of polynomials.
   */

  /*
    Here follows restriction to subspace cut out by the binomials.
   */
  log2 cerr<<"Restricting to subspace determined by binomials and computing tropical hypersurfaces.";

  int n=g.getRing().getNumberOfVariables();
  IntegerVectorList equations;

  if(restrictingCone)//add those from the restricting cone
    {
      equations=restrictingCone->shrink(restrictingCone->lifted.getEquations());
    }

  int numberOfNonBinomials=0;
  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    {
      if(i->numberOfTerms()==2)
	equations.push_back(binomialToIntegerVector(*i));
      else
	numberOfNonBinomials++;
    }
  FieldMatrix A=integerMatrixToFieldMatrix(rowsToIntegerMatrix(equations,n),Q);

  int npivots=A.reduceAndComputeRank();
  // Coordinates to remove are those corresponding to pivots
  // We want to rewrite inequalities without the non-pivots
  // Unfortunatetly we need to compute the fans first


  //  AsciiPrinter(Stderr)<<g;

  list<int> nonPivots=A.nonPivotColumns();

  //  log0 cerr<<"Coordinates left "<<nonPivots.size()<<endl;

  vector<vector<HalfOpenCone> > fullDimFanList;
  vector<vector<HalfOpenCone> > coDimOneFanList;
  vector<IntegerVector> parents(numberOfNonBinomials);

    list<int> pivots=A.pivotColumns();

  int I=0;
  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    if(i->numberOfTerms()!=2)
      {
	log2 cerr << I;
	HalfOpenConeList l,lf;
	log2 cerr<< "Building fan"<<endl;
	buildFanFromPolynomial(*i, &lf, &l, &(parents[I]),restrictingCone);
	log2 cerr<< "Number of cones:" << lf.size()<<","<< l.size()<<endl;
	log2 cerr<< "rewriting"<<endl;

	vector<HalfOpenCone> L;
	for(HalfOpenConeList::const_iterator i=l.begin();i!=l.end();i++)
	  L.push_back(i->rewrite(A,nonPivots));
	vector<HalfOpenCone> F;
	for(HalfOpenConeList::const_iterator i=lf.begin();i!=lf.end();i++)
	  F.push_back(i->rewrite(A,nonPivots));

	fullDimFanList.push_back(F);
	coDimOneFanList.push_back(L);
	log2 cerr<< "Done Building fan"<<endl;
	I++;
      }
  log2 cerr<<"Done computing tropical hypersurfaces.";

  /* Now we must create a new HalfOpenConeProcessor and tell it how to expand a cone.
     We must insert the pivot columns from A and add in the equtions gotten from the binomials.
   */
  A.removeZeroRows();
  IntegerVectorList newEquations=fieldMatrixToIntegerMatrixPrimitive(A).getRows();

  HalfOpenConeProcessorAffineChange myProcessor(processor,pivots,newEquations,B);

  RecursionData data(coDimOneFanList,myProcessor);
  data.extractInformationFromFullFan(fullDimFanList,parents);

//  data.completeTable();//HERE ENABLE HERE
//  data.table.print();//HERE


  IntegerVectorList empty;
  log2 cerr<<"Doing intersection.";
  if(intervalHigh!=-1)data.setInterval(intervalLow,intervalHigh);
  data.rek(0, HalfOpenCone(nonPivots.size(),empty,empty,empty));
  log2 cerr<<"Done doing intersection.";
  //  data.rek(0, HalfOpenCone(dimension,empty,empty,empty));
  //  log2 fprintf(Stderr,"LPs solved:%i for relation table\n",data.table.numberOfSolvedLPs);

  /*  HalfOpenConeList ret;

  for(HalfOpenConeList::const_iterator i=data.ret.begin();i!=data.ret.end();i++)
    {
      //      ret.push_back(i->rewriteExpand(pivots,newEquations));
      ret.push_back(liftEquations(i->rewriteExpand(pivots,newEquations),B));
    }

  return ret;
  */
}


PolyhedralFan faceComplexOfCone(HalfOpenCone &c)
{
  PolyhedralFan ret(c.dimension);
  list<HalfOpenCone> h;
  c.splitIntoRelativelyOpenCones(h);
  for(list<HalfOpenCone>::iterator i=h.begin();i!=h.end();i++)
    {
      PolyhedralCone temp=i->closure();
      temp.canonicalize();
      ret.insert(temp);
    }
  return ret;
}


bool operator<(HalfOpenCone const &a, HalfOpenCone const &b)
{
  return a.lifted<b.lifted;
}


class SingleSubspaceProcessor : public HalfOpenConeProcessor
{
  int desiredDimension;
  bool desiredDimensionFound;
public:
  SingleSubspaceProcessor(int desiredDimension_):
    desiredDimension(desiredDimension_),
    desiredDimensionFound(false)
  {
  }
  void process(HalfOpenCone const &c)
  {
    if(c.lifted.dimension()-1==desiredDimension)
      {
        desiredDimensionFound=true;
//        return false;
      }
//    return true;
  }
  bool success()
  {
    return desiredDimensionFound;
  }
};

bool nonEmptyStableIntersection(PolynomialSet const &g)
{
  vector<vector<HalfOpenCone> > coDimOneFanList;

  {
    int I=0;
    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
      {
        HalfOpenConeList l;
        PolynomialSet A(g.getRing());
        A.push_back(*i);
        A.markAndScale(LexicographicTermOrder());
        PolyhedralCone h=homogeneitySpace(A);
        h.canonicalize();
        IntegerVectorList equations=h.getEquations();
        for(IntegerVectorList::const_iterator j=equations.begin();j!=equations.end();j++)
          {
            IntegerVectorList eq2;
            eq2.push_back(*j);
            HalfOpenCone C(j->size(),eq2,IntegerVectorList(),IntegerVectorList());
            l.push_back(C);
          }


        vector<HalfOpenCone> L;
        for(HalfOpenConeList::const_iterator i=l.begin();i!=l.end();i++)
          {
            L.push_back(*i);
          }
        coDimOneFanList.push_back(L);
        I++;
      }
  }

  SingleSubspaceProcessor p(g.getRing().getNumberOfVariables()-g.size());
  RecursionData data(coDimOneFanList,p);

  IntegerVectorList empty;
  data.rek(0, HalfOpenCone(g.getRing().getNumberOfVariables(),empty,empty,empty));

  log2 fprintf(Stderr,"LPs solved:%i for relation table\n",data.table.numberOfSolvedLPs);

  return p.success();
}
