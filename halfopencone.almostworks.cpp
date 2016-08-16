#include "halfopencone.h"

#include "buchberger.h"
#include "enumeration.h"
#include "reversesearch.h"
#include "wallideal.h"

#include "printer.h"

#include "lp.h"


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


HalfOpenCone::HalfOpenCone(int dimension_, IntegerVectorList const &equations, IntegerVectorList const &nonstrict, IntegerVectorList const &strict):
  dimension(dimension_),
  liftedDimension(dimension_+1),
  lifted(dimension_+1)
{
  IntegerVectorList equationList,inequalityList;

  appendList(equationList,equations,0);
  appendList(inequalityList,nonstrict,0);
  appendList(inequalityList,strict,-1);

  //      AsciiPrinter(Stderr).printVectorList(inequalityList);
  //      AsciiPrinter(Stderr).printVectorList(equationList);
  //      AsciiPrinter(Stderr).printInteger(liftedDimension);
  lifted=PolyhedralCone(inequalityList,equationList,liftedDimension);
  //  lifted.findFacets();
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
				   swapFirstLast(lifted.getLinealitySpace())
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
  AsciiPrinter(Stderr).printVectorList(lifted.getLinealitySpace());
  fprintf(Stderr,"hasSolution=%i\n",ret1);
  */
  //  assert(ret1==ret2);

  return ret1;
}


bool haveEmptyIntersection(const HalfOpenCone &a, const HalfOpenCone &b)
{
  assert(a.dimension==b.dimension);
  IntegerVectorList inequalityList=a.lifted.getHalfSpaces();
  IntegerVectorList equationList=a.lifted.getLinealitySpace();

  IntegerVectorList inequalityList2=b.lifted.getHalfSpaces();
  IntegerVectorList equationList2=b.lifted.getLinealitySpace();

  inequalityList.splice(inequalityList.begin(),inequalityList2);
  equationList.splice(equationList.begin(),equationList2);

  bool ret1=!hasHomogeneousSolution(a.liftedDimension,swapFirstLast(inequalityList),swapFirstLast(equationList));

  /*  HalfOpenCone c=intersection(a,b);
  if(c.isEmpty()!=ret1)
    {
      AsciiPrinter(Stderr).printVectorList(inequalityList);
      AsciiPrinter(Stderr).printVectorList(equationList);

      AsciiPrinter(Stderr).printVectorList(c.lifted.getHalfSpaces());
      AsciiPrinter(Stderr).printVectorList(c.lifted.getLinealitySpace());

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

HalfOpenCone intersection(const HalfOpenCone &a, const HalfOpenCone &b)
{
  assert(a.dimension==b.dimension);
  /*  
  fprintf(Stderr,"-----------------------------------------------------------\n");
  fprintf(Stderr,"Intersecting:\n");
  AsciiPrinter P(Stderr);
  printHalfOpenCone(P,a);
  printHalfOpenCone(P,b);
  */
  HalfOpenCone ret=HalfOpenCone(a.dimension,intersection(a.lifted,b.lifted));

  {
    static int t;
    t++;
    if(!(t&7))ret.lifted.findFacets();

    //1 4:53
    //3 3:38
    //7
  }

  /*
  fprintf(Stderr,"Result:\n");
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
  return PolyhedralCone(shrink(lifted.getHalfSpaces()),shrink(lifted.getLinealitySpace()),dimension);
}


/*
 */
HalfOpenConeList orientedBoundary(PolyhedralCone C, TermOrder const &t)
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

  // Let's make the non-strict inequalities strict one at a time and add a cone for each iteration
  while(!nonStrictList.empty())
    {
      IntegerVector v=nonStrictList.front();
      nonStrictList.pop_front();


      IntegerVectorList equationList;
      equationList.push_back(v);
      ret.push_back(HalfOpenCone(dimension,equationList,nonStrictList,strictList));
      strictList.push_back(v);
    }
  return ret;
}


HalfOpenConeList tropicalHyperSurface(Polynomial const &p1)
{
  Polynomial p=p1.homogenization();
  HalfOpenConeList ret;
  PolynomialSet g;
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


HalfOpenConeList tropicalHyperSurfaceIntersection(int dimension, PolynomialSet const &g)
{
  HalfOpenConeList intersection;
  intersection.push_back(HalfOpenCone(dimension,IntegerVectorList(),IntegerVectorList(),IntegerVectorList()));

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


PolyhedralFan tropicalHyperSurfaceIntersectionClosed(int dimension, PolynomialSet const &g)
{
  HalfOpenConeList intersection=tropicalHyperSurfaceIntersection(dimension,g);


  AsciiPrinter P(Stderr);
  printHalfOpenConeList(intersection, P);

  PolyhedralFan ret(dimension);
  for(HalfOpenConeList::iterator i=intersection.begin();i!=intersection.end();i++)
    {
      PolyhedralCone c=i->closure();
      c.canonicalize();
      ret.insert(c);
    }

  return ret;
}


void HalfOpenCone::splitIntoRelativelyOpenCones(list<HalfOpenCone> &l)
{
  fprintf(Stderr,"BEGIN\n");
  AsciiPrinter P(Stderr);
  print(P);
  lifted.findFacets();
  print(P);
  
  IntegerVectorList inequalityList=lifted.getHalfSpaces();
  IntegerVectorList equationList=lifted.getLinealitySpace();

  IntegerVectorList strict,nonstrict;
  

  for(IntegerVectorList::const_iterator i=inequalityList.begin();i!=inequalityList.end();i++)
    if((*i)[i->size()-1]<0)
      strict.push_back(*i);
    else
      nonstrict.push_back(*i);

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
  fprintf(Stderr,"END\n");
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
      fprintf(Stderr,"---------------------------------------------------------------\n");
      temp.print(P);
      fprintf(Stderr,"---------------------------------------------------------------\n");
      temp.splitIntoRelativelyOpenCones(tempSplit);

      fprintf(Stderr,"Splits into:");
      for(HalfOpenConeList::const_iterator i=tempSplit.begin();i!=tempSplit.end();i++)
	i->print(P);
      fprintf(Stderr,"Splits into End.");

      ret.splice(ret.begin(),tempSplit);

      fprintf(Stderr,"B\n");
    }

  return ret;
}


void printHalfOpenConeList(HalfOpenConeList const &l, class Printer & p)
{
  HalfOpenConeList L=splitIntoRelativelyOpenCones(l);

  list<PolyhedralCone> cones;
  for(HalfOpenConeList::iterator i=L.begin();i!=L.end();i++)
    cones.push_back(i->closure());

  int homog=1000000;
  int largest=0;
  for(list<PolyhedralCone>::const_iterator i=cones.begin();i!=cones.end();i++)
    {
      if(i->dimension()<homog)homog=i->dimension();
      if(i->dimension()>largest)largest=i->dimension();
    }
  assert(homog!=1000000);

  fprintf(stderr,"homog:%i largest:%i\n",homog,largest);

  IntegerVectorList rays;
  for(list<PolyhedralCone>::iterator i=cones.begin();i!=cones.end();i++)
    if(i->dimension()==homog+1)rays.push_back(i->getRelativeInteriorPoint());
  
  p.printString("Rays:\n");

  p.printVectorList(rays);

  for(int d=homog;d<=largest;d++)
    {
      list<PolyhedralCone> cones2;
      for(list<PolyhedralCone>::iterator i=cones.begin();i!=cones.end();i++)
	if(i->dimension()==d)
	  cones2.push_back(*i);

      p.printString("Printing ");p.printInteger(d);p.printString("dimensional cones.\n");
      p.printString("Number of cones:");p.printInteger(d);p.printString("\n");
      
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
	  p.printVector(v);
	}
    }
}
