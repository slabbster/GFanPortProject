#include "halfopencone.h"

#include "buchberger.h"
#include "enumeration.h"
#include "reversesearch.h"
#include "wallideal.h"

#include "printer.h"
#include "parser.h"

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


bool HalfOpenCone::contains(IntegerVector const &v)const
{
  IntegerVectorList inequalityList=lifted.getHalfSpaces();
  IntegerVectorList equationList=lifted.getLinealitySpace();

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


HalfOpenCone::HalfOpenCone(int dimension_, IntegerVectorList const &equations, IntegerVectorList const &nonstrict, IntegerVectorList const &strict, bool findFacets):
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
    if(!(t&7))ret.lifted.findFacets();

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
      ret.push_back(HalfOpenCone(dimension,equationList,nonStrictList,strictList,true));
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


HalfOpenConeList tropicalHyperSurfaceIntersection2(int dimension, PolynomialSet const &g)
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
  //  fprintf(Stderr,"BEGIN\n");
  //  AsciiPrinter P(Stderr);
  //  print(P);
  lifted.findFacets();
  //  print(P);
  
  /*	{
	  IntegerVector v=StringParser("(3,0,3,2,0,3)").parseIntegerVector();
	  if(contains(v))fprintf(stderr,"??????????????????????????????????????????\n");
	  }*/

  IntegerVectorList inequalityList=lifted.getHalfSpaces();
  IntegerVectorList equationList=lifted.getLinealitySpace();

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
      assert(i->dimensionOfLargestContainedSubspace()==homog);
    }
  fprintf(stderr,"Ambient dimension: %i, maximal dimension: %i, dimension of lineality space: %i\n",ambientDimension,largest,homog);

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
	//	fprintf(stderr,"%i\n",i);
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
    for(int c2=0;c2<fanList[fan2].size();c2++)
      intersectTriviallyInIntersection(fan1,cone1,fan2,c2);    

    return knownEmptyIntersectionInIntersection.nonCandidates(fan1,cone1,fan2);
  }
  void markNoIntersectionInIntersection(int fan1, int cone1, int fan2, int cone2)
  {
    knownEmptyIntersectionInIntersection.set(fan1,cone1,fan2,cone2);
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
public:
  RelationTable table;
  RecursionData(vector<vector<HalfOpenCone> > const &fans_):
    table(fans_),
    fans(fans_),
    chosen(fans_.size()),
    chosenFans(fans_.size()),
    usedFans(fans_.size()),
    iterators(fans_.size()),
    nCandidates(fans_.size()),
    numberOfUsefulCalls(0),
    totalNumberOfCalls(0)
  {
  }

  HalfOpenConeList ret;

  BitSet computeCandidates(int index, int fanNumber)
  {
    BitSet nonCandidates(fans[fanNumber].size());
    for(int i=0;i<index;i++)
      {
	nonCandidates.add(table.getNonCandidates(chosenFans[i],chosen[i],fanNumber));
      }
    return nonCandidates.negated();
  }
  bool rek(int index, HalfOpenCone const &current)
  {
    totalNumberOfCalls++;

    bool success=false;

    if(index == fans.size())
      {
	fprintf(Stderr,"ADDING CONE\n");
	AsciiPrinter P(Stderr);
	//  current.print(P);
	ret.push_back(current);
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


	nCandidates[index]=bestNumberOfCandidates;//just for printing

	static int iterationNumber;
	if(!(iterationNumber++ & 31))
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

		  HalfOpenCone next=intersection(current,fans[chosenFans[index]][i],false);
		  if((fans.size()!=10)&&(fans.size()!=9))
		    {
		      bool s=rek(index+1,next);
		      success|=s;
		    }
		  else
		    {
		      vector<vector<HalfOpenCone> > L2;
		      {
			for(int i=0;i<fans.size();i++)
			  if(!usedFans[i])
			    {
			      vector<HalfOpenCone> L;
			      for(vector<HalfOpenCone>::const_iterator j=fans[i].begin();j!=fans[i].end();j++)
				{
				  if(!haveEmptyIntersection(next,*j))
				    {
				      L.push_back(intersection(next,*j,true));
				    }
				}
			      fprintf(stderr,"New fan size:%i\n",L.size()); 
			      L2.push_back(L);
			    }
		      }
		      RecursionData data(L2);
		      data.completeTable();
		      success|=data.rek(0, HalfOpenCone(next.dimension,IntegerVectorList(),IntegerVectorList(),IntegerVectorList()));
		      ret.splice(ret.begin(),data.ret);
		    }
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
		      fprintf(stderr,"ONE FOR FREE\n");
		    }
		}
		}*/

	nCandidates[index]=-1;//just for printing
	iterators[index]=0;//just for printing

	usedFans[chosenFans[index]]=false;
	chosenFans[index]=-1;
      }
    if(success)numberOfUsefulCalls++;
    return success;
  }
  void transitiveClosure()
  {
    bool rep;
    do{
      rep=false;
      int a=0;
    for(int f1=0;f1<fans.size();f1++)
      {
	//	fprintf(stderr,"%i\n",f1);
	for(int f2=f1+1;f2<fans.size();f2++)
	  for(int c1=0;c1<fans[f1].size();c1++)
	    for(int c2=0;c2<fans[f2].size();c2++)
	      {
		if(!table.intersectTriviallyInIntersection(f1,c1,f2,c2))
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
		      }
		    if(dontintersect)
		      {
			table.markNoIntersectionInIntersection(f1,c1,f2,c2);
			rep=true;
		      }
		  }
	      }
      }
    fprintf(stderr,"%i FOR FREE\n",a);
    }
    while(rep);
  }

  void completeTable()
  {
    for(int f1=0;f1<fans.size();f1++)
      {
	fprintf(stderr,"%i\n",f1);
	for(int f2=f1+1;f2<fans.size();f2++)
	for(int c1=0;c1<fans[f1].size();c1++)
	  for(int c2=0;c2<fans[f2].size();c2++)
	    table.intersectTriviallyInIntersection(f1,c1,f2,c2);    
      }
    transitiveClosure();
  }
};

HalfOpenConeList tropicalHyperSurfaceIntersection(int dimension, PolynomialSet const &g)
{
  vector<vector<HalfOpenCone> > L2;
  {
    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
      {
	HalfOpenConeList l=tropicalHyperSurface(*i);
	vector<HalfOpenCone> L;
	for(HalfOpenConeList::const_iterator i=l.begin();i!=l.end();i++)
	  {
	    L.push_back(*i);
	  }
	L2.push_back(L);
      }
  }

  RecursionData data(L2);
  //  data.completeTable();
  data.rek(0, HalfOpenCone(dimension,IntegerVectorList(),IntegerVectorList(),IntegerVectorList()));
  fprintf(stderr,"LPs solved:%i for relation table\n",data.table.numberOfSolvedLPs);
  return data.ret;
}
