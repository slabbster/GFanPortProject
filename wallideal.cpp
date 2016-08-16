#include "wallideal.h"

#include "division.h"
#include "buchberger.h"
#include "timer.h"
#include "printer.h"
#include "lp.h"
#include "polyhedralcone.h"
#include "tropical2.h"
#include "linalg.h"
#include "log.h"

#include <set>
#include <iostream>

static Timer flipTimer("Flip",10);
static Timer flipTimer1("Flip1",10);
static Timer flipTimer2("Flip2",10);
static Timer flipTimer3("Flip3",10);
static Timer flipTimer4("Flip4",10);
static Timer flipTimer5("Flip5",10);
static Timer coneTimer("fajskfda",10);

Polynomial wallPolynomial(Polynomial const &p, IntegerVector const &wallNormal)
{
  Polynomial r(p.getRing());
  IntegerVector markedExponent=p.getMarked().m.exponent;

  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      IntegerVector dif=markedExponent-i->first.exponent;
      if(dependent(dif,wallNormal))
        r+=Polynomial(Term(i->second,i->first));
    }

  r.mark(Monomial(r.getRing(),markedExponent));

  return r;
}

static Polynomial wallPolynomial(Polynomial const &p, IntegerVectorList const &wallEqualities)
{
  Polynomial r(p.getRing());
  IntegerVector markedExponent=p.getMarked().m.exponent;

  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      IntegerVector dif=markedExponent-i->first.exponent;
      bool dep=false;
      for(IntegerVectorList::const_iterator j=wallEqualities.begin();j!=wallEqualities.end();j++)
	{
	  if(dependent(dif,*j))
	    {
	      dep=true;
	      break;
	    }
	}

      if(dep || dif.isZero())
	r+=Polynomial(Term(i->second,i->first));
    }

  r.mark(Monomial(p.getRing(),markedExponent));

  return r;
}

PolynomialSet wallIdeal(PolynomialSet const &groebnerBasis, IntegerVector const &wallNormal)
{
  PolynomialRing theRing=groebnerBasis.getRing();
  PolynomialSet r(theRing);

  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    {
      r.push_back(wallPolynomial(*i,wallNormal));
    }
  return r;
}

PolynomialSet lowerDimensionalWallIdeal(PolynomialSet const &groebnerBasis, IntegerVectorList const &wallEqualities)
{
  PolynomialRing theRing=groebnerBasis.getRing();
  PolynomialSet r(theRing);

  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    {
      r.push_back(wallPolynomial(*i,wallEqualities));
    }
  return r;
}

IntegerVectorList wallRemoveScaledInequalities(IntegerVectorList const &l)
{
  IntegerVectorList ret;

  /*for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      bool found=false;

      assert(!i->isZero());

      for(IntegerVectorList::const_iterator k=ret.begin();k!=ret.end();k++)
	if(dependent(*i,*k)&&dotLong(*i,*k)>0)found=true;

      if(!found)ret.push_back(*i);
    }
      */
  set<IntegerVector> temp;
  //	    log0 fprintf(Stderr,"C\n");
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)temp.insert(normalized(*i));
  //	    log0 fprintf(Stderr,"D\n");
  for(set<IntegerVector>::const_iterator i=temp.begin();i!=temp.end();i++)ret.push_back(*i);

  return ret;
}


IntegerVectorList algebraicTest(IntegerVectorList const &l, PolynomialSet const &groebnerBasis) // TO DO: FIGURE OUT IF THIS TEST WORKS IN THE NON-HOMOGENEOUS CASE
{
  PolynomialRing theRing=groebnerBasis.getRing();
  LexicographicTermOrder T;
  PolynomialSet LT(theRing);

  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    {
      LT.push_back(Polynomial(i->getMarked()));
    }


  //fprintf(Stderr,"In Size:%i\n",l.size());
  IntegerVectorList ret;
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      bool accept=true;
      PolynomialSet g2=wallIdeal(groebnerBasis,*i);
      for(PolynomialSet::const_iterator j=g2.begin();j!=g2.end() && accept;j++)
	for(PolynomialSet::const_iterator k=j;k!=g2.end();k++)
	  {
	    //	    fprintf(Stderr,"test\n");
	    if((!j->isMonomial()) || (!k->isMonomial()))
	    if(!relativelyPrime(j->getMarked().m.exponent,k->getMarked().m.exponent))
	      {
		Polynomial s=sPolynomial(*j, *k);
		if(!s.isZero())
		  {
		    bool witness=true;
		    for(TermMap::const_iterator i=s.terms.begin();i!=s.terms.end();i++)//TODO: prove that it suffices to check just the leading term of the S-polynomial and adjust the code accordingly.
		      {
			bool inideal=false;
			for(PolynomialSet::const_iterator j=LT.begin();j!=LT.end();j++)
			  if(j->getMarked().m.exponent.divides(i->first.exponent))
			    {
			      inideal=true;
			      break;
			    }
			if(inideal)
			  {
			    witness=false;
			    break;
			  }
		      }
		    if(witness)
		      {
			accept=false;
			break;
		      }

/*
		    s.mark(T); // with respect to some termorder
		    s.scaleMarkedCoefficientToOne();

		    if(1)
		      {
			Polynomial t=division(s,LT,T);
			if((t-s).isZero())
			  {
			    accept=false;
			    break;
			  }
		      }
		    	    else
		      {
			s=division(s,g2,T);
			if(!s.isZero())
			  {
			    accept=false;
			    break;
			  }
			  }*/
		  }
	      }

	  }
      if(accept)ret.push_back(*i);
    }
  //  fprintf(Stderr,"Out Size:%i\n",ret.size());
  return ret;
}


IntegerVectorList exponentDifferences(PolynomialSet const &groebnerBasis)
{
	IntegerVectorList ret;
	set<IntegerVector> temp;
	for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
	    {
		IntegerVector markedExponent=i->getMarked().m.exponent;

		for(TermMap::const_iterator j=i->terms.begin();j!=i->terms.end();j++)
		{
			IntegerVector dif=markedExponent-j->first.exponent;

			if(!dif.isZero())
			{
				/* bool found=false; //These four lines are now done in wallRemoveScaledInequalities()

	              for(IntegerVectorList::const_iterator k=ret.begin();k!=ret.end();k++)
	                if(dependent(dif,*k))found=true;

	              if(!found)
		      */
				temp.insert(dif);
			}
		}
	    }
	for(set<IntegerVector>::const_iterator i=temp.begin();i!=temp.end();i++)ret.push_back(*i);
	return ret;
}

IntegerVectorList wallInequalities(PolynomialSet const &groebnerBasis)
{
  return wallRemoveScaledInequalities(exponentDifferences(groebnerBasis));
  //  return algebraicTest(wallRemoveScaledInequalities(ret),groebnerBasis);
}


PolynomialSet flip(PolynomialSet const &groebnerBasis, IntegerVector const &wallNormal, TermOrder *autoReduceHint)
{
  PolynomialRing theRing=groebnerBasis.getRing();
  //  fprintf(Stderr,"ENTERING flip\n");
  //  fprintf(Stderr,"flip - start\n");
  //  AsciiPrinter(Stderr).printPolynomialSet(groebnerBasis);
  //  AsciiPrinter(Stderr).printVector(wallNormal);

  TimerScope ts(&flipTimer);
  //Subroutine 3.7 in [Sturmfels]

  // Step 1
  //  fprintf(Stderr,"flip - step1\n");
  flipTimer1.on();
  PolynomialSet wall=wallIdeal(groebnerBasis,wallNormal);
  wall.markAndScale(WeightTermOrder(wallNormal));// This marking will be used later on when we lift
  //  fprintf(Stderr,"Changed order:\n");
  // AsciiPrinter(Stderr).printPolynomialSet(wall);
  flipTimer1.off();

  // Step 2
  //  fprintf(Stderr,"flip - step2\n");
  // AsciiPrinter(Stderr).printPolynomialSet(wall);
  flipTimer2.on();
  PolynomialSet oldWall=wall;
  WeightTermOrder flipOrder(-wallNormal);
  buchberger(&wall,flipOrder);
  flipTimer2.off();

  // Step 3
  //  fprintf(Stderr,"flip - step3\n");
  flipTimer3.on();
  PolynomialSet newBasis(theRing);
  flipTimer3.off();

  // Step 4 Lift
  //  fprintf(Stderr,"flip - lifting\n");
  //  fprintf(Stderr,"flip - step4\n");

  {
    // liftBasis() could be used for this!!!!

    TimerScope ts(&flipTimer4);
    for(PolynomialSet::const_iterator j=wall.begin();j!=wall.end();j++)
      {
	newBasis.push_back(divisionLift(*j, oldWall, groebnerBasis, LexicographicTermOrder()));
	/*{
	    // The following should also work:
	    Polynomial q=*j-division(*j,groebnerBasis,LexicographicTermOrder());
	    assert(!q.isZero());

	    newBasis.push_back(q);
	    }*/
      }
  }

  // Step 5 Autoreduce
  //  fprintf(Stderr,"flip - autoreduce\n");
  //  AsciiPrinter(Stderr).printPolynomialSet(wall);
  //  AsciiPrinter(Stderr).printPolynomialSet(newBasis);
  //  fprintf(Stderr,"flip - step5\n");
  {
    TimerScope ts(&flipTimer5);
    PolynomialSet::const_iterator k=wall.begin();
    for(PolynomialSet::iterator j=newBasis.begin();j!=newBasis.end();j++)
      {
	j->mark(k->getMarked().m);
	k++;
      }

    //  fprintf(Stderr,"Marked order:\n");
    //  AsciiPrinter(Stderr).printPolynomialSet(wall);
    //  AsciiPrinter(Stderr).printPolynomialSet(newBasis);
    //  fprintf(Stderr,"Not reduced lifted basis:\n");
    //  AsciiPrinter(Stderr).printPolynomialSet(newBasis);

    static int t;
    t++;
    t&=0;
    if(t==0)
      {
	// fprintf(Stderr,"autoreducing ..");
    	if(autoReduceHint)
    		autoReduce(&newBasis,*autoReduceHint);
    	else
    		autoReduce(&newBasis,LexicographicTermOrder());
	// fprintf(Stderr,".. done\n");
      }

    //autoReduce(&newBasis,StandardGradedLexicographicTermOrder());
  }

  //  fprintf(Stderr,"flip - done\n");
  //  fprintf(Stderr,"LEAVING flip\n");

  //fprintf(Stderr,"%i",newBasis.size());

  return newBasis;
}


bool wallContainsPositiveVector(IntegerVector const &wallNormal)
{
  //This is not right I think
  int n=wallNormal.size();
  for(int i=0;i<n;i++)if(wallNormal[i]<0)return true;

  return false;
}

void wallAddCoordinateWalls(IntegerVectorList &normals)
{
  assert(!normals.empty());
  int n=normals.begin()->size();
  for(int i=0;i<n;i++)normals.push_back(IntegerVector::standardVector(n,i));
}


bool isIdealHomogeneous(PolynomialSet const &groebnerBasis)
{
  int n=groebnerBasis.numberOfVariablesInRing();
  IntegerVectorList a;
  PolyhedralCone p=intersection(PolyhedralCone(a,wallInequalities(groebnerBasis),n),PolyhedralCone::positiveOrthant(n));

  return p.containsPositiveVector();
}


/* This routine is a preprocessing step for redudancy removal.
   The routine normalizes the list of vectors in gcd sense.
   It removes duplicates.
   It removes direction that are sums of other directions.
   The input must satisfy:
     Input must be pointed, meaning that there must exist a codimension one subspace with all the input vectors strictly on one side.
     It particular, there is no zero-vector in the list (It would be easy to change the routine to handle this case)
   These requirements guarantee that *i and *j are not removed in the line with the comment.
 */
/*
  There are two versions of this routine.
 */
//400 -> 80
IntegerVectorList normalizedWithSumsAndDuplicatesRemoved2(IntegerVectorList const &a)
{
  IntegerVectorList ret;
  set<IntegerVector> b;

  for(IntegerVectorList::const_iterator i=a.begin();i!=a.end();i++)
    {
      assert(!(i->isZero()));
      b.insert(normalized(*i));
    }

  for(set<IntegerVector>::const_iterator i=b.begin();i!=b.end();i++)
    for(set<IntegerVector>::const_iterator j=i;j!=b.end();j++)
	if(i!=j)b.erase(normalized(*i+*j));//this can never remove *i or *j

  for(set<IntegerVector>::const_iterator i=b.begin();i!=b.end();i++)
    ret.push_back(*i);

  return ret;
}


class MyHashMap
{
	typedef vector<set<IntegerVector> > Container;
	Container container;
	int tableSize;
public:
	class iterator
	{
		class MyHashMap &hashMap;
		int index; // having index==-1 means that we are before/after the elements.
		set<IntegerVector>::iterator i;
	public:
		bool operator++()
		{
			if(index==-1)goto search;
			i++;
			while(i==hashMap.container[index].end())
			{
				search:
				index++;
				if(index>=hashMap.tableSize){
					index=-1;
					return false;
				}
				i=hashMap.container[index].begin();
			}
			return true;
		}
		IntegerVector const & operator*()const
		{
			return *i;
		}
		IntegerVector operator*()
		{
			return *i;
		}
/*		bool operator()()
		{
			return index!=-1;
		}*/
		iterator(MyHashMap &hashMap_):
			hashMap(hashMap_)
			{
				index=-1;
			}
	};
	unsigned int function(const IntegerVector &v)
	{
		unsigned int ret=0;
		int n=v.size();
		for(int i=0;i<n;i++)
			ret=(ret<<3)+(ret>>29)+v.UNCHECKEDACCESS(i);
		return ret%tableSize;
	}
	MyHashMap(int tableSize_):
		container(tableSize_),
		tableSize(tableSize_)
		{
			assert(tableSize_>0);
		}
	void insert(const IntegerVector &v)
	{
		container[function(v)].insert(v);
	}
	void erase(IntegerVector const &v)
	{
		container[function(v)].erase(v);
	}
	iterator begin()
	{
		iterator ret(*this);
		++ret;
		return ret;
	}
	int size()
	{
		iterator i=begin();
		int ret=0;
		do{ret++;}while(++i);
		return ret;
	}
	void print()
	{
		for(int i=0;i<container.size();i++)
		{
			debug.printInteger(i);
			debug<<":\n";
			for(set<IntegerVector>::const_iterator j=container[i].begin();j!=container[i].end();j++)debug<<*j;
		}
	}
};


IntegerVectorList normalizedWithSumsAndDuplicatesRemoved(IntegerVectorList const &a)
{
	if(a.empty())return a;
	int n=a.front().size();
	IntegerVector temp1(n);
	IntegerVector temp2(n);
	IntegerVectorList ret;
	MyHashMap b(a.size());

//	  log0 debug<<"a.size"<<(int)a.size()<<"\n";
//	  log0 debug<<"b.size"<<(int)b.size()<<"\n";

	  //	log0 cerr << "N:"<<a.size();
	for(IntegerVectorList::const_iterator i=a.begin();i!=a.end();i++)
	{
		assert(!(i->isZero()));
		b.insert(normalized(*i));
    }
//	b.print();
//debug<<"TEST";
  //  log0 cerr << "N:"<<b.size();

    {
    	MyHashMap::iterator i=b.begin();

    	do
    	{
    		MyHashMap::iterator j=i;
    		while(++j)
    		{
    			//    			b.erase(normalized(*i+*j));//this can never remove *i or *j
    			IntegerVector const &I=*i;
    			IntegerVector const &J=*j;
    			for(int k=0;k<n;k++)temp1[k]=I.UNCHECKEDACCESS(k)+J.UNCHECKEDACCESS(k);
    			normalizedLowLevel(temp1,temp2);
    			b.erase(temp2);//this can never remove *i or *j
    		}
    		}
    	while(++i);
//    	log0 cerr << "N:"<<b.size();
    }
  vector<IntegerVector> original;
  {
  	MyHashMap::iterator i=b.begin();

  	do
  		{
    original.push_back(*i);
  		}
  	while(++i);
  		}

  for(vector<IntegerVector>::const_iterator i=original.begin();i!=original.end();i++)
    for(list<IntegerVector>::const_iterator j=a.begin();j!=a.end();j++)
      /*  for(set<IntegerVector>::const_iterator i=b.begin();i!=b.end();i++)
	  for(set<IntegerVector>::const_iterator j=i;j!=b.end();j++)*/
      //	if(*i!=*j)b.erase(normalized(*i+*j));//this can never remove *i or *j
      if(!dependent(*i,*j))
    	  {
//    	  b.erase(normalized(*i+*j));//this can never remove *i or *j
			IntegerVector const &I=*i;
			IntegerVector const &J=*j;
			for(int k=0;k<n;k++)temp1[k]=I.UNCHECKEDACCESS(k)+J.UNCHECKEDACCESS(k);
			normalizedLowLevel(temp1,temp2);
    	  b.erase(temp2);//this can never remove *i or *j
    	  }

//    log0 cerr << "N:"<<b.size();

//  log0 debug<<"b.size"<<(int)b.size()<<"\n";
    {
    	MyHashMap::iterator i=b.begin();

    	do
    	{
        	IntegerVector temp=*i;
        	ret.push_back(temp);
        }
        while(++i);
    }
  return ret;
}




// 400 -> 20
IntegerVectorList normalizedWithSumsAndDuplicatesRemovedNOHASH(IntegerVectorList const &a)
{
  /*  IntegerVectorList a;
  {
    set<IntegerVector> b;
    for(IntegerVectorList::const_iterator i=A.begin();i!=A.end();i++)b.insert(*i);
    for(set<IntegerVector>::const_iterator i=b.begin();i!=b.end();i++)a.push_back(*i);
    }*/

  IntegerVectorList ret;
  set<IntegerVector> b;

  log0 cerr << "N:"<<a.size();
  for(IntegerVectorList::const_iterator i=a.begin();i!=a.end();i++)
    {
      assert(!(i->isZero()));
      b.insert(normalized(*i));
    }

    log0 cerr << "N:"<<b.size();

  for(set<IntegerVector>::const_iterator i=b.begin();i!=b.end();i++)
    for(set<IntegerVector>::const_iterator j=i;j!=b.end();j++)
	if(i!=j)b.erase(normalized(*i+*j));//this can never remove *i or *j

  log0 cerr << "N:"<<b.size();

  vector<IntegerVector> original;
  for(set<IntegerVector>::const_iterator i=b.begin();i!=b.end();i++)
    original.push_back(*i);

  for(vector<IntegerVector>::const_iterator i=original.begin();i!=original.end();i++)
    for(list<IntegerVector>::const_iterator j=a.begin();j!=a.end();j++)
      /*  for(set<IntegerVector>::const_iterator i=b.begin();i!=b.end();i++)
	  for(set<IntegerVector>::const_iterator j=i;j!=b.end();j++)*/
      //	if(*i!=*j)b.erase(normalized(*i+*j));//this can never remove *i or *j
      if(!dependent(*i,*j))b.erase(normalized(*i+*j));//this can never remove *i or *j

    log0 cerr << "N:"<<b.size();

  for(set<IntegerVector>::const_iterator i=b.begin();i!=b.end();i++)
    ret.push_back(*i);

  return ret;
}


//400 -> 16
IntegerVectorList normalizedWithSumsAndDuplicatesRemoved3(IntegerVectorList const &a)
{
  IntegerVectorList ret;
  set<IntegerVector> b;

  for(IntegerVectorList::const_iterator i=a.begin();i!=a.end();i++)
    {
      assert(!(i->isZero()));
      b.insert(normalized(*i));
    }

  vector<IntegerVector> original;
  for(set<IntegerVector>::const_iterator i=b.begin();i!=b.end();i++)
    original.push_back(*i);

  for(vector<IntegerVector>::const_iterator i=original.begin();i!=original.end();i++)
    for(vector<IntegerVector>::const_iterator j=i;j!=original.end();j++)
      /*  for(set<IntegerVector>::const_iterator i=b.begin();i!=b.end();i++)
	  for(set<IntegerVector>::const_iterator j=i;j!=b.end();j++)*/
	if(i!=j)b.erase(normalized(*i+*j));//this can never remove *i or *j

  for(set<IntegerVector>::const_iterator i=b.begin();i!=b.end();i++)
    ret.push_back(*i);

  return ret;
}


IntegerVectorList wallFlipableNormals(PolynomialSet const &groebnerBasis, bool isKnownToBeHomogeneous)
{
  // New optimised version using PolyhedralCone
  int n=groebnerBasis.numberOfVariablesInRing();
  IntegerVectorList a;
  //AsciiPrinter(Stderr).printVectorList(wallInequalities(groebnerBasis));
  //PolyhedralCone p=intersection(PolyhedralCone(wallInequalities(groebnerBasis),a),PolyhedralCone::positiveOrthant(n));

  IntegerVectorList normals=algebraicTest(wallInequalities(groebnerBasis),groebnerBasis);
    coneTimer.on();
    //    PolyhedralCone p=PolyhedralCone(wallInequalities(groebnerBasis),a);

    /*    AsciiPrinter(Stderr).printVectorList(normals);
    AsciiPrinter(Stderr).printInteger(normals.size());
    */
    normals=normalizedWithSumsAndDuplicatesRemoved(normals);

    //    AsciiPrinter(Stderr).printVectorList(normals);
    //   AsciiPrinter(Stderr).printInteger(normals.size());

    PolyhedralCone p=PolyhedralCone(normals,a,n);
    //  fprintf(Stderr,"Finding facets\n");
    p.findFacets();
    //  fprintf(Stderr,"Done Finding facets\n");
    coneTimer.off();

  IntegerVectorList ilist=p.getHalfSpaces();
  IntegerVectorList ret;
  for(IntegerVectorList::iterator i=ilist.begin();i!=ilist.end();i++)
    {
      bool facetOK=true;
      if(!isKnownToBeHomogeneous)
	{
	  *i=-1 * (*i);
	  PolyhedralCone temp=intersection(PolyhedralCone(ilist,a),PolyhedralCone::positiveOrthant(n));
	  facetOK=(temp.dimension()==n);
	  *i=-1 * (*i);
	}

      if(facetOK)
	ret.push_back(*i);
    }

  // New version using PolyhedralCone
  /*  int n=groebnerBasis.numberOfVariablesInRing();
  IntegerVectorList a;
  //  AsciiPrinter(Stderr).printVectorList(wallInequalities(groebnerBasis));
  PolyhedralCone p=intersection(PolyhedralCone(wallInequalities(groebnerBasis),a),PolyhedralCone::positiveOrthant(n));
  //  PolyhedralCone p=PolyhedralCone(wallInequalities(groebnerBasis),a);
  p.findFacets();
  IntegerVectorList ilist=p.getHalfSpaces();
  IntegerVectorList ret;
  for(IntegerVectorList::const_iterator i=ilist.begin();i!=ilist.end();i++)
    if(wallContainsPositiveVector(*i))ret.push_back(*i);
  */

  /* Old code not using PolyhedralCone
  IntegerVectorList ilist=wallInequalities(groebnerBasis);
  wallAddCoordinateWalls(ilist);
  IntegerVectorList ret;
  for(IntegerVectorList::const_iterator i=ilist.begin();i!=ilist.end();i++)
    if(isFacet(ilist,i)&&wallContainsPositiveVector(*i))
      ret.push_back(*i);
  */
  return ret;
}


Polynomial flipMinkowski(Polynomial p, IntegerVector const &wallNormal)
{
  IntegerVector best=p.getMarked().m.exponent;

  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      IntegerVector e=i->first.exponent;
      IntegerVector diff=e-best;
      if(dependent(diff,wallNormal)&&dot(diff,wallNormal)<0)best=e;
    }
  p.mark(Monomial(p.getRing(),best));

  return p;
}


PolynomialSet flipMinkowski(PolynomialSet const &groebnerBasis, IntegerVector const &wallNormal)
{
  PolynomialSet r(groebnerBasis.getRing());

  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    r.push_back(flipMinkowski(*i,wallNormal));

  return r;
}


PolyhedralCone homogeneitySpace(PolynomialSet const &reducedGroebnerBasis)
{
  IntegerVectorList l=wallInequalities(reducedGroebnerBasis);
  IntegerVectorList a;
  PolyhedralCone c(a,l,reducedGroebnerBasis.getRing().getNumberOfVariables());
  c.findImpliedEquations();
//  c.canonicalize();
  return c;
}

PolyhedralCone groebnerCone(PolynomialSet const &reducedGroebnerBasis, bool useAlgebraicTest)
{
  int n=reducedGroebnerBasis.getRing().getNumberOfVariables();
  IntegerVectorList l=wallInequalities(reducedGroebnerBasis);
  if(useAlgebraicTest)l=algebraicTest(l,reducedGroebnerBasis);
  l=normalizedWithSumsAndDuplicatesRemoved(l);
  IntegerVectorList a;
  PolyhedralCone c(l,a,n);
  c.canonicalize();
  return c;
}

int dimensionOfHomogeneitySpace(PolynomialSet const &reducedGroebnerBasis)
{
  return homogeneitySpace(reducedGroebnerBasis).dimension();
}


PolynomialSet liftBasis(PolynomialSet const &toBeLifted, PolynomialSet const &originalBasisForFullIdeal)
{
  PolynomialRing theRing=toBeLifted.getRing();
  assert(toBeLifted.isValid());
  assert(originalBasisForFullIdeal.isValid());

  PolynomialSet newBasis(theRing);

  //  fprintf(Stderr,"LIFTING:");
  //  AsciiPrinter(Stderr).printPolynomialSet(toBeLifted);

  for(PolynomialSet::const_iterator j=toBeLifted.begin();j!=toBeLifted.end();j++)
    {
      assert(!j->isZero());
      //AsciiPrinter(Stderr).printVector(j->getMarked().m.exponent);

      Polynomial q=*j-division(*j,originalBasisForFullIdeal,LexicographicTermOrder());
      assert(!q.isZero());
      q.mark(j->getMarked().m);
      newBasis.push_back(q);
    }
  autoReduce(&newBasis,LexicographicTermOrder());
  //  fprintf(Stderr,"TO:");
  //  AsciiPrinter(Stderr).printPolynomialSet(newBasis);

  return newBasis;
}


bool isMarkingConsistent(PolynomialSet const &g)
{
  IntegerVectorList empty;
  PolyhedralCone c(wallInequalities(g),empty,g.getRing().getNumberOfVariables());
  c=intersection(c,PolyhedralCone::positiveOrthant(c.ambientDimension()));
  log1 AsciiPrinter(Stderr).printPolyhedralCone(c);
  return c.dimension()==c.ambientDimension();
}


/*
  Heuristic for checking if inequality of full dimensional cone is a
  facet. If the routine returns true then the inequality is a
  facet. If it returns false it is unknown.
 */
static bool fastIsFacetCriterion(IntegerVectorList const &normals, IntegerVectorList::const_iterator i)
{
  int n=i->size();
  for(int j=0;j<n;j++)
    if((*i)[j]>0 || (*i)[j]<0)
      {
	int sign=(*i)[j]>0?1:-1;
	bool isTheOnly=true;
	for(IntegerVectorList::const_iterator k=normals.begin();k!=normals.end();k++)
	  if(k!=i)
	    {
	      if((*k)[j]*sign>0)
		{
		  isTheOnly=false;
		  break;
		}
	    }
	if(isTheOnly)return true;
      }
  return false;
}

bool fastIsFacet(IntegerVectorList const &normals, IntegerVectorList::const_iterator i)
{
  if(fastIsFacetCriterion(normals,i))return true;
  //  log0 fprintf(Stderr,"LP\n");
  return isFacet(normals,i);
}


/*
ONLY WORKS AFFINELY AND NOT PROJECTIVELY!!!
bool fastIsFacet(IntegerVectorList const &normals, IntegerVectorList::const_iterator i)
{
  int n=i->size();

  bool rightAnswer=isFacet(normals,i);

  IntegerVectorList tempNormals=normals;
  bool doLoop=true;
  log0 AsciiPrinter(Stderr).printVector(*i);
  do
    {
      log0 fprintf(Stderr,"TempNormal size:%i\n",tempNormals.size());
      log0 AsciiPrinter(Stderr).printVectorList(tempNormals);
      IntegerVector maxVector=*i;
      IntegerVector minVector=*i;
      for(IntegerVectorList::const_iterator k=tempNormals.begin();k!=tempNormals.end();k++)
	{
	  maxVector=max(maxVector,*k);
	  minVector=min(minVector,*k);
	}
      IntegerVector maxAttained(n);
      IntegerVector minAttained(n);
      for(IntegerVectorList::const_iterator k=tempNormals.begin();k!=tempNormals.end();k++)
	{
	  for(int j=0;j<n;j++)
	    {
	      if((*k)[j]==maxVector[j])maxAttained[j]++;
	      if((*k)[j]==minVector[j])minAttained[j]++;
	    }
	}
      int bestEntry=-1;
      int bestCount=2000000000;
      int bestValue=0;
      for(int j=0;j<n;j++)
	{
	  if((*i)[j]==maxVector[j])
	    {
	      if(maxAttained[j]<bestCount)
		{
		  bestEntry=j;
		  bestValue=maxVector[j];
		  bestCount=maxAttained[j];
		}
	    }
	  if((*i)[j]==minVector[j])
	    {
	      if(minAttained[j]<bestCount)
		{
		  bestEntry=j;
		  bestValue=minVector[j];
		  bestCount=minAttained[j];
		}
	    }
	}
      log0 fprintf(Stderr,"Best entry:%i, Best Count:%i bestValue : %i\n",bestEntry,bestCount,bestValue);

      IntegerVectorList newList;
      for(IntegerVectorList::const_iterator k=tempNormals.begin();k!=tempNormals.end();k++)
	{
	  if((*k)[bestEntry]==bestValue)newList.push_back(*k);
	}
      if(newList.size()==1)
	{
	  assert(rightAnswer==true);
	  return true;
	}
      doLoop=tempNormals.size()>newList.size();
      tempNormals=newList;
    }
  while(doLoop);


  log0 fprintf(Stderr,"LP\n");
  return isFacet(normals,i);
}
*/
 /*
bool fastIsFacet(IntegerVectorList const &normals, IntegerVectorList::const_iterator i)
{
  int n=i->size();

  bool rightAnswer=isFacet(normals,i);

  IntegerVectorList tempNormals=normals;
  bool doLoop=true;
  //log0 AsciiPrinter(Stderr).printVector(*i);
  do
    {
      //      log0 fprintf(Stderr,"TempNormal size:%i\n",tempNormals.size());
      //log0 AsciiPrinter(Stderr).printVectorList(tempNormals);
      IntegerVector maxVector=*i;
      IntegerVector minVector=*i;
      for(IntegerVectorList::const_iterator k=tempNormals.begin();k!=tempNormals.end();k++)
	{
	  maxVector=max(maxVector,*k);
	  minVector=min(minVector,*k);
	}



      /*      IntegerVector maxAttained(n);
      IntegerVector minAttained(n);
      for(IntegerVectorList::const_iterator k=tempNormals.begin();k!=tempNormals.end();k++)
	{
	  for(int j=0;j<n;j++)
	    {
	      if((*k)[j]==maxVector[j])maxAttained[j]++;
	      if((*k)[j]==minVector[j])minAttained[j]++;
	    }
	}
      */
/*    int bestEntry=-1;
      int bestCount=2000000000;
      int bestValue=0;
      for(int j=0;j<n;j++)
	{
	  if((*i)[j]==0)
	    {
	      if(maxVector[j]==0)
	    }
	  else
	    {
	      int sign=(*i)[j]>0?1:-1;

	      int numberOfVectorsWithSameSign=0;
	      for(IntegerVectorList::const_iterator k=normals.begin();k!=normals.end();k++)
		{
		  if((*k)[j]*sign>0)
		    {
		      numberOfVectorsWithSameSign++;
		    }
		}
	      if(numberOfVectorsWithSameSign==1)
		return true;
	    }


	  if((*i)[j]==maxVector[j])
	    {
	      if(maxAttained[j]<bestCount)
		{
		  bestEntry=j;
		  bestValue=maxVector[j];
		  bestCount=maxAttained[j];
		}
	    }
	  if((*i)[j]==minVector[j])
	    {
	      if(minAttained[j]<bestCount)
		{
		  bestEntry=j;
		  bestValue=minVector[j];
		  bestCount=minAttained[j];
		}
	    }
	}
      //log0 fprintf(Stderr,"Best entry:%i, Best Count:%i bestValue : %i\n",bestEntry,bestCount,bestValue);

      IntegerVectorList newList;
      for(IntegerVectorList::const_iterator k=tempNormals.begin();k!=tempNormals.end();k++)
	{
	  if((*k)[bestEntry]==bestValue)newList.push_back(*k);
	}
      if(newList.size()==1)
	{
	  assert(rightAnswer==true);
	  return true;
	}
      doLoop=tempNormals.size()>newList.size();
      tempNormals=newList;
    }
  while(doLoop);


  //  log0 fprintf(Stderr,"LP\n");
  return isFacet(normals,i);
}
*/

IntegerVectorList fastNormals(IntegerVectorList const &inequalities)
{
//    log0 cerr << "Fast normals begin" << endl;
  //  log0 fprintf(Stderr,"E");
  //log0 fprintf(Stderr,"Number of inequalities:%i\n",inequalities.size());
  IntegerVectorList normals=normalizedWithSumsAndDuplicatesRemoved(inequalities);
   // log0 fprintf(Stderr,"Number of inequalities:%i\n",normals.size());
  //  log0 fprintf(Stderr,"F");


  //     log0 AsciiPrinter(Stderr).printVectorList(normals);
  for(IntegerVectorList::iterator i=normals.begin();i!=normals.end();i++)
    //if(!isFacet(normals,i))
    if(!fastIsFacet(normals,i))
      {
	IntegerVectorList::iterator temp=i;
	temp++;
	normals.erase(i);
	temp--;
	i=temp;
      }

  //  log0 fprintf(Stderr,"Number of inequalities:%i\n",normals.size());

  //  log0 fprintf(Stderr,"G");
  //log2 cerr << "Fast normals end" << endl;
  return normals;
}

/* Virker ikke:
IntegerVectorList fastNormals(IntegerVectorList const &inequalities)
{
  log0 fprintf(Stderr,"E");
   log0 fprintf(Stderr,"Number of inequalities:%i\n",inequalities.size());
  IntegerVectorList normals=normalizedWithSumsAndDuplicatesRemoved(inequalities);
  log0 fprintf(Stderr,"Number of inequalities:%i\n",normals.size());
  log0 fprintf(Stderr,"F");


  //     log0 AsciiPrinter(Stderr).printVectorList(normals);

  bool didRemoveSomething;
  do
    {
      didRemoveSomething=false;
      for(IntegerVectorList::iterator i=normals.begin();i!=normals.end();i++)
	//if(!isFacet(normals,i))
	if(!fastIsFacetCriterion(normals,i))
	  {
	    IntegerVectorList::iterator temp=i;
	    temp++;
	    normals.erase(i);
	    temp--;
	    i=temp;
	    didRemoveSomething=true;
	  }
    }
  while(didRemoveSomething);
  for(IntegerVectorList::iterator i=normals.begin();i!=normals.end();i++)
    {
      if(!fastIsFacet(normals,i))
	{
	  IntegerVectorList::iterator temp=i;
	  temp++;
	  normals.erase(i);
	  temp--;
	  i=temp;
	  didRemoveSomething=true;
	}
    }
   log0 fprintf(Stderr,"Number of inequalities:%i\n",normals.size());

  log0 fprintf(Stderr,"G");
  return normals;
}
*/


PolyhedralCone normalConeOfMinkowskiSum(PolynomialSet const &polynomials, IntegerVector const &w)
{
  //  log0 cerr << "A";
  WeightReverseLexicographicTermOrder T(w);
  PolynomialSet g=polynomials;
  //log0 cerr << "I";
  g.markAndScale(T);
  //log0 cerr << "I";
  IntegerVectorList inequalities=fastNormals(wallInequalities(g));
  g=initialFormsAssumeMarked(g,w);
  IntegerVectorList equations=wallInequalities(g);
  //log0 cerr << "B";
  PolyhedralCone ret(inequalities,equations,g.getRing().getNumberOfVariables());
  //log0 cerr << "C";
  ret.canonicalize();
  //log0 cerr << "D";
  return ret;
}


IntegerVectorList commonHomogeneitySpaceEquations(PolynomialSet const &polynomials)
{
  LexicographicTermOrder T;
  PolynomialSet g=polynomials;
  g.markAndScale(T);
  IntegerVectorList l=wallInequalities(g);
  FieldMatrix A=integerMatrixToFieldMatrix(rowsToIntegerMatrix(l,g.getRing().getNumberOfVariables()),Q);
  A.reduce();
  A.removeZeroRows();
  return fieldMatrixToIntegerMatrixPrimitive(A).getRows();
}


IntegerVectorList commonHomogeneitySpaceGenerators(PolynomialSet const &polynomials)
{
  LexicographicTermOrder T;
  PolynomialSet g=polynomials;
  g.markAndScale(T);
  IntegerVectorList l=wallInequalities(g);
  FieldMatrix A=integerMatrixToFieldMatrix(rowsToIntegerMatrix(l,g.getRing().getNumberOfVariables()),Q);
  A=A.reduceAndComputeKernel();
  return fieldMatrixToIntegerMatrixPrimitive(A).getRows();
}
