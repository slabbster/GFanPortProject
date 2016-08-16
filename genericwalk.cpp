#include "genericwalk.h"

#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "printer.h"

//-----------------------------------------------------------------
// Generic Groebner Walk
//-----------------------------------------------------------------

static bool isCandidate(IntegerVector const &v, const TermOrder &source, const TermOrder &target)
{
  return !source(v,v,1,0) && target(v,v,1,0);
}

static bool compareInequalities(IntegerVector const &a, IntegerVector const &b, TermOrder const &source, TermOrder const &target)
{
  int n=a.size();

  for(int i=0;i<n;i++)
    {
      int tdb=target.rowDot(i,b);
      int tda=target.rowDot(i,a);
      if(tdb==0 && tda==0)continue;
      if(tda==0)
	{
	  if(tdb>0)return true;
	  if(tdb<0)return false;
	  continue;
	}
      if(tdb==0)
	{
	  if(tda<0)return true;
	  if(tda>0)return false;
	  continue;
	}
      if(source(b,a,tda,tdb))return true;
      if(source(a,b,tdb,tda))return false; 
    }
  return false;
}

IntegerVectorList::const_iterator shootGenericRay(IntegerVectorList const &g, const TermOrder &source, const TermOrder &target)
{
  IntegerVectorList::const_iterator candidate=g.end();
  
  for(IntegerVectorList::const_iterator i=g.begin();i!=g.end();i++)
    if(isCandidate(*i,source,target))
      {
	if(candidate==g.end()||compareInequalities(*i,*candidate,source,target))
	  candidate=i;
      }  
  
  return candidate;
}

PolynomialSet genericWalk(PolynomialSet const &start, const TermOrder &source, const TermOrder &target)
{
  int nflips=0;
  PolynomialSet g=start;

  assert(start.checkMarkings(source));

  while(1)
  {
    IntegerVectorList w=wallInequalities(g);
    IntegerVectorList::const_iterator i=shootGenericRay(w,source,target);
    if(i==w.end())break;
    g=flip(g,*i);

    nflips++;
  }
  fprintf(Stderr,"Number of flips:%i\n",nflips);
  return g;
}


class PreOrder
{
  TermOrder const &source;
  TermOrder const &target;
  int sourceDegree;
  int targetDegree;
public:
  PreOrder(TermOrder const &source_, TermOrder const &target_, int sourceDegree_, int targetDegree_):
    source(source_),
    target(target_),
    sourceDegree(sourceDegree_),
    targetDegree(targetDegree_)
  {
  }
  bool isCandidate(IntegerVector const &v)
  {
    /*    AsciiPrinter(Stderr).printVector(v);
    AsciiPrinter(Stderr).printNewLine();
    AsciiPrinter(Stderr).printTermOrder(source);
    AsciiPrinter(Stderr).printTermOrder(target);*/
    return !source(v,v,1,0) && target(v,v,1,0); // Should this depend on the perturbation degree??
  }
  bool operator()(IntegerVector const &a, IntegerVector const &b)
  {
    int n=a.size();

    for(int i=0;i<targetDegree;i++)
      {
	int tdb=target.rowDot(i,b);
	int tda=target.rowDot(i,a);
	if(tdb==0 && tda==0)continue;
	if(tda==0)
	  {
	    if(tdb>0)return true;
	    if(tdb<0)return false;
	    continue;
	  }
	if(tdb==0)
	  {
	    if(tda<0)return true;
	    if(tda>0)return false;
	    continue;
	  }
	if(source(b,a,tda,tdb,sourceDegree))return true;
	if(source(a,b,tdb,tda,sourceDegree))return false; 
      }
    return false;
  }
};


PolynomialSet genericWalkPerturbation(PolynomialSet const &start, const TermOrder &source, const TermOrder &target, int sourceDegree, int targetDegree)
{
  PolynomialRing theRing=start.getRing();
  PreOrder p(source,target,sourceDegree,targetDegree);
  int nflips=0;
  PolynomialSet g=start;

  assert(start.checkMarkings(source));

  while(1)
  {

    IntegerVectorList inequalities=wallInequalities(g);
    assert(!inequalities.empty());
    IntegerVectorList::const_iterator optimal=inequalities.end();

    for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
      {
	//	fprintf(Stderr,"%i %i\n",!source(*i,*i,1,0),target(*i,*i,1,0));
	if(p.isCandidate(*i))
	  {
	    if(optimal==inequalities.end() || p(*i,*optimal))
	      optimal=i;
	  }
      }
    
    if(optimal==inequalities.end())
      break;
    else
      AsciiPrinter(Stderr).printVector(*optimal);

    PolynomialSet faceIdeal(theRing);
    for(PolynomialSet::const_iterator j=g.begin();j!=g.end();j++)
      {
	IntegerVector markedExponent=j->getMarked().m.exponent;
	Polynomial r(theRing);
	for(TermMap::const_iterator i=j->terms.begin();i!=j->terms.end();i++)
	  {
	    IntegerVector dif=markedExponent-i->first.exponent;
	    if( ((!p(dif,*optimal))&&(!p(*optimal,dif)))|| dif.isZero())
	      r+=Polynomial(Term(i->second,i->first));
	  }	
	faceIdeal.push_back(r);
      }

    
    faceIdeal.markAndScale(source);
    //    AsciiPrinter(Stderr).printPolynomialSet(faceIdeal);

    PolynomialSet oldFaceIdeal=faceIdeal;
    buchberger(&faceIdeal,target);

    
    fprintf(Stderr,"Lifting\n");
    PolynomialSet newG(theRing);
    for(PolynomialSet::const_iterator j=faceIdeal.begin();j!=faceIdeal.end();j++)
      {
	newG.push_back(divisionLift(*j, oldFaceIdeal, g, LexicographicTermOrder()));
      }

    {
      PolynomialSet::const_iterator k=faceIdeal.begin();
      for(PolynomialSet::iterator j=newG.begin();j!=newG.end();j++)
	{
	  j->mark(k->getMarked().m);
	  k++;
	}
    }     

    fprintf(Stderr,"Autoreducing\n");
    autoReduce(&newG,StandardGradedLexicographicTermOrder());

    g=newG;

    nflips++;
    fprintf(Stderr,"Flip %i, new size %i\n",nflips,g.size());
  }
  fprintf(Stderr,"Number of flips:%i\n",nflips);

  autoReduce(&g,StandardGradedLexicographicTermOrder());

  return g;
}
