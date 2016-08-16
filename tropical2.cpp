#include "tropical2.h"

#include <stdlib.h>
#include <iostream>

#include "buchberger.h"
#include "division.h"
#include "tropical.h"
#include "wallideal.h"
#include "dimension.h"
#include "halfopencone.h"
#include "breadthfirstsearch.h"

#include "timer.h"
#include "log.h"

static Timer tropicalPrincipalIntersectionTimer("Tropical principal intersection",1);

static void startingConeError()
{
  fprintf(Stderr,"UNABLE TO COMPUTE STARTING CONE.\n");
  fprintf(Stderr,"THE STARTING CONE ALGORITHM IN GFAN IS BASED ON HEURISTICS WHICH HAVE FAILED ON THIS EXAMPLE.\n");
  assert(0);
}

PolynomialSet initialIdeal(PolynomialSet const &g, IntegerVector const &weight)
//Assume homogeneous
{
  PolynomialSet ret=g;
  WeightReverseLexicographicTermOrder T(weight);
  buchberger(&ret,T);
  return initialForms(ret,weight);
}


Polynomial initialFormAssumeMarked(Polynomial const &p, IntegerVector const &weight)
{
  Polynomial r(p.getRing());
  IntegerVector markedExponent=p.getMarked().m.exponent;

  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      IntegerVector dif=markedExponent-i->first.exponent;

      if(dot(dif,weight)==0)
	r+=Polynomial(Term(i->second,i->first));
    }
  r.mark(Monomial(p.getRing(),markedExponent));

  return r;
}


PolynomialSet initialFormsAssumeMarked(PolynomialSet const &groebnerBasis, IntegerVector const &weight)
{
  PolynomialRing theRing=groebnerBasis.getRing();
  PolynomialSet r(theRing);

  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    {
      r.push_back(initialFormAssumeMarked(*i,weight));
    }
  return r;
}


Polynomial initialForm(Polynomial const &p, IntegerVector const &weight)
{
  if(p.isZero())return p;
  int64 a=dotLong(p.terms.begin()->first.exponent,weight);
  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      int64 b=dotLong(i->first.exponent,weight);
      if(b>a)a=b;
    }

  Polynomial r(p.getRing());
  bool ismarked=p.isMarked();
  IntegerVector markedExponent;
  if(ismarked)markedExponent=p.getMarked().m.exponent;

  bool markedFound=false;

  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      if(dotLong(i->first.exponent,weight)==a)
	{
	  r+=Polynomial(Term(i->second,i->first));
	  if(ismarked)if((markedExponent-(i->first.exponent)).isZero())markedFound=true;
	}
    }
  if(markedFound)
    r.mark(Monomial(p.getRing(),markedExponent));

  return r;
}


PolynomialSet initialForms(PolynomialSet const &groebnerBasis, IntegerVector const &weight)
{
  PolynomialRing theRing=groebnerBasis.getRing();
  PolynomialSet r(theRing);
  if(theRing.getNumberOfVariables()!=weight.size())
    {
      cerr << "Error: Number of varaibles in polynomial ring "<<theRing.getNumberOfVariables()<< " length of weight vector " << weight.size() <<endl;
      assert(0);
    }

  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    {
      r.push_back(initialForm(*i,weight));
    }
  return r;
}



PolyhedralFan tropicalPrincipalIntersection(int n, PolynomialSet const &g, int linealitySpaceDimension)
{
  //return tropicalHyperSurfaceIntersection(n, g);////////////////////////////////////////
  log2 fprintf(Stderr,"Intersecting\n");
  log3 AsciiPrinter(Stderr).printPolynomialSet(g);

  TimerScope ts(&tropicalPrincipalIntersectionTimer);
  PolyhedralFan ret=PolyhedralFan::fullSpace(n);

  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    {
      ret=refinement(ret,PolyhedralFan::bergmanOfPrincipalIdeal(*i),linealitySpaceDimension,true);
    }
  log2 fprintf(Stderr,"Done intersecting\n");
  return ret;
}



static PolynomialSet checkList(IntegerVectorList const &l, PolynomialSet const &groebnerBasis, PolynomialSet *fullNeighbourBasis, int h, bool &result, bool onlyCheckRays)
{
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      WeightReverseLexicographicTermOrder t(*i);
      log2 fprintf(Stderr,"Computing Gr\"obner basis with respect to:");
      log2 AsciiPrinter(Stderr).printVector(*i);
      log2 fprintf(Stderr,"\n");
      PolynomialSet h2=groebnerBasis;
      buchberger(&h2,t);
      log2 fprintf(Stderr,"Done computing Gr\"obner basis.\n");

      log3 AsciiPrinter(Stderr).printPolynomialSet(h2);
      PolynomialSet wall=initialFormsAssumeMarked(h2,*i);

      log3 AsciiPrinter(Stderr).printString("Initial ideal:\n");
      log3 AsciiPrinter(Stderr).printPolynomialSet(wall);

      int hdim2=dimensionOfHomogeneitySpace(wall);
      if(hdim2>h)
	{
	  if(!containsMonomial(wall))
	    {
	      log1 fprintf(Stderr,"Iterating recursively.\n");
	      //PolynomialSet initialIdeal=guessInitialIdealWithoutMonomial(wall,0);
	      PolynomialSet initialIdeal=guessInitialIdealWithoutMonomial(wall,fullNeighbourBasis,onlyCheckRays);

	      if(fullNeighbourBasis)
		{
		  //*fullNeighbourBasis=liftBasis(initialIdeal,h2);
		  *fullNeighbourBasis=liftBasis(*fullNeighbourBasis,h2);
		}


	      result=true;
	      return initialIdeal;
	    }
	}
    }
  result=false;
  return groebnerBasis;
}

PolynomialSet guessInitialIdealWithoutMonomial(PolynomialSet const &groebnerBasis, PolynomialSet *fullNeighbourBasis, bool onlyCheckRays) //ideal must be homogeneous
  // fullNeighbourBasis is set to a Groebner basis of the full ideal. The returned basis and fullNeighbourBasis have at least one termorder in common
{
  //  log0 fprintf(Stderr,"A\n");
  assert(groebnerBasis.isValid());
  //  log0 fprintf(Stderr,"B\n");
  if(fullNeighbourBasis)
    {
      assert(fullNeighbourBasis->isValid());
    }
  //  log0 fprintf(Stderr,"C\n");

  int n=groebnerBasis.numberOfVariablesInRing();
  //  log0 fprintf(Stderr,"D\n");
  int h=dimensionOfHomogeneitySpace(groebnerBasis);
  //  log0 fprintf(Stderr,"E\n");
  int d=krullDimension(groebnerBasis);
  //  log0 fprintf(Stderr,"F\n");

  if(d==h)
    {
      if(fullNeighbourBasis)*fullNeighbourBasis=groebnerBasis;
      return groebnerBasis;
    }

  {
    log2 fprintf(Stderr,"Computing extreme rays.\n");
    //IntegerVectorList a;
    PolyhedralCone p=coneFromMarkedBasis(groebnerBasis);
    //PolyhedralCone p=PolyhedralCone(wallInequalities(groebnerBasis),a);
    IntegerVectorList extreme=p.extremeRays();
    log2 fprintf(Stderr,"Extreme rays of Groebner cone:\n");
    log2 AsciiPrinter(Stderr).printVectorList(extreme);

    bool result;
    PolynomialSet r=checkList(extreme,groebnerBasis,fullNeighbourBasis,h,result, onlyCheckRays);
    if(result)return r;
  }
  if(onlyCheckRays)startingConeError();
  PolyhedralFan f=PolyhedralFan::fullSpace(n);
  /*  for(int i=0;i<d-1;i++)
    {
      IntegerVector v(n);
      for(int j=0;j<n;j++)v[j]=rand()&1;
      IntegerVectorList a,b;
      b.push_back(v);
      PolyhedralFan F(n);
      F.insert(PolyhedralCone(a,b));
      f=refinement(f,F);
    }
  AsciiPrinter P(Stderr);
  f.print(&P);
  */
  int hypersurfacesToGo=groebnerBasis.size();
  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    {
      fprintf(Stderr,"Hypersurfaces to go:%i\n",hypersurfacesToGo--);
      fprintf(Stderr,"Max dimension: %i\n",f.getMaxDimension());
      f=refinement(f,PolyhedralFan::bergmanOfPrincipalIdeal(*i));
      f.removeAllExcept(3);

      IntegerVectorList l=f.getRelativeInteriorPoints();

      bool result;
      PolynomialSet r=checkList(l,groebnerBasis,fullNeighbourBasis,h,result, onlyCheckRays);
      if(result)return r;
    }
  startingConeError();
  return groebnerBasis;
}

static PolynomialSet checkListStably(IntegerVectorList const &l, PolynomialSet const &groebnerBasis, PolynomialSet *fullNeighbourBasis, int h, bool &result, bool onlyCheckRays)
{
  debug<< "Checklist called on"<<groebnerBasis;
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      WeightReverseLexicographicTermOrder t(*i);
      log2 fprintf(Stderr,"Taking initial forms with respect to:");
      log2 AsciiPrinter(Stderr).printVector(*i);
      log2 fprintf(Stderr,"\n");
      PolynomialSet h2=groebnerBasis;
      log2 fprintf(Stderr,"Done computing Gr\"obner basis.\n");

      log3 AsciiPrinter(Stderr).printPolynomialSet(h2);
      PolynomialSet wall=initialForms(h2,*i);

      log3 AsciiPrinter(Stderr).printString("Initial ideal:\n");
      log3 AsciiPrinter(Stderr).printPolynomialSet(wall);

      int hdim2=dimensionOfHomogeneitySpace(wall);
      if(hdim2>h)
        {
          if(nonEmptyStableIntersection(wall))
            {
              log1 fprintf(Stderr,"Iterating recursively.\n");
              //PolynomialSet initialIdeal=guessInitialIdealWithoutMonomial(wall,0);
              PolynomialSet initialIdeal=guessInitialIdealWithoutMonomialStably(wall,fullNeighbourBasis,onlyCheckRays);

              if(fullNeighbourBasis)
                {
                  //*fullNeighbourBasis=liftBasis(initialIdeal,h2);
//                  *fullNeighbourBasis=liftBasis(*fullNeighbourBasis,h2);
                  *fullNeighbourBasis=groebnerBasis;
                  fullNeighbourBasis->copyMarkings(initialIdeal);
                }


              result=true;
              return initialIdeal;
            }
        }
    }
  result=false;
  return groebnerBasis;
}

PolynomialSet guessInitialIdealWithoutMonomialStably(PolynomialSet const &groebnerBasis, PolynomialSet *fullNeighbourBasis, bool onlyCheckRays) //ideal must be homogeneous
  // fullNeighbourBasis is set to a Groebner basis of the full ideal. The returned basis and fullNeighbourBasis have at least one termorder in common
{
  int n=groebnerBasis.numberOfVariablesInRing();
  int h=dimensionOfHomogeneitySpace(groebnerBasis);
  int d=n-groebnerBasis.size();//krullDimension(groebnerBasis);

  debug<</*"d"<<d<<*/"h"<<h<<"n"<<n<<"\n";


  if(d==h)
    {
      if(fullNeighbourBasis)*fullNeighbourBasis=groebnerBasis;
      return groebnerBasis;
    }

  {
    log2 fprintf(Stderr,"Computing extreme rays.\n");
    //IntegerVectorList a;
    PolyhedralCone p=coneFromMarkedBasis(groebnerBasis);
    //PolyhedralCone p=PolyhedralCone(wallInequalities(groebnerBasis),a);
    IntegerVectorList extreme=p.extremeRays();
    log2 fprintf(Stderr,"Extreme rays of Groebner cone:\n");
    log2 AsciiPrinter(Stderr).printVectorList(extreme);

    bool result;
    PolynomialSet r=checkListStably(extreme,groebnerBasis,fullNeighbourBasis,h,result, onlyCheckRays);
    if(result)return r;
  }
  if(onlyCheckRays)startingConeError();

  PolyhedralFan f=PolyhedralFan::fullSpace(n);

  int hypersurfacesToGo=groebnerBasis.size();
  for(PolynomialSet::const_iterator i=groebnerBasis.begin();i!=groebnerBasis.end();i++)
    {
      fprintf(Stderr,"Hypersurfaces to go:%i\n",hypersurfacesToGo--);
      fprintf(Stderr,"Max dimension: %i\n",f.getMaxDimension());
      f=refinement(f,PolyhedralFan::bergmanOfPrincipalIdeal(*i));
      f.removeAllExcept(3);

      IntegerVectorList l=f.getRelativeInteriorPoints();

      bool result;
      PolynomialSet r=checkListStably(l,groebnerBasis,fullNeighbourBasis,h,result, onlyCheckRays);
      if(result)return r;
    }
  startingConeError();
  return groebnerBasis;
}

