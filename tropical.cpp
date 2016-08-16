#include "tropical.h"
#include "vektor.h"
#include "lp.h"
#include "wallideal.h"
#include "printer.h"
#include "buchberger.h"
#include "division.h"
#include "matrix.h"
#include "linalg.h"
#include "log.h"

bool isFullColored(IntegerVectorList const &inequalityColors, IntegerVector const &v)
{
  assert(inequalityColors.size()==v.size());

  int nColors=inequalityColors.begin()->size();
  for(int i=0;i<nColors;i++)
    {
      bool colored=false;
      int j=0;
      for(IntegerVectorList::const_iterator J=inequalityColors.begin();J!=inequalityColors.end();J++)
	{
	  if(v[j])
	    {
	      if((*J)[i])colored=true;
	    }

	  j++;
	}
      if(!colored)return false;
    }
  return true;
}

static bool increase(IntegerVector &v)
{
  int i=0;
  while(i<v.size() && v[i]==1)
    {
      v[i]=0;
      i++;
    }
  if(i==v.size())return false;
  v[i]=1;
  return true;
}

static bool increaseSkip(IntegerVector &v)
{
  int i=0;
  while(i<v.size() && v[i]==0)
    {
      i++;
    }
  if(i==v.size())return false;
  while(i<v.size() && v[i]==1)
    {
      v[i]=0;
      i++;
    }
  if(i==v.size())return false;
  v[i]=1;
  return true;
}

PolynomialSetList fullColoredIdeals(PolynomialSet const &g, bool skipColorTest)
{
  // by computing the cone faces in the most stupid way
  PolynomialSetList ret;
  IntegerVectorList inequalities=wallInequalities(g);

  {
    IntegerVectorList normals=wallRemoveScaledInequalities(inequalities);
    IntegerVectorList facets;
    for(IntegerVectorList::const_iterator i=normals.begin();i!=normals.end();i++)
      //if(wallContainsPositiveVector(*i))
      if(isFacet(normals,i))
        {
	  facets.push_back(*i);
        }
    inequalities=facets;
  }
  fprintf(Stderr,"INEQUALITIES\n");
  AsciiPrinter(Stderr).printVectorList(inequalities);



  assert(!inequalities.empty());

  IntegerVectorList inequalityColors;

  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
    {
      IntegerVector colors(g.size());
      int j=0;
      for(PolynomialSet::const_iterator J=g.begin();J!=g.end();J++)
	{
	  PolynomialSet temp(g.getRing());
	  temp.push_back(*J);
	  IntegerVectorList tempInequalities=wallInequalities(temp);
	  bool colored=false;

	  for(IntegerVectorList::const_iterator k=tempInequalities.begin();k!=tempInequalities.end();k++)
	    {
	      if(dependent(*k,*i))
		{
		  colored=true;
		  break;
		}
	    }
	  colors[j]=colored;
	  j++;
	}
      //      AsciiPrinter(Stderr).printVector(colors);
      inequalityColors.push_back(colors);
    }

  //  AsciiPrinter(Stderr).printVectorList(inequalityColors);

  IntegerVector equalitySet(inequalities.size());
  bool skip;
  do
    {
      AsciiPrinter(Stderr).printVector(equalitySet);fprintf(Stderr,"\n");
      skip=false;
      if(skipColorTest || isFullColored(inequalityColors,equalitySet))
	{
	  if(hasInteriorPoint(inequalities,true,&equalitySet))
	    {
	      fprintf(Stderr,"Adding: ");
	      AsciiPrinter(Stderr).printVector(equalitySet);

	      IntegerVectorList es;
	      int i=0;
	      for(IntegerVectorList::const_iterator I=inequalities.begin();I!=inequalities.end();I++)
		{
		  if(equalitySet[i])es.push_back(*I);
		  i++;
		}
	      PolynomialSet w=lowerDimensionalWallIdeal(g,es);
	      //AsciiPrinter(Stderr).printPolynomialSet(w);
	      ret.push_back(w);
	    }
	  else
	    skip=true;
	}
    }
  while(increase(equalitySet));
  //  while(skip ? increaseSkip(equalitySet) : increase(equalitySet));

  return ret;
}


/*
Computes a monomial in the ideal. The monomial will have variables x_m...x_n-1
Recursive. The monomial is assumed to exist.
The ideal is assumed to be homogeneous
*/
Term computeTermInIdeal(PolynomialSet const &ideal, int m)
{
  assert(!ideal.empty());

  int n=ideal.begin()->numberOfVariablesInRing();

  //  fprintf(Stderr,"computTermInIdeal(%i,%i)\n",m,n);
  //  AsciiPrinter(Stderr).printPolynomialSet(ideal);


  if(m>=n)
    {
      Polynomial p=*ideal.begin();
      assert(!p.isZero());
      p.mark(LexicographicTermOrder());
      Term t=p.getMarked();
      t.m.exponent=IntegerVector(n);
      //fprintf(Stderr,"returning\n");
      return t;
    }
  PolynomialSet g=ideal;
  buchberger(&g,ReverseLexicographicTermOrder(m+1)); //is this the right order?
  //fprintf(Stderr,"GB:\n");
  //AsciiPrinter(Stderr).printPolynomialSet(g);

  PolynomialSet g2=g;
  for(PolynomialSet::iterator j=g.begin();j!=g.end();j++)
    j->saturate(m);
  //AsciiPrinter(Stderr).printPolynomialSet(g);
  buchberger(&g,ReverseLexicographicTermOrder(m)); //is this the right order?
  //AsciiPrinter(Stderr).printPolynomialSet(g);
  Term a=computeTermInIdeal(g,m+1);
  //fprintf(Stderr,"Returned term:");
  //AsciiPrinter(Stderr).printPolynomial(Polynomial(a));
  //fprintf(Stderr,"\n");
  for(int i=0;true;i++)
    {
      a.m.exponent[m]=i;
      Polynomial p=Polynomial(a);

      if(division(p,g2,LexicographicTermOrder()).isZero())return a;
    }
  assert(0);
  return a;
}

bool containsMonomialTraditional(PolynomialSet const &ideal)
{
  // Assuming the ideal is homogeneous

  if(ideal.empty())return false;

  for(PolynomialSet::const_iterator i=ideal.begin();i!=ideal.end();i++)
    if(i->isMonomial())return true;

  int n=ideal.begin()->numberOfVariablesInRing();

  PolynomialSet g=ideal;
  for(int i=0;i<n;i++)
    {
      buchberger(&g,ReverseLexicographicTermOrder(i));
     //  fprintf(Stderr,"before:");
      //  AsciiPrinter(Stderr).printPolynomialSet(g);
      for(PolynomialSet::iterator j=g.begin();j!=g.end();j++)
	j->saturate();
      //  fprintf(Stderr,"after:");
      //  AsciiPrinter(Stderr).printPolynomialSet(g);
    }
  buchberger(&g,ReverseLexicographicTermOrder(0));
  //  fprintf(Stderr,"after:");
  //  AsciiPrinter(Stderr).printPolynomialSet(g);

  return g.size()==1 && g.begin()->isMonomial();
}

bool containsMonomialDehomogenize(PolynomialSet const &ideal)
{
  PolynomialSet ideal2=ideal.multiDeHomogenization();
  PolynomialRing R3=ideal2.getRing().withVariablesAppended("H");
  return containsMonomialTraditional(ideal2.homogenization(R3));
}

bool containsMonomial(PolynomialSet const &ideal)
{
	log2 debug<<"containsMonomial() called on input ideal:\n"<<ideal;

	return containsMonomialDehomogenize(ideal);
}


PolynomialSet saturatedIdeal(PolynomialSet const &ideal)
{
  // Assuming the ideal is homogeneous

  assert(!ideal.empty());

  int n=ideal.begin()->numberOfVariablesInRing();

  PolynomialSet g=ideal;
  for(int i=0;i<n;i++)
    {
      log2 fprintf(Stderr,"Saturating with respect to variable %i.\n",i);
      buchberger(&g,ReverseLexicographicTermOrder(i));
      for(PolynomialSet::iterator j=g.begin();j!=g.end();j++)
	j->saturate();
    }
  buchberger(&g,ReverseLexicographicTermOrder(0));

  return g;
}
