#include "dimension.h"
#include "buchberger.h"
#include "log.h"
#include "printer.h"

PolynomialSet radicalOfMonomialIdeal(PolynomialSet const &monomialGenerators)
{
  PolynomialRing theRing=monomialGenerators.getRing();
  PolynomialSet temp=monomialGenerators;

  temp.markAndScale(LexicographicTermOrder()); //just to make sure that some term is marked

  PolynomialSet ret(theRing);
  for(PolynomialSet::const_iterator i=temp.begin();i!=temp.end();i++)
    {
      IntegerVector e=i->getMarked().m.exponent;
      e=e.supportVector();
      ret.push_back(Polynomial(Term(i->getMarked().c,Monomial(theRing,e))));
    }
  return ret;
}

static bool increase(IntegerVector &v, int &numberOfOnes)
{
  int i=0;
  while(i<v.size() && v[i]==1)
    {
      v[i]=0;
      numberOfOnes--;
      i++;
    }
  if(i==v.size())return false;
  v[i]=1;
  numberOfOnes++;
  return true;
}


static void rek(IntegerVector &ones, IntegerVector &zeros, int nOnes, int nZeros, IntegerVectorList const &vectors, int &best)
{
  if(nOnes>best)best=nOnes;
  if(ones.size()-nZeros<best)return;
  if(nOnes+nZeros==ones.size())return;


  log3
    {
      fprintf(Stderr,"Ones:\n");
      AsciiPrinter(Stderr).printVector(ones);
      fprintf(Stderr,"Zeros:\n");
      AsciiPrinter(Stderr).printVector(zeros);

      AsciiPrinter(Stderr).printVectorList(vectors);
    }

  int index=0;
  for(int i=0;i<ones.size();i++,index++)
    if((!ones[i])&&(!zeros[i]))break;

  assert(index<ones.size());

  ones[index]=1;
  bool good=true;
  for(IntegerVectorList::const_iterator i=vectors.begin();i!=vectors.end();i++)
    if(i->divides(ones))
      {
	good=false;
	break;
      }

  if(good)
    {
      rek(ones,zeros,nOnes+1,nZeros,vectors,best);
    }
  ones[index]=0;

  IntegerVectorList vectorsSubset;

  for(IntegerVectorList::const_iterator i=vectors.begin();i!=vectors.end();i++)
    {
      if((*i)[index]==0)vectorsSubset.push_back(*i);
    }
  zeros[index]=1;
  rek(ones,zeros,nOnes,nZeros+1,vectorsSubset,best);
  //  rek(ones,zeros,nOnes,nZeros+1,vectors,best);
  zeros[index]=0;
}


int krullDimensionOfMonomialIdeal(PolynomialSet const &monomialGenerators)
{
  PolynomialSet temp=radicalOfMonomialIdeal(monomialGenerators);
  minimize(&temp);
  IntegerVectorList vectors;
  for(PolynomialSet::const_iterator i=temp.begin();i!=temp.end();i++)
    vectors.push_back(i->getMarked().m.exponent);

  int best=0;
  
  int n=monomialGenerators.getRing().getNumberOfVariables();
  IntegerVector zeros(n);
  IntegerVector ones(n);

  rek(ones,zeros,0,0,vectors,best);

  return best;
}

/*int krullDimensionOfMonomialIdeal(PolynomialSet const &monomialGenerators)
{
  PolynomialSet temp=radicalOfMonomialIdeal(monomialGenerators);
  minimize(&temp);
  IntegerVectorList vectors;
  for(PolynomialSet::const_iterator i=temp.begin();i!=temp.end();i++)
    vectors.push_back(i->getMarked().m.exponent);

  assert(!vectors.empty());
  int n=vectors.begin()->size();
  IntegerVector subset(n);
  int numberOfOnes=0;
  int dimension=0;

  do
    {
      //  AsciiPrinter(Stderr).printVector(subset);

      if(numberOfOnes>dimension)
	{
	  bool ok=true;
	  for(IntegerVectorList::const_iterator i=vectors.begin();i!=vectors.end();i++)
	    if(i->divides(subset))
	      {
		ok=false;
		break;
	      }
	  if(ok)
	    dimension=numberOfOnes;
	}
    }
  while(increase(subset,numberOfOnes));

  return dimension;
}
*/
int krullDimension(PolynomialSet const &groebnerBasis)
{
  return krullDimensionOfMonomialIdeal(groebnerBasis.markedTermIdeal());
}
