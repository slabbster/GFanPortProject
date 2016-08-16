#include "dimension.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "division.h"
#include "field_rationals.h"
#include "buchberger.h"

class RealRootsApplication : public GFanApplication
{
public:
PolynomialSet sturmPolynomials(Polynomial  f1)
{
    Polynomial f2=f1.derivative();

    PolynomialRing theRing=f1.getRing();
    PolynomialSet result(theRing);
    result.push_back(f1);
    while(!f2.isZero())
      {
	result.push_back(f2);
	PolynomialSet g(theRing);
	Polynomial temp=f2;
	g.push_back(f2);
	g.markAndScale(LexicographicTermOrder());
	f2=(f1-f1)-division(f1,g,LexicographicTermOrder());
	f1=temp;
      }
    return result;
}



int numberOfSignChangesAtMinusInfinity(PolynomialSet const &l)
{
	int ret=0;
	int sign=0;
	for(PolynomialSet::const_iterator i=l.begin();i!=l.end();i++)
	{
		Polynomial p=*i;
		p.mark(LexicographicTermOrder());
		int newSign=p.getMarked().c.sign()*(1-2*(p.getMarked().m.exponent[0]&1));
		if(newSign && (newSign!=sign))
		{
			ret++;
			sign=newSign;
		}
	}
	return ret;
}


int numberOfSignChangesAtInfinity(PolynomialSet const &l)
{
	int ret=0;
	int sign=0;
	for(PolynomialSet::const_iterator i=l.begin();i!=l.end();i++)
	{
		Polynomial p=*i;
		p.mark(LexicographicTermOrder());
		int newSign=p.getMarked().c.sign();
		if(newSign && (newSign!=sign))
		{
			ret++;
			sign=newSign;
		}
	}
	return ret;
}


int numberOfSignChanges(PolynomialSet const &l, FieldElement const &a)
{
	int ret=0;
	int sign=0;
	for(PolynomialSet::const_iterator i=l.begin();i!=l.end();i++)
	{
		FieldElement v=i->evaluate(a);
		int newSign=v.sign();
		if(newSign && (newSign!=sign))
		{
			ret++;
			sign=newSign;
		}
	}
	return ret;
}

/**
 * Returns a negative number less than all roots of the polynomial whose Sturm sequence is given.
 */
FieldElement lowerBoundForRoots(PolynomialSet const &l)
{
	FieldElement ret=Q.zHomomorphism(-1);
	while(numberOfSignChangesAtMinusInfinity(l)!=numberOfSignChanges(l,ret))ret*=Q.zHomomorphism(2);
	return ret-Q.zHomomorphism(1);
}

FieldElement upperBoundForRoots(PolynomialSet const &l)
{
	FieldElement ret=Q.zHomomorphism(1);
	while(numberOfSignChangesAtInfinity(l)!=numberOfSignChanges(l,ret))ret*=Q.zHomomorphism(2);
	return ret+Q.zHomomorphism(1);
}

list<FieldElement> intervals(PolynomialSet const &l)
{
	list<FieldElement> ret;
	FieldElement lo(lowerBoundForRoots(l));

//	ret.push_back(a);
	while(1)
	{
		FieldElement hi=upperBoundForRoots(l);
	int cLo=numberOfSignChanges(l,lo);
	int cHi=numberOfSignChanges(l,hi);
	if(cLo<=cHi)break;
	{
		while(cHi!=cLo-1)
		{
			FieldElement med=(hi+lo)*Q.zHomomorphism(2).inverse();
			int cMed=numberOfSignChanges(l,med);
			if(cMed==cLo)
			{
				lo=med;
				cLo=cMed;
			}
			else
			{
				hi=med;
				cHi=cMed;
			}
		}
		ret.push_back(lo);
		ret.push_back(hi);
	}
	lo=hi;
	}
	return ret;
}

list<FieldElement> narrow(PolynomialSet const &l, list<FieldElement> intervals, FieldElement epsilon)
{
	list<FieldElement> ret;

	for(list<FieldElement>::const_iterator i=intervals.begin();i!=intervals.end();i++)
	{
		FieldElement lo=*i;i++;
		FieldElement hi=*i;

		while((hi-lo-epsilon).sign()>0)
		{
			FieldElement med=(hi+lo)*Q.zHomomorphism(2).inverse();
			if(numberOfSignChanges(l,med)==numberOfSignChanges(l,hi))
				hi=med;
			else
				lo=med;
		}
		ret.push_back(lo);
		ret.push_back(hi);
	}
	return ret;
}

//Returns a length that fist in between any two consequtive intervals
FieldElement smallestDistance(list<FieldElement> intervals)
{
	FieldElement ret=Q.zHomomorphism(1000);
	if(intervals.size()>=4)
	{
		list<FieldElement>::const_iterator i=intervals.begin();
		i++;
		for(;i!=intervals.end();i++)
		{
			FieldElement a=*i;i++;
			FieldElement diff=*i-a;
			if((diff-ret).sign()<0)ret=diff;
		}
	}
	return ret;
}
/*int numberOfRootsBetweenMinusInfinityAndHere(FieldElement )
{

}*/

bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "Not working yet. Given generators for a zero-dimensional ideal this program will compute all real points on the variety.\n";
  }
  RealRootsApplication() {
    registerOptions();
  }
  const char *name()
  {
    return "_realroots";
  }
  int main()
  {
    FileParser P(Stdin);

    PolynomialSet g=P.parsePolynomialSetWithRing();
    PolynomialRing r=g.getRing();
    int n=r.getNumberOfVariables();
//    Polynomial f1=P.parsePolynomialWithRing();
    WeightReverseLexicographicTermOrder T(IntegerVector::allOnes(n));
    buchberger(&g,T);
    debug<<g;
    int d=krullDimension(g);
    if(0!=d)
    {
		debug<<"Input ideal is not zero-dimensional.\n";
    	assert(0);
    }

    PolynomialRing r2(r.getField(),1);
    PolynomialSet projectionPolys(r2);

    for(int i=0;i<n;i++)
    {
    	LexicographicTermOrder T((i+1)%n);

    	buchberger(&g,T);

    	debug<<g;

    	list<int> l;
    	l.push_back(i);
    	PolynomialSet intersection=g.polynomialRingIntersection(r2,&l);
    	assert(intersection.size()==1);
    	projectionPolys.push_back(*intersection.begin());
    }
    debug<<projectionPolys;

    PolynomialSetList sturmPolys;
    for(PolynomialSet::const_iterator i=projectionPolys.begin();i!=projectionPolys.end();i++)
    {
			sturmPolys.push_back(sturmPolynomials(*i));
    	}
    debug.printPolynomialSetList(sturmPolys);


    FieldElement bound=Q.zHomomorphism(1);
    FieldElement distance=Q.zHomomorphism(1);
    for(PolynomialSetList::const_iterator i=sturmPolys.begin();i!=sturmPolys.end();i++)
    {
    	debug<<"Lower "<< lowerBoundForRoots(*i)<<"\n";
    	debug<<"Upper "<< upperBoundForRoots(*i)<<"\n";
    	debug<<"Sign changes "<< numberOfSignChangesAtMinusInfinity(*i)<<"\n";
    	FieldElement s=Q.zHomomorphism(-10000);
    	debug<<"Sign changes "<< numberOfSignChanges(*i,s)<<"\n";
    	FieldElement t=Q.zHomomorphism(100000);
    	debug<<"Sign changes "<< numberOfSignChanges(*i,t)<<"\n";
    	debug<<"Sign changes "<< numberOfSignChangesAtInfinity(*i)<<"\n";

    	list<FieldElement> l=intervals(*i);
    	FieldElement epsilon=(upperBoundForRoots(*i)-lowerBoundForRoots(*i))*Q.zHomomorphism(100).inverse();
    	l=narrow(*i,l,epsilon);


    	for(list<FieldElement>::const_iterator j=l.begin();j!=l.end();j++)
    		debug<<*j<<" "<< fieldElementToFloatingPoint(*j)<<"\n";

    	FieldElement b=upperBoundForRoots(*i)-lowerBoundForRoots(*i);
//    	if(bound<b)bound=b;
//    	FieldElement d=smallestDistance(l);
//    	if(d<distance)distance=d;
    }
/*
    FieldElement epsilon2=distance/bound; //<------------------fix this
	l=narrow(*i,l,epsilon);

	FieldElement delta=epsilon/bound;//<--------------------fix this

	PolynomialRing r3=r.withVariablesAppended("W");
	Polynomial f=-r3.ithVariable(n);
	FieldElement multiplier=Q.zHomomorphism(1);
	for(int i=0;i<n;i++)
	{
		f+=multiplier*r3.ithVariable(i);
		multiplier*=delta;
	}
	PolynomialSet g2=g.embeddedInto(r3);
	g2.push_back(f);
*/


/*    FieldRationalFunctions k(r.getField(),"t");
    PolynomialRing r2(k,n+1);
    Polynomial generic(r2)=-r2.ithVariable(n);
    for(int i=0;i<n;i++)
    	generic+=k.exponent(i)*r2.ithVariable(i);

    PolynomialSet g2(r2);
    g2.push_back(generic);

    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    {
    	Polynomial f(r2);
    	for(TermMap::const_iterator j=g.begin();j!=g.end();j++)
    	{
    		f+=Term(k.fromCoefficientField(),Monomial(r2,j->first.m.v));
    	}
    	g2.push_back(f);
    }
    pout << g2;
    LexicographicTermOrder T2;
    buchberger(g2,T2);
    Polynomial p(r2);
    for(PolynomialSet::const_iterator i=g2.begin();i!=g2.end();i++)
    {
    	if(i->numberOfVariablesInUseConsecutive()==1)p=*i;
    }

    PolynomialRing r3(k,1);

*/



//AsciiPrinter(Stdout).printPolynomialSet(result);

/*	if(evaluateOption.getValue())
      {
	IntegerVector v=P.parseIntegerVector();
	for(int i=0;i<v.size();i++)
	  {
	    FieldElement x=Q.zHomomorphism(v[i]);
	    AsciiPrinter(Stdout).printString("Evaluating in ");
	    AsciiPrinter(Stdout).printFieldElement(x);
	    AsciiPrinter(Stdout).printNewLine();
	    for(PolynomialSet::const_iterator j=result.begin();j!=result.end();j++)
	      {
		AsciiPrinter(Stdout).printFieldElement(j->evaluate(x));
		AsciiPrinter(Stdout).printNewLine();
	      }
	  }
      }
*/
    return 0;
  }
};

static RealRootsApplication theApplication;
