#include <sstream>
#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "lp.h"
#include "gfanapplication.h"
#include "polyhedralcone.h"
#include "tropical2.h"
#include "nbody.h"

#include "polymakefile.h"
#include "determinant.h"
#include "subspace.h"
#include "triangulation.h"
#include "groebnerengine.h"

class SmalesSixth2Application : public GFanApplication
{
  StringOption inputOption;
  SimpleOption singleOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program solves ........\n";
  }
  SmalesSixth2Application():
    inputOption("-i","Specify the name of the input file.","examples/nbody/5ACAC.out"),
    singleOption("-s","Chekc just a single initial system.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_smalessixth2";
  }

  void printIntList(list<int> const &v)
  {
    FILE *f=Stderr;
    fprintf(f,"{");
    for(list<int>::const_iterator i=v.begin();i!=v.end();i++)
      {
	if(i!=v.begin())fprintf(f," ");
	fprintf(f,"%i",*i);
      }
    fprintf(f,"}\n");
  }
  void printIntListList(list<list<int> > const &l)
  {
    for(list<list<int> >::const_iterator i=l.begin();i!=l.end();i++)
      printIntList(*i);
  }

  PolynomialSet saturationHomog(PolynomialSet const &g, int i)
  {
    WeightReverseLexicographicTermOrder T(-1*IntegerVector::standardVector(g.getRing().getNumberOfVariables(),i));
    PolynomialSet g2=GE_groebnerBasis(g,T,false);
    //    g2.saturate(i);
    g2.saturate();
    return g2;
  }

  PolynomialSet removePositiveFactors(PolynomialSet const &g)
  {
    PolynomialSet ret(g.getRing());

    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
      {
	Polynomial p=*i;
	if(p.numberOfTerms()==2)
	  {
	    FieldElement c1=p.terms.begin()->second;
	    FieldElement c2=p.terms.rbegin()->second;
	    IntegerVector e1=p.terms.begin()->first.exponent;
	    IntegerVector e2=p.terms.rbegin()->first.exponent;

	    if((c1+c2).isZero())
	      {
		int d=gcdOfVector(concatenation(e1,e2));
		if(d&1)
		  p=Polynomial(Term(c1,Monomial(g.getRing(),e1/d)))+Polynomial(Term(c2,Monomial(g.getRing(),e2/d)));
		else
		  p=Polynomial(Term(c1,Monomial(g.getRing(),2*e1/d)))+Polynomial(Term(c2,Monomial(g.getRing(),2*e2/d)));
	      }
	  }
	ret.push_back(p);
      }

    return ret;
  }

  PolynomialSet withVariablesEqual(PolynomialSet const &g, Polynomial const &b)
  {
    PolynomialSet ret(g.getRing());
    PolynomialSet B(g.getRing());
    B.push_back(b);
    LexicographicTermOrder T;
    B.markAndScale(T);

    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
      {
	ret.push_back(division(*i,B,T));
      }
    //    ret.push_back();
    return ret;
  }

  PolynomialSet binomialReduce(PolynomialSet g)
  {
    PolynomialSet binomials(g.getRing());
    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
      if(i->numberOfTerms()==2)binomials.push_back(*i);

    for(PolynomialSet::const_iterator i=binomials.begin();i!=binomials.end();i++)
      {
	g=withVariablesEqual(g,*i);
	g.saturate();
      }
    for(PolynomialSet::const_iterator i=binomials.begin();i!=binomials.end();i++)
      {
	g.push_back(*i);
      }
    return g;
  }
  void processSingle(PolynomialSet g, int N)
  {
	  	  {//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    PolynomialSet g2(g.getRing());
	    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
	      if(i->totalDegree()<20
		&& i->numberOfTerms()!=6 && i->numberOfTerms()<10
		 )g2.push_back(*i);
	    g=g2;
	    }


    AsciiPrinter P(Stderr);
	g.saturate();
	P<<g<<"\n";

	ReverseLexicographicTermOrder T;
	//LexicographicTermOrder T;
	//g.markAndScale(T);
	//g=initialForms(g,IntegerVector(15));

	//	g=removePositiveFactors(g);
	P<<g<<"\n";
	g.saturate();
	P<<g<<"\n";
	g=binomialReduce(g);
	P<<g<<"\n";
	//g=removePositiveFactors(g);
	P<<g<<"\n";
	g.saturate();
	P<<g<<"\n";
	g=binomialReduce(g);
	P<<g<<"\n";

	PolynomialRing r3(Q,N);
	PolynomialSet G=g.polynomialRingIntersection(r3);
	PolynomialSet G2(G.getRing());
	for(PolynomialSet::const_iterator i=G.begin();i!=G.end();i++)if(i->numberOfTerms()!=0)G2.push_back(*i);
	P<<"Intersection"<<G2<<"\n";

	if(G2.size()==0)
	{
	  PolynomialRing r(g.getRing().getField(),g.getRing().getNumberOfVariables()+1);



	  PolynomialSet g2=g.homogenization(r);

	  AsciiPrinter(Stdout)<<g2;



	  //    for(int i=0;i<g2.getRing().getNumberOfVariables();i++)
	  //	  for(int i=0;i<g.getRing().getNumberOfVariables();i++)
	    for(int i=g2.getRing().getNumberOfVariables()-1;i>=0;i--)
	    {
	      cerr << i << endl;
	      g2=saturationHomog(g2,i);


	      {	PolynomialSet G=g2.polynomialRingIntersection(r3);
	PolynomialSet G2(G.getRing());
	for(PolynomialSet::const_iterator i=G.begin();i!=G.end();i++)if(i->numberOfTerms()!=0)G2.push_back(*i);
	P<<"Intersection"<<G2<<"\n";
	      }

	    }
	  WeightReverseLexicographicTermOrder T(concatenation(IntegerVector::allOnes(N),40*IntegerVector::allOnes(1+N*(N-1)/2)));
	  g2=GE_groebnerBasis(g2,T,false);

	  PolynomialRing r3(Q,N);
	  PolynomialSet g3=g2.polynomialRingIntersection(r3);


	  P<<"Pseudo intersection:\n";
	  P<<g3;
	  //	if(g3.size())break;
	}
  }

  bool isDifficult(int I)
  {
    if(I>429)return true;
    if(I==41 || I==48 || I==180 || I==223 || I==235 || I==236 || I==264 || I==418 || I==419 || I==420 || I==421)return true;
    if(I==41 || I==48 || I==110 || I==143 || I==175 || I==180 || I==189 || I==192
       || I==193 || I==198 || I==228 || I==232 || I==235 || I==236 || I==254 || I==264
       || I==270 || I==298 || I==342 || I==399 || I==400 || I==418 || I==419 || I==420
       || I==421 || (I==423) || (I==424) || (I==425) || (I==426) || I==427 || I==428 || I==428 || I==429)return true;

    if(I==22 || I==23 || I==30 || I==36 || I==54 || I==61 || I==64 || I==108
       || I==149 || I==161 || I==162 || I==172 || I==173 || I==174 || I==177
       || I==178 || (I==187) || I==190 || I==191 || I==194 || I==195 || I==204
       || I==206 || I==207 || I==209 || I==211 || I==212 || I==213 || I==214
       || I==217 || I==219 || I==222 || I==224 || (I==225) || I==230 || (I==233)
       || I==234 || I==241 || (I==244) || I==246 || I==249 || I==255 || I==257
       || I==263 || I==271 || I==283 || I==338 || I==339 || I==345 || I==360
       || I==396 || I==397 || I==398 || I==401 || I==404 || I==409 || I==410 || I==415)return true;
    return false;
  }

  int main()
  {
    int N=5;
    LpSolver::printList(Stderr);
    lpSetSolver("cddgmp");


    if(singleOption.getValue())
      {
	PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();
	processSingle(g,N);
	return 0;
      }

    AsciiPrinter P(Stderr);

    //    cerr<< "Simplified"<<endl;

    //    P<<AlbouyChencinerEquationsSimplified(N);

    PolymakeFile inFile;
    inFile.open(inputOption.getValue());

    int n=inFile.readCardinalProperty("AMBIENT_DIM");
    int nRays=inFile.readCardinalProperty("N_RAYS");

    IntegerMatrix rays=inFile.readMatrixProperty("RAYS",nRays,n);

    /*    {// Make sure that representatives of rays are in the same subspace.
      IntegerVectorList rays2;
      IntegerMatrix linealitySpace=inFile.readMatrixProperty("LINEALITY_SPACE",nRays,n);
      Subspace l(linealitySpace.getRows(),linealitySpace.getWidth());
      for(int i=0;i<rays.getHeight();i++)
	rays2.push_back(l.canonicalizeVector(rays[i]));
      rays=rowsToIntegerMatrix(rays2,n);
    }*/

    vector<list<int> > cones=inFile.readMatrixIncidenceProperty("CONES_COMPRESSED");


    //    P<<AlbouyChencinerEquations(5);
    PolynomialSet ac=AlbouyChencinerEquations(N,true,false);
    PolynomialSet ac2=AlbouyChencinerEquations(N,true,true);
    //P<<DziobekEquations(ac.getRing());


    PolynomialSet all=ac;
    for(PolynomialSet::const_iterator i=ac2.begin();i!=ac2.end();i++)all.push_back(i->embeddedInto(ac.getRing()));

    P<<"ALL\n"<<all<<"\n";
    vector<int> skipped;
    int I=cones.size()-1;
    for(int i=0;i<cones.size();i++,I--)
      {
	IntegerVector v(n);
	bool hasNegativeSum=false;
	IntegerVectorList raySubset;
	for(list<int>::const_iterator j=cones[i].begin();j!=cones[i].end();j++)
	  {
	    v+=rays[*j];
	    if(dotLong(IntegerVector::allOnes(n),rays[*j])>=0)hasNegativeSum=true;
	    raySubset.push_back(rays[*j]);
	  }
	AsciiPrinter(Stderr)<<raySubset<<"\n";
	if(!hasNegativeSum)
	  {
	    cerr<<"SKIP"<<endl;
	    continue;
	  }
	//	r.push_front(v);

	cerr<<"CHECKING INDEX "<<i<<endl;
	if(isDifficult(I) || raySubset.size()==1)
	  skipped.push_back(i);
	else
	{
	  PolynomialSet g=initialForms(all,concatenation(IntegerVector(N),v));

	  P<<"---------------------------------------\n";
	  P<<"Index(backwards)";
	  cerr <<I;
	  P<<" weight "<<v<<"\n";

	  processSingle(g,N);
	}
      }
    IntegerVector skipp(skipped.size());
    for(int i=0;i<skipp.size();i++)skipp[i]=skipped[i];
    AsciiPrinter(Stdout) << skipp;

    return 0;
  }
};

static SmalesSixth2Application theApplication;

