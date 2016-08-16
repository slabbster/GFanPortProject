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

class SmalesSixthApplication : public GFanApplication
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
  SmalesSixthApplication():
    //    inputOption("-i","Specify the name of the input file.","examples/5body10sym.out"/*"examples/5bodyNEW.out2"*/)
    //    inputOption("-i","Specify the name of the input file.","examples/4bodyMoreEquations.out"/*"examples/5bodyNEW.out2"*/)
    //    inputOption("-i","Specify the name of the input file.","examples/nbody/5bodyACAC.out")
    //    inputOption("-i","Specify the name of the input file.","examples/nbody/5bodyACACDZ.out")
    //    inputOption("-i","Specify the name of the input file.","examples/nbody/5bodyACACDZDE.out"),
    inputOption("-i","Specify the name of the input file.","examples/nbody/5bodyACACDZDEMA.out"),
    singleOption("-s","Chekc just a single initial system.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_smalessixth";
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

  /*  PolynomialSet DziobekEquations(PolynomialRing const &r)
  {
    StringParser P("{r243*r353 - r233*r243*r353 - r233*r453 + r233*r243*r453 +r233*r353*r453 - r243*r353*r453,"
		   "   r253*r343 - r233*r253*r343 -  r233*r453 + r233*r253*r453 + r233*r343*r453 - r253*r343*r453,"
		   " -r253*r343 + r243*r253*r343 + r243*r353 - r243*r253*r353 - r243*r343*r353 + r253*r343*r353,"
		   " r143*r353 - r133*r143*r353 - r133*r453 + r133*r143*r453 + r133*r353*r453 - r143*r353*r453,"
		   " r153*r343 - r133*r153*r343 - r133*r453 + r133*r153*r453 + r133*r343*r453 - r153*r343*r453,"
		   " -r153*r343 + r143*r153*r343 +  r143*r353 - r143*r153*r353 - r143*r343*r353 + r153*r343*r353,"
		   " r143*r253 - r123*r143*r253 - r123*r453 + r123*r143*r453 +  r123*r253*r453 - r143*r253*r453,"
		   " r153*r243 - r123*r153*r243 -  r123*r453 + r123*r153*r453 + r123*r243*r453 - r153*r243*r453,"
		   " -r153*r243 + r143*r153*r243 + r143*r253 - r143*r153*r253 -  r143*r243*r253 + r153*r243*r253,"
		   " r133*r253 - r123*r133*r253 -  r123*r353 + r123*r133*r353 + r123*r253*r353 - r133*r253*r353,"
		   " r153*r233 - r123*r153*r233 - r123*r353 + r123*r153*r353 +  r123*r233*r353 - r153*r233*r353,"
		   " -r153*r233 + r133*r153*r233 +  r133*r253 - r133*r153*r253 - r133*r233*r253 + r153*r233*r253,"
		   " r133*r243 - r123*r133*r243 - r123*r343 + r123*r133*r343 +  r123*r243*r343 - r133*r243*r343,"
		   " r143*r233 - r123*r143*r233 -  r123*r343 + r123*r143*r343 + r123*r233*r343 - r143*r233*r343,"
		   " -r143*r233 + r133*r143*r233 + r133*r243 - r133*r143*r243 -  r133*r233*r243 + r143*r233*r243}");

    return P.parsePolynomialSet(r);
    }*/
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
    /*
      Om facstd() i Singular:
 Hvis to trinomiumsfaktorer afsluttes beregningen efter 6min med 97 baser
     Hvis foerste er binomium og andet trinomium afsluttes beregningen efter 5 min med 46 baser
     Hvis foerste er trionomium og andet binomium afsluttes efter 12 sek med 34 baser.
     Hvis begger er binomier afsluttes efter 1 sek med 29 baser.
*/

    //g.push_back(StringParser("r35^1-r25^1").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //    g.push_back(StringParser("r35^2+r35*r25+r25^2").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //g.push_back(StringParser("r45^1-r15^1").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //    g.push_back(StringParser("r45^2+r45*r15+r15^2").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //       	  g.push_back(StringParser("r34^1-r24^1").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  //    	  g.push_back(StringParser("r342+r24r34+r242").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //  	  g.push_back(StringParser("r35^1-r25^1").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		  //  g.push_back(StringParser("r352+r25r35+r252").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  /*g.push_back(StringParser("r452+r45r35+r352").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  g.push_back(StringParser("r342+r34r23+r232").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  	  g.push_back(StringParser("r352+r25r35+r252").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */


	  /*	  PolynomialSet additional=StringParser(
     "{r25^3*r34^3-r25^3*r34^3*r35^3-r24^3*r35^3+r24^3*r34^3*r35^3+r24^3*r25^3*r35^3-r24^3*r25^3*r34^3,"
     "-r13^3*r15^3*r34^3*r45^3+r13^3*r15^3*r34^3*r35^3*r45^3+r13^3*r14^3*r35^3*r45^3-r13^3*r14^3*r34^3*r35^3*r45^3-r13^3*r14^3*r15^3*r35^3*r45^3+r13^3*r14^3*r15^3*r34^3*r45^3,"
     "-r12^3*r15^3*r24^3*r45^3+r12^3*r15^3*r24^3*r25^3*r45^3+r12^3*r14^3*r25^3*r45^3-r12^3*r14^3*r24^3*r25^3*r45^3-r12^3*r14^3*r15^3*r25^3*r45^3+r12^3*r14^3*r15^3*r24^3*r45^3,"
     "-r12^3*r15^3*r23^3*r35^3+r12^3*r15^3*r23^3*r25^3*r35^3+r12^3*r13^3*r25^3*r35^3-r12^3*r13^3*r23^3*r25^3*r35^3-r12^3*r13^3*r15^3*r25^3*r35^3+r12^3*r13^3*r15^3*r23^3*r35^3,"
     "-r12^3*r14^3*r23^3*r34^3+r12^3*r14^3*r23^3*r24^3*r34^3+r12^3*r13^3*r24^3*r34^3-r12^3*r13^3*r23^3*r24^3*r34^3-r12^3*r13^3*r14^3*r24^3*r34^3+r12^3*r13^3*r14^3*r23^3*r34^3}").parsePolynomialSet(g.getRing());
	  for(PolynomialSet::const_iterator i=additional.begin();i!=additional.end();i++)
	    g.push_back(*i);
	  */
	  //"-r14^3+r14^3*r34^3+r13^3-r13^3*r34^3"
	  //	  g.push_back(StringParser("r14^3-r13^3").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  //	  g.push_back(StringParser("r343-1").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  //g.push_back(StringParser("r143-1").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  //	  g.push_back(StringParser("r34-1").parsePolynomial(g.getRing()));//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    /*{//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    PolynomialSet g2(g.getRing());
	    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
	      if(i->totalDegree()<20
		&& i->numberOfTerms()!=6 && i->numberOfTerms()<10
		 )g2.push_back(*i);
	    g=g2;
	    }*/


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
    IntegerVectorList r;

    for(int i=0;i<cones.size();i++)
      {
    	bool hasPos=false;
    	bool hasNeg=false;
    	IntegerVector v(n);
	for(list<int>::const_iterator j=cones[i].begin();j!=cones[i].end();j++)
	  {
		int sign=rays[*j].sum();
		if(sign<0)hasNeg=true;
		if(sign>0)hasPos=true;
	    v+=rays[*j];
	  }
	//r.push_front(v);
	if(hasPos ||((!hasPos)&&(!hasNeg)))
		{
		r.push_back(v);
		P.printInteger(r.size()-1);
		P<<cones[i];
		}
      }


    P<<"Weights to check";
    P.printVectorList(r,true);


    bool massEquations=false;//was true in spatial case
    bool dziobek=false;//was true in spatial case

    //    P<<AlbouyChencinerEquations(5);
    PolynomialSet ac=AlbouyChencinerEquations(N,true,false);
    PolynomialSet ac2=AlbouyChencinerEquations(N,true,true);
    //P<<DziobekEquations(ac.getRing());


    PolynomialSet all=ac;
    for(PolynomialSet::const_iterator i=ac2.begin();i!=ac2.end();i++)all.push_back(i->embeddedInto(ac.getRing()));

if(dziobek)
{
	PolynomialSet dziobek=DziobekEquations(ac.getRing(),N,true);

    for(PolynomialSet::const_iterator i=dziobek.begin();i!=dziobek.end();i++)all.push_back(*i);
}
    PolynomialSet d=nbodyDeterminants(ac.getRing(),N,true,5);
    for(PolynomialSet::const_iterator i=d.begin();i!=d.end();i++)all.push_back(*i);

    if(massEquations)all.push_back(massEquation(ac.getRing(),N,true));

    P<<"ALL\n"<<all<<"\n";

    int I=0;


    /*
backward index 51
{
x0^8*x1^9+x2^9*x3^8+3*x1^3*x2^6*x3^8+3*x1^6*x2^3*x3^8+x1^9*x3^8-x0*x2^9*x3^7-3*x0*x1^3*x2^6*x3^7-3*x0*x1^6*x2^3*x3^7-x0*x1^9*x3^7+x0^2*x2^9*x3^6+3*x0^2*x1^3*x2^6*x3^6+3*x0^2*x1^6*x2^3*x3^6+x0^2*x1^9*x3^6+2*x0^3*x2^9*x3^5-10*x0^3*x1^3*x2^6*x3^5-10*x0^3*x1^6*x2^3*x3^5+2*x0^3*x1^9*x3^5-2*x0^4*x2^9*x3^4+10*x0^4*x1^3*x2^6*x3^4+10*x0^4*x1^6*x2^3*x3^4-2*x0^4*x1^9*x3^4+2*x0^5*x2^9*x3^3-10*x0^5*x1^3*x2^6*x3^3-10*x0^5*x1^6*x2^3*x3^3+2*x0^5*x1^9*x3^3+x0^6*x2^9*x3^2+3*x0^6*x1^3*x2^6*x3^2+3*x0^6*x1^6*x2^3*x3^2+x0^6*x1^9*x3^2-x0^7*x2^9*x3-3*x0^7*x1^3*x2^6*x3-3*x0^7*x1^6*x2^3*x3-x0^7*x1^9*x3+x0^8*x2^9+3*x0^8*x1^3*x2^6+3*x0^8*x1^6*x2^3}


backward index 57
{
x2^4-14544/1015*x3^4-324/145*x2*x3^3-916/203*x2^2*x3^2-20/29*x2^3*x3+180/29*x1*x3^3+2964/1015*x1*x2*x3^2+64/145*x1*x2^2*x3-2743/203*x1^2*x3^2,
x1*x2^3-11520/4147*x3^4+8181/4147*x2*x3^3-2437/4147*x2^2*x3^2+2525/4147*x2^3*x3-305/143*x2^4-22725/4147*x1*x3^3+2960/319*x1*x2*x3^2-1616/4147*x1*x2^2*x3+3065/319*x1^2*x3^2,
x1^2*x2+144/29*x3^3-99/29*x2*x3^2-25/29*x2^2*x3+50/29*x2^3-189/29*x1*x3^2+41/29*x1*x2*x3-61/29*x1*x2^2-25/29*x1^2*x3,
x1^3+1296/145*x3^3-1152/145*x2*x3^2-45/29*x2^2*x3+61/29*x2^3-288/29*x1*x3^2+369/145*x1*x2*x3-288/145*x1*x2^2-45/29*x1^2*x3}
    */

	IntegerVector skipVector=FileParser(Stdin).parseIntegerVector();

	for(IntegerVectorList::const_iterator i=r.begin();i!=r.end();i++,I++)
      {
		P<<"Processing ";
		P.printInteger(I);
		P<<"\n";
		bool skip=false;
		for(int j=0;j<skipVector.size();j++)if(skipVector[j]==I)skip=true;
		if(skip)
		{
			P<<"Skipping \n";
			continue;
		}

	/*(0,1,4,5,6,8,9,10,15,16,17,18,19,24,25,29,30,39,40,41,43,51,52,54,55,58,59,63,65,70,73,75,92,96,99,104,107)*/
//	if(I!=98)continue;
	//if(I!=21) continue;
	//if(I!=47) continue;
	/*		if((I==10  || I==13 || I==18 || I==21 || I==25 || I==42 || I==44 || I==47 || I==52
	    || I==54 || I==58 || I==59 || I==62 || I==63 || I==65 || I==66 || I==74 || I==76
	    || I==77 ||  I==78 || I==87 || I==88 || I==92 || I==93 || I==98 ||I==99 || I==100 || I==101
	    || I==102 || I==107 || I==108 || I==109 || I==111 || I==112 || I==113 || I==116 || I==117)) continue;
	*/
// || I==13 || I==21 || I==56 || I==57 || I==60 || I==61 || I==63 || I==64 || I==72 || I==74 || I==103))continue;

	//	if(I<104) continue;
	//if(I==50 || I==54 || I==59 || I==61 || I==87 || I==92 || I==93 || I==94 || I==95 || I==100 || I==102 || I==104)continue;
	//	if(I<58)continue;
	//if(!(/*I==89||*/I==93))continue;

	/*	if((I==10 || I==13 || I==21 || I==56 || I==57 || I==60 || I==61 || I==63 || I==64 || I==72 || I==74 || I==103))continue;

	if(I==50 || I==52 || I==89 || I==93 || I==95 || I==96 || I==97 || I==102 || I==104 || I==106 || I==107 || I==111)continue;//ACACDZDE
	*/

	//if(I==42 || I==44 || I==56)continue;//ACACDZ


	//ACAC
	//	if(I==41 || I==48 || I==180 || I==223 || I==235 || I==236 || I==264 || I==418 || I==419 || I==420 || I==421)continue;
	/*	if(I==41 || I==48 || I==110 || I==143 || I==175 || I==180 || I==189 || I==192
	   || I==193 || I==198 || I==228 || I==232 || I==235 || I==236 || I==254 || I==264
	   || I==270 || I==298 || I==342 || I==399 || I==400 || I==418 || I==419 || I==420
	   || I==421 || (I==423) || (I==424) || (I==425) || (I==426) || I==427 || I==428 || I==428 || I==429)continue;
	*/
	/*	if(I==22 || I==23 || I==30 || I==36 || I==54 || I==61 || I==64 || I==108
	   || I==149 || I==161 || I==162 || I==172 || I==173 || I==174 || I==177
	   || I==178 || (I==187) || I==190 || I==191 || I==194 || I==195 || I==204
	   || I==206 || I==207 || I==209 || I==211 || I==212 || I==213 || I==214
	   || I==217 || I==219 || I==222 || I==224 || (I==225) || I==230 || (I==233)
	   || I==234 || I==241 || (I==244) || I==246 || I==249 || I==255 || I==257
	   || I==263 || I==271 || I==283 || I==338 || I==339 || I==345 || I==360
	   || I==396 || I==397 || I==398 || I==401 || I==404 || I==409 || I==410 || I==415)continue;
	*/
	PolynomialSet g=initialForms(all,concatenation(IntegerVector(N),*i));

	/*	if(N==4)
	  g.push_back(StringParser("1-r12r13r14r23r24r34").parsePolynomial(g.getRing()));
	if(N==5)
	  g.push_back(StringParser("1-r12r13r14r15r23r24r25r34r35r45").parsePolynomial(g.getRing()));
	*/

	P<<"---------------------------------------\n";
//	P<<"Index(backwards)";
	P<<"Index";
	cerr <<I;
	P<<" weight "<<*i<<"\n";

	processSingle(g,N);

      }

	/*
	WeightReverseLexicographicTermOrder T(concatenation(IntegerVector::allOnes(N),2*IntegerVector::allOnes(N*(N-1)/2)));

	PolynomialSet gSubset(g.getRing());
	for(PolynomialSet::reverse_iterator j=g.rbegin();j!=g.rend();j++)
	  {
	    gSubset.push_back(*j);
	    P<<gSubset;
	PolynomialSet g2=GE_groebnerBasis(gSubset,T,false);
	PolynomialRing r3(Q,N);
	PolynomialSet g3=g2.polynomialRingIntersection(r3);
	*/


    return 0;
  }
};

static SmalesSixthApplication theApplication;

