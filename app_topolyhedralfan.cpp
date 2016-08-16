#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "wallideal.h"
#include "lp.h"
#include "polyhedralcone.h"
#include "gfanapplication.h"
#include "polyhedralfan.h"
#include "halfopencone.h"
#include "symmetry.h"

class ToPolyhedralFanApplication : public GFanApplication
{
  SimpleOption optionRestrict;
  //  SimpleOption optionPair;
  SimpleOption optionSymmetry;
public:
  const char *helpText()
  {
    return "This program takes a list of reduced Groebner bases and produces the fan of all faces of these. In this way by giving the complete list of reduced Groebner bases, the Groebner fan can be computed as a polyhedral complex. The option --restrict lets the user choose between computing the Groebner fan or the restricted Groebner fan.\n";
  }
  ToPolyhedralFanApplication():
    optionSymmetry("--symmetry",
		   "Tell the program to read in generators for a group of symmetries (subgroup of $S_n$) after having read in the ring. The output is grouped according to these symmetries. Only one representative for each orbit is needed on the input.\n"),
    optionRestrict("--restrict","Add an inequality for each coordinate, so that the the cones are restricted to the non-negative orthant.")
    //    optionPair("--pair","The Groebner cone is given by a pair of compatible Groebner bases. The first basis is for the initial ideal and the second for the ideal itself. See the tropical section of the manual.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_topolyhedralfan";
  }

  int main()
  {
    PolynomialRing r=FileParser(Stdin).parsePolynomialRing();

    int n=r.getNumberOfVariables();

    SymmetryGroup s(n);
    if(optionSymmetry.getValue())
      s.computeClosure(FileParser(Stdin).parseIntegerVectorList());

    PolyhedralFan F(n);

    PolynomialSetList l=FileParser(Stdin).parsePolynomialSetList(r);

    for(PolynomialSetList::const_iterator i=l.begin();i!=l.end();i++)
      {
	PolynomialSet g=*i;
	PolynomialSet m=g.markedTermIdeal();

	IntegerVectorList equalities=wallInequalities(m);
	IntegerVectorList normals=algebraicTest(wallInequalities(g),g);
	if(optionRestrict.getValue())
	  {
	    for(int i=0;i<n;i++)
	      normals.push_back(IntegerVector::standardVector(n,i));
	  }

	PolyhedralCone c(normals,equalities,n);
	c.canonicalize();
	//      AsciiPrinter(Stdout).printPolyhedralCone(c);

	/*	IntegerVectorList empty;
	HalfOpenCone C(n,c.getEquations(),c.getHalfSpaces(),empty,true);
	PolyhedralFan F=faceComplexOfCone(C);
	*/
	F.insert(c);
      }

    AsciiPrinter P(Stdout);
    F.printWithIndices(&P,
		       FPF_default|
		       (optionSymmetry.getValue()?FPF_group|FPF_conesCompressed:0),
		       &s);
    //    F.printWithIndices(&P,false,&s,optionSymmetry.getValue());
    return 0;
  }
};

static ToPolyhedralFanApplication theApplication;
