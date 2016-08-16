#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "termorder.h"
#include "gfanapplication.h"
#include "wallideal.h"
#include "polyhedralcone.h"

class WeightVectorApplication : public GFanApplication
{
  SimpleOption optionMultiple;
  SimpleOption optionEcho;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program computes a strictly positive interior point of the cone given by the marked Groebner basis on the input.\n";
  }
  WeightVectorApplication():
    optionMultiple("-m","Do the same thing for a list of Groebner bases. The output is the list of interior points.\n"),
    optionEcho("-e","Echo the input. The program will output the Groebner basis before outputting its weight vector.\n")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_weightvector";
  }

  void process(PolynomialSet const &g)
  {
    IntegerVectorList normals=wallInequalities(g);

    if(optionEcho.getValue())
      {
	AsciiPrinter(Stdout).printString("(");
	AsciiPrinter(Stdout).printPolynomialSet(g);
	AsciiPrinter(Stdout).printString(",");
      }
    PolyhedralCone C(normals,IntegerVectorList());
    AsciiPrinter(Stdout).printVector(intersection(PolyhedralCone::positiveOrthant(C.ambientDimension()),C).getRelativeInteriorPoint());

    if(optionEcho.getValue())
      AsciiPrinter(Stdout).printString(")");
  }
  int main()
  {
    LpSolver::printList(Stderr);
    lpSetSolver("cddgmp");

    FileParser P(Stdin);



    if(optionMultiple.getValue())
      {
	fprintf(Stdout,"{\n");
	PolynomialSetList l=P.parsePolynomialSetListWithRing();
	for(PolynomialSetList::const_iterator i=l.begin();i!=l.end();i++)
	  {
	    if(i!=l.begin())fprintf(Stdout,",\n");
	    process(*i);
	  }
	fprintf(Stdout,"}\n");
      }
    else
      {
	PolynomialSet p=P.parsePolynomialSetWithRing();
	process(p);
      }

    //interiorPoint(normals,Stdout,true);

    return 0;
  }
};

static WeightVectorApplication theApplication;
