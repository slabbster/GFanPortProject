/*
  TO DO:
  Remove Timer
  Remove Subspace option and maybe add an option for restricting to an arbitrary polyhedral cone instead.
  Remove Minkowski code
 */

#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "breadthfirstsearch.h"
#include "termorder.h"
#include "ep_standard.h"
#include "ep_xfig.h"
#include "field_rationalfunctions2.h"
#include "gfanapplication.h"
#include "timer.h"
#include "log.h"

#include "traverser_groebnerfan.h"

static Timer globalTimer("Global-timer",1);

class GCats : public GFanApplication
{
  const char *name1;
  SimpleOption optionInputIsGroebnerBasis;
  SimpleOption optionSymmetry;
  SimpleOption optionEchoSymmetry;
  SimpleOption optionSubspace;
  //  SimpleOption optionPerformanceTimer;
  SimpleOption optionDisableSymmetryTest;
  IntegerOption optionParameters;
public:
  bool includeInDefaultInstallation()
  {
    return name1[0]!=0;//Make sure that "gfan" itself does not go to the documentation section of the manual.
  }
  const char *helpText()
  {
#define HELP "This is a program for computing all reduced Groebner bases of a polynomial ideal. It takes the ring and a generating set for the ideal as input. By default the enumeration is done by an almost memoryless reverse search. If the ideal is symmetric the symmetry option is useful and enumeration will be done up to symmetry using a breadth first search. The program needs a starting Groebner basis to do its computations. If the -g option is not specified it will compute one using Buchberger's algorithm.\n"
//    if(name1[0])
//      return "This is the Gfan program for computing Groebner fans and tropical vaieties"HELP;
    return HELP;
  }
  GCats(const char *name2):
    name1(name2),
    //    optionPerformanceTimer("-T",
    // 			   "Enable performance timer.\n"),
    optionInputIsGroebnerBasis("-g",
			       "Tells the program that the input is already a Groebner basis (with the initial term of each polynomial being "
			       "the first ones listed). Use this option if it takes too much time to compute "
			       "the starting (standard degree lexicographic) Groebner basis and the input is already a Groebner basis.\n"),
    optionSymmetry("--symmetry",
		   "Tells the program to read in generators for a group of symmetries (subgroup of $S_n$) after having read in the ideal. The program checks that the ideal stays fixed when permuting the variables with respect to elements in the group. The program uses breadth first search to compute the set of reduced Groebner bases up to symmetry with respect to the specified subgroup.\n"),
    optionEchoSymmetry("-e","Echo. Output the generators for the symmetry group.\n"),
    optionSubspace("--subspace",
		   "Only do breadth first search on cones with their interior intersecting a specified subspace. The subspace is given by a list of hyperplane normals at the end of the input. The intersection of the hyperplanes is the subspace being specified. Note that the set of Groebner cones intersecting the subspace could be disconnected and that only one connected component is computed. Works only together with --symmetry.\n"),
    optionDisableSymmetryTest("--disableSymmetryTest","When using --symmetry this option will disable the check that the group read off from the input actually is a symmetry group with respect to the input ideal.\n"),
    optionParameters("--parameters","With this option you can specify how many variables to treat as parameters instead of variables. This makes it possible to do computations where the coefficient field is the field of rational functions in the parameters. This does not work well at the moment.",0)
  {
    registerOptions();
    optionSubspace.hide();//Is not supported anymore
  }

  const char *name()
  {
    return name1;
    //    return "";
  }

  int main()
  {
    if(name1[0]==0)
      {
        debug<<"This is the Gfan program for computing Groebner fans and tropical varieties.\n"
        "Use the command \"gfan list\" to view all subcommands.\n"
        "The command \"gfan\" is deprecate for computing all Groebner bases of an ideal.\n"
        "Please use subcommand \"gfan _bases\" instead. Awaiting input. <Ctrl>-D to end.\n";
      }

    LexicographicTermOrder myOrder;

    PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();
    if(optionParameters.getValue())g=makeVariablesParameters(makeVariablesParameters(g.getRing(),optionParameters.getValue()),g);

    log3 AsciiPrinter(Stderr).printPolynomialSet(g);
    AsciiPrinter(Stdout).printPolynomialRing(g.getRing());

    printf("\n");

    EnumerationFilePrinter *ep;

    {
      ep=new StandardEnumerationPrinter();
    }


    bool outputLatex=true;

    Printer *P;
    LatexPrinter *Q;
    FILE *latexFile;

    globalTimer.on();
    {
      if(optionInputIsGroebnerBasis.getValue())
	{
	  log1 fprintf(Stderr,"Minimizing and autoreducing input...\n");
	  minimize(&g);
	  autoReduce(&g, LexicographicTermOrder());
	}
      else
	{
	  log1 fprintf(Stderr,"Computing Groebner Basis...\n");
	  buchberger(&g,StandardGradedLexicographicTermOrder());
	  log2 AsciiPrinter(Stderr).printPolynomialSet(g);
	}
      log1 fprintf(Stderr,"A reduced Groebner basis has been computed\n");
    }

    SymmetryGroup s(g.numberOfVariablesInRing());

    IntegerVectorList generators;
    {
      EnumerationAlgorithm *rs;
      if(optionSymmetry.getValue())
	{
	  generators=FileParser(Stdin).parseIntegerVectorList();
	    if(!optionDisableSymmetryTest.getValue())
	      {
		for(IntegerVectorList::iterator i=generators.begin();i!=generators.end();i++)
		  {
		    //      fprintf(Stderr,"testing\n");
		    assert(areIdealsEqual(g,SymmetryGroup::permutePolynomialSet(g,*i)));
		  }
	      }

	  s.computeClosure(generators);
	  log1 s.print(Stderr);

	  if(0)
	  {//using old breadth first traversal
		  BreadthFirstSearch *bs=new BreadthFirstSearch(s,/*minkowski*/0);
		  /*	  if(optionSubspace.getValue())
			bs->setSubspace(FileParser(Stdin).parseIntegerVectorList());
		   */
		  rs=bs;

		  if(optionEchoSymmetry.getValue())AsciiPrinter(Stdout).printVectorList(generators);
		  ep->open(Stdout);

		  rs->setEnumerationTarget(ep);
		  rs->enumerate(g);
		  delete rs;
	  }
	  else
	  {//using new traversal
		GroebnerFanTraverser traverser(g);
		TargetGlue target(*ep);
		if(optionEchoSymmetry.getValue())AsciiPrinter(Stdout).printVectorList(generators);
		ep->open(Stdout);
		symmetricTraverse(traverser,target,&s);
	  }

	}
      else
      {
    	  rs=new ReverseSearch(myOrder);
      ep->open(Stdout);
      rs->setEnumerationTarget(ep);
      rs->enumerate(g);
      delete rs;
      }
    }

    ep->close();
    delete ep;

    printf("\n");

    globalTimer.off();
	//    if(optionPerformanceTimer.getValue())Timer::printList();

    return 0;
  }
};

static GCats theApplication("_bases");
static GCats theApplication2("");

