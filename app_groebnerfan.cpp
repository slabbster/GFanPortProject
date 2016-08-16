#include <iostream>

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

class GroebnerFanApplication : public GFanApplication
{
  SimpleOption optionInputIsGroebnerBasis;
  SimpleOption optionSymmetry;
  SimpleOption optionDisableSymmetryTest;
  StringOption optionRestrictingFan;
  //StringOption optionRestrictingFanSymmetry;
  IntegerOption optionParameters;
  SimpleOption optionIgnoreCones;
  public:
  const char *helpText()
  {
    return "This is a program for computing the Groebner fan of a polynomial ideal as a polyhedral complex. It takes a generating set for the ideal as input. If the ideal is symmetric the symmetry option is useful and enumeration will be done up to symmetry. The program needs a starting Groebner basis to do its computations. If the -g option is not specified it will compute one using Buchberger's algorithm.\n";
  }
  GroebnerFanApplication():
    optionInputIsGroebnerBasis("-g",
			       "Tells the program that the input is already a Groebner basis (with the initial term of each polynomial being "
			       "the first ones listed). Use this option if it takes too much time to compute "
			       "the starting (standard degree lexicographic) Groebner basis and the input is already a Groebner basis.\n"),
    optionSymmetry("--symmetry",
		   "Tells the program to read in generators for a group of symmetries (subgroup of $S_n$) after having read in the ideal. The program checks that the ideal stays fixed when permuting the variables with respect to elements in the group. The program uses breadth first search to compute the set of reduced Groebner bases up to symmetry with respect to the specified subgroup.\n"),
    optionDisableSymmetryTest("--disableSymmetryTest","When using --symmetry this option will disable the check that the group read off from the input actually is a symmetry group with respect to the input ideal.\n"),
    optionRestrictingFan("--restrictingfan","Specify the name of a file containing a polyhedral fan in Polymake format. The computation of the Groebner fan will be restricted to this fan. If the --symmetry option is used then this restricting fan must be invariant under the symmetry and the orbits in the file must be with respect to the specified group of symmetries. The orbits of maximal cones of the file are then read in rather than the maximal cones.\n",0),
    //optionRestrictingFan("--restrictingfansymmetry","Use this option to indicate that the restricting fan also has the same symmetries.\n",0),
    optionParameters("--parameters","With this option you can specify how many variables to treat as parameters instead of variables. This makes it possible to do computations where the coefficient field is the field of rational functions in the parameters. This does not work well at the moment.",0),
    optionIgnoreCones("--nocones","Tells the program not to output the CONES and MAXIMAL_CONES sections, but still output CONES_COMPRESSED and MAXIMAL_CONES_COMPRESSED if --symmetry is used.")

  {
    registerOptions();
  }

  const char *name()
  {
    return "_groebnerfan";
  }

  int main()
  {
    LexicographicTermOrder myOrder;

    PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();
    if(optionParameters.getValue())g=makeVariablesParameters(makeVariablesParameters(g.getRing(),optionParameters.getValue()),g);
    int n=g.numberOfVariablesInRing();

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

      SymmetryGroup s(n);

    IntegerVectorList generators;
      if(optionSymmetry.getValue())
    	  generators=FileParser(Stdin).parseIntegerVectorList();
	    if(!optionDisableSymmetryTest.getValue())
	      {
		for(IntegerVectorList::iterator i=generators.begin();i!=generators.end();i++)
		  {
		    assert(areIdealsEqual(g,SymmetryGroup::permutePolynomialSet(g,*i)));
		  }
	      }

	  s.computeClosure(generators);
	  s.createTrie();
	  log3 s.print(Stderr);

		SymmetricTargetFanBuilder target(n,s);

	if(!optionRestrictingFan.getValue())
	{
		GroebnerFanTraverser traverser(g);
		symmetricTraverse(traverser,target,&s);
	}
	else
	{
	    //PolyhedralFan f1=PolyhedralFan::readFan(optionRestrictingFan.getValue());
	    //static PolyhedralFan readFan(string const &filename, bool onlyMaximal=true, IntegerVector *w=0, set<int> const *conesIndice=0, SymmetryGroup const *sym=0, bool readCompressedIfNotSym=false);
	    PolyhedralFan f1=PolyhedralFan::readFan(optionRestrictingFan.getValue(),true,0,0,/*optionSymmetry.getValue()?&s:0*/0,true);


	    for(PolyhedralFan::coneIterator i=f1.conesBegin();i!=f1.conesEnd();i++)
		{
	    	static int t;
	    	log2 cerr<<"Processing Cone "<<t++<<" which has dimension "<<i->dimension()<<endl;
	    	GroebnerFanTraverser traverser(groebnerBasisWithFullDimensionalIntersection(g,*i),*i);
			symmetricTraverse(traverser,target,&s);
		}
	}
	AsciiPrinter Q(Stdout);

  	target.getFanRef().printWithIndices(&Q,
				    (optionSymmetry.getValue()?FPF_group|FPF_conesCompressed:0)|
				    (optionIgnoreCones.getValue()?0:FPF_conesExpanded)|
				    FPF_maximalCones|FPF_cones,
				    &s);
    return 0;
  }
};

static GroebnerFanApplication theApplication;
