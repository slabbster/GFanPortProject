#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "minkowskisum.h"
#include "newtonpolytope.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "tropical.h"
#include "division.h"
#include "bergman.h"
#include "tropical2.h"
#include "dimension.h"
#include "timer.h"
#include "log.h"
#include "linalg.h"
#include "tropicaltraverse.h"
#include "traverser_tropical.h"
#include "symmetrictraversal.h"
#include "traverser_stableintersection.h"

class TropicalTraverseApplication : public GFanApplication
{
  //  SimpleOption optionPerformanceTimer;
  SimpleOption optionSymmetry;
  //  SimpleOption optionNoIncidencePrinting;
  SimpleOption optionTorusSymmetry;
  SimpleOption optionIgnoreCones;
  SimpleOption optionDisableSymmetryTest;
  SimpleOption optionStableIntersection;
  //  SimpleOption optionInitialIdealPrinting;
  //  SimpleOption optionLargeDimensionalPrinting;
public:
  bool includeInDefaultInstallation()
  {
    return true;
  }
  const char *helpText()
  {
    return "This program computes a polyhedral fan representation of the tropical variety of a homogeneous prime ideal $I$. Let $d$ be the Krull dimension of $I$ and let $\\omega$ be a relative interior point of $d$-dimensional Groebner cone contained in the tropical variety. The input for this program is a pair of marked reduced Groebner bases with respect to the term order represented by $\\omega$, tie-broken in some way. The first one is for the initial ideal $in_\\omega(I)$ the second one for $I$ itself. The pair is the starting point for a traversal of the $d$-dimensional Groebner cones contained in the tropical variety. If the ideal is not prime but with the tropical variety still being pure $d$-dimensional the program will only compute a codimension $1$ connected component of the tropical variety.\n";
  }
  TropicalTraverseApplication():
    //    optionPerformanceTimer("-T","Enable performance timer"),
    optionSymmetry("--symmetry","Do computations up to symmetry and group the output accordingly. If this option is used the program will read in a list of generators for a symmetry group after the pair of Groebner bases have been read. Two advantages of using this option is that the output is nicely grouped and that the computation can be done faster."),
    //    optionNoIncidencePrinting("--noincidence","Disables printing of incidence information. This will avoid extracting the full fan from its orbits and avoid computation of the lower dimensional faces. The output is a list of orbits of maximal cones and a list of orbits of codimension 1 cones."),
    // optionInitialIdealPrinting("--initialideals","Print the initial ideals during incidence printing ONLY without --noincidence."),
    //    optionLargeDimensionalPrinting("--largedimensional","Print the full dimensional cones and the codimension 1 cones computed in the traversal.")
    optionDisableSymmetryTest("--disableSymmetryTest","When using --symmetry this option will disable the check that the group read off from the input actually is a symmetry group with respect to the input ideal.\n"),
    optionTorusSymmetry("--symsigns","Specify for each generator of the symmetry group an element of ${-1,+1}^n$ which by its multiplication on the variables together with the permutation will keep the ideal fixed. The vectors are given as the rows of a matrix."),
    optionIgnoreCones("--nocones","Tells the program not to output the CONES and MAXIMAL_CONES sections, but still output CONES_ORBITS and MAXIMAL_CONES_ORBITS if --symmetry is used."),
    optionStableIntersection("--stable","Traverse the stable intersection or, equivalently, pretend that the coefficients are genereric.")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_tropicaltraverse";
  }
  int main()
  {
    PolynomialSetList tropical;

    //    lpSetSolver("cddgmp");
    FileParser P(Stdin);

    AsciiPrinter p(Stdout);
    PolynomialSet coneGroebnerBasis=P.parsePolynomialSetWithRing();
    PolynomialSet idealGroebnerBasis=P.parsePolynomialSet(coneGroebnerBasis.getRing());
    coneGroebnerBasis.changeNumberOfVariables(idealGroebnerBasis.numberOfVariablesInRing());

    SymmetryGroup s(idealGroebnerBasis.numberOfVariablesInRing());
    if(optionSymmetry.getValue())
      {
	IntegerVectorList generators=P.parseIntegerVectorList();
	FieldMatrix *torusAction=0;
	FieldMatrix m(idealGroebnerBasis.getRing().getField(),0,0);
	if(optionTorusSymmetry.getValue())
	  {
	    m=integerMatrixToFieldMatrix(rowsToIntegerMatrix(P.parseIntegerVectorList()),idealGroebnerBasis.getRing().getField());
	    torusAction=&m;
	  }
    if(!optionDisableSymmetryTest.getValue())
    	assertSymmetriesMatch(generators,idealGroebnerBasis,torusAction);
	s.computeClosure(generators);
	  s.createTrie();
	log2 s.print(Stderr);
	log2 fprintf(Stderr,"\n");
      }
    if(optionTorusSymmetry.getValue() && !optionSymmetry.getValue())
      {
	fprintf(Stderr,"Option --torus can only be used in combination with option --symmetry.\n");
      }

    int n=idealGroebnerBasis.numberOfVariablesInRing();
    PolyhedralCone hspace=homogeneitySpace(idealGroebnerBasis);
    IntegerVectorList hv=hspace.dualCone().getEquations();
    int h=hv.size();
    int d=krullDimension(idealGroebnerBasis);

    log1 fprintf(Stderr,"Ambient dimension: %i\n",n);
    log1 fprintf(Stderr,"Dimension of homogeneity space: %i\n",h);
    log1 fprintf(Stderr,"Dimension of tropical variety: %i\n",d);
    //    fprintf(Stdout,"Homogeneity space:\n");


    {
      IntegerVectorList ridges=StringParser("{(7,-8,-8,2,7,-8,7,7,2,-8,1,1,1,-4,1),"
					    //"(9,-11,-11,-1,14,-3,7,7,2,-13,-6,4,4,-1,-1),"
					    //"(7,-8,-8,2,7,4,4,4,-1,-11,-11,4,4,-1,4),"
					    //"(11,-4,-14,-4,11,2,2,7,2,-13,-13,2,7,2,2),"
					    /////"(9,-11,-11,-1,14,-3,7,7,2,-13,-6,4,4,-1,-1),"
					    //"(2,-8,-8,2,12,-1,4,4,-1,-6,-1,4,4,-1,-6),"
					    //"(6,-4,-14,-4,16,-3,2,7,2,-8,-3,2,7,2,-8),"
					    //"(9,-11,-11,-1,14,-3,7,7,2,-13,-6,4,4,-1,-1),"
					    "(11,-4,-14,-4,11,-10,5,10,5,-10,-1,-1,4,-1,-1)}").parseIntegerVectorList();

      IntegerVectorList rays=StringParser("{(3,0,0,-6,3,0,-1,-1,2,0,-3,1,1,4,-3),"
					  //"(5,0,0,0,-5,9,-3,-3,-3,0,-14,3,3,3,5),"
					  //"(1,0,0,-2,1,-1,0,0,1,0,0,0,0,1,-1),"
					  //"(0,-2,0,2,0,0,1,0,-1,0,0,1,0,-1,0),"
					  /////"(-8,0,0,0,8,3,-1,-1,-1,0,5,1,1,1,-8),"
					  //"(2,0,0,-2,0,-1,0,0,1,0,-1,0,0,1,0),"
					  //"(0,-2,0,2,0,0,1,0,-1,0,0,1,0,-1,0),"
					  //"(3,0,0,0,-3,-12,4,4,4,0,9,-4,-4,-4,3),"
					  "(-8,-2,0,2,8,3,0,-1,-2,0,5,2,1,0,-8)}").parseIntegerVectorList();


      //      coneChangeDebugger(coneGroebnerBasis,idealGroebnerBasis,ridges,rays);

      //return 0;
    }

    PolyhedralFan p1(0);

    if(1)
    {
    	if(optionStableIntersection.getValue())
    	{
    		StableIntersectionTraverser traverser(coneGroebnerBasis,idealGroebnerBasis);
    		SymmetricTargetFanBuilder target(n,s);
    		symmetricTraverse(traverser,target,&s);
    		p1=target.getFanRef();
    	}
    	else
    	{
    		TropicalTraverser traverser(coneGroebnerBasis,idealGroebnerBasis);
    		SymmetricTargetFanBuilder target(n,s);
    		symmetricTraverse(traverser,target,&s);
    		p1=target.getFanRef();
    	}
    }
    	else
    if(1)
      {
	BergmanFan f=bergman(coneGroebnerBasis,idealGroebnerBasis,&s);
	f.computeMultiplicities();
	/*	log1 fprintf(Stderr,"Is simplicial: %s\n",f.isSimplicial()?"true":"false");*/
	log1 fprintf(Stderr,"Order of input symmetry group: %i\n",s.elements.size());
	log1 fprintf(Stderr,"Number of maximal cones: %i\n",f.numberOfMaximalCones());
	log1 fprintf(Stderr,"Modulo the homogeneity space:\n");
	log1 AsciiPrinter(Stderr).printVectorList(hv);
	log1 fprintf(Stderr,"Converting representation to a polyhedral complex modulo symmetry...\n");
	p1=f.toPolyhedralFan();
	log1 fprintf(Stderr,"Done converting representation to a polyhedral complex modulo symmetry.\n");
      }
    else
      p1=tropicalTraverse(coneGroebnerBasis,idealGroebnerBasis,&s);


    /*    if(optionNoIncidencePrinting.getValue())
      f.print(p);
      else*/
      {

	p1.removeAllLowerDimensional();

	AsciiPrinter Q(Stdout);

	p1.printWithIndices(&Q,
			    FPF_multiplicities|
			    (optionSymmetry.getValue()?FPF_group|FPF_conesCompressed:0)|
			    (optionIgnoreCones.getValue()?0:FPF_conesExpanded)|
			    FPF_maximalCones|FPF_cones,
			    &s);
	//	p1.printWithIndices(&Q,true,&s,optionSymmetry.getValue(),optionIgnoreCones.getValue());
      }
    /*    if(optionInitialIdealPrinting.getValue())
      {
	//	    AsciiPrinter p(Stderr);
	PolyhedralFan p1=f.toPolyhedralFan();
	int n=p1.getAmbientDimension();
	IncidenceList a=p1.getIncidenceList(&s);
	fprintf(Stderr,"Computing rays...\n");
	IntegerVectorList rays=p1.getRaysInPrintingOrder(&s);
	fprintf(Stderr,"Done computing rays.\n");
	p.printString("Rays:\n");
	p.printVectorList(rays);
	vector<IntegerVector> rays2(rays.size());
	int K=0;
	for(IntegerVectorList::const_iterator k=rays.begin();k!=rays.end();k++)
	  rays2[K++]=*k;

	for(IncidenceList::const_iterator j=a.begin();j!=a.end();j++)
	  {
	    p.printInteger(j->first);
	    for(IntegerVectorList::const_iterator i=j->second.begin();i!=j->second.end();i++)
	      {
		p.printVector(*i);
		IntegerVector v(n);
		for(int t=0;t<i->size();t++)
		  {
		    v+=rays2[(*i)[t]];
		  }
		p.printVector(v);
		PolynomialSet g2=idealGroebnerBasis;
		buchberger(&g2,WeightReverseLexicographicTermOrder(v));
		g2=initialFormsAssumeMarked(g2,v);
		g2=saturatedIdeal(g2);
		p.printPolynomialSet(g2);

	      }
	  }
	  }*/

    //    if(optionPerformanceTimer.getValue())Timer::printList();
    return 0;
  }
};

static TropicalTraverseApplication theApplication;
