#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "termorder.h"
#include "ep_standard.h"
#include "ep_xfig.h"
#include "gfanapplication.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "tropical.h"
#include "tropical2.h"
#include "symmetry.h"
#include "halfopencone.h"
#include "symmetrictraversal.h"
#include "traverser_groebnerfan.h"
#include "tropical_weildivisor.h"
#include "log.h"

class SymmetricTargetTropicalBasisTester : public SymmetricTarget
{
public:
	PolynomialSet g;
	SymmetricTargetTropicalBasisTester(PolynomialSet const &g_):
		g(g_)
		{

		}
	bool process(ConeTraverser &traverser)
	{
		IntegerVector w=traverser.refToPolyhedralCone().getRelativeInteriorPoint();
		log2 AsciiPrinter(Stderr) << "Testing weight vector:\n"<<w<<"\n";
		WeightReverseLexicographicTermOrder T(w);
		buchberger(&g,T);
		PolynomialSet temp=initialForms(g,w);
		if(containsMonomial(temp))
		{
			AsciiPrinter(Stdout)<<"The following vector is in intersection, but initial ideal contains a monomial:\n"<<w;
			assert(0);
		}
	}
};

class HalfOpenConeProcessorTropicalBasisTester :public HalfOpenConeProcessor
{
	PolynomialSet g;
public:
	  void process(HalfOpenCone const &c)
	  {
		  HalfOpenCone c2=c;
		  PolyhedralCone C=c2.closure();
	    	GroebnerFanTraverser traverser(groebnerBasisWithFullDimensionalIntersection(g,C),C);
	    	SymmetricTargetTropicalBasisTester target(g);
	    	symmetricTraverse(traverser,target);
	  }
	  HalfOpenConeProcessorTropicalBasisTester(PolynomialSet const &g_):
		  g(g_)
	  {

	  }
};

class TropicalIntersectionApplication : public GFanApplication
{

  SimpleOption optionTestIfTropicalBasis;
  SimpleOption optionTPlane;
  //  SimpleOption optionIncidencePrinting;
  SimpleOption optionParseSymmetry;
  SimpleOption optionExploitSymmetry;
  //  SimpleOption optionMinkowskiRefinement;
  SimpleOption optionIgnoreCones;
  SimpleOption optionRestrict;
//  SimpleOption optionXml;
  IntegerOption optionLow;
  IntegerOption optionHigh;
  SimpleOption optionStableIntersection;
public:
  const char *helpText()
  {
    return "This program computes the set theoretical intersection of a set of tropical hypersurfaces (or to be precise, their common refinement as a fan). The input is a list of polynomials with each polynomial defining a hypersurface. Considering tropical hypersurfaces as fans, the intersection can be computed as the common refinement of these. Thus the output is a fan whose support is the intersection of the tropical hypersurfaces.\n";
    //"The fan will be presented as a list of some of its closed cones. If a cone is a face of another cone in the fan it is not guaranteed to be listed. But the support of the fan will be the union of the listed cones.\n";
  }
  TropicalIntersectionApplication():
  //  optionXml("--xml","Produce a polymake file in XML format.\n"),
    optionTestIfTropicalBasis("-t","Note that the input polynomials generate an ideal. This option will make the program choose a relative interior point for each listed output cone and check if its initial ideal contains a monomial. The actual check is done on a homogenization of the input ideal, but this does not affect the result.\n"),
    optionTPlane("--tplane","This option intersects the resulting fan with the plane x_0=-1, where x_0 is the first variable. To simplify the implementation the output is actually the common refinement with the non-negative half space. This means that \"stuff at infinity\" (where x_0=0) is not removed."),
    optionRestrict("--restrict","Restrict the computation to a full-dimensional cone given by a list of marked polynomials. The cone is the closure of all weight vectors choosing these marked terms."),
    //optionIncidencePrinting("--incidence","Print incidence information of the fan. Only faces of maximal dimensional cones will be printed, so this works best if the fan is pure.")
    optionParseSymmetry("--symmetryPrinting","Parse a group of symmetries after the input has been read. Used when printing with --incidence."),
    //    optionMinkowskiRefinement("--minkowski","Compute the normal fan of the  Minkowski sum of the Newton polytopes instead.")
    optionExploitSymmetry("--symmetryExploit","Restrict computation to the closed lexicographic fundamental domain of the specified symmetry group. This overwrites --restrict."),
    optionIgnoreCones("--nocones","Tells the program not to output the CONES and MAXIMAL_CONES sections, but still output CONES_COMPRESSED and MAXIMAL_CONES_COMPRESSED if --symmetry is used."),
    optionHigh("--endcone","Specify interval [start,end[ of indices of the first fan (after appropriate reordering) which the computation is restricted to. Only with --restrict or --symmetryExploit. Useful for parallelizing a computation manually.",-1),
    optionLow("--startcone","Specify interval [start,end[ of indices of the first fan (after appropriate reordering) which the computation is restricted to. Only with --restrict or --symmetryExploit. Useful for parallelizing a computation manually.",-1),
    optionStableIntersection("--stable","Find the stable intersection of the input polynomials using tropical intersection theory. This can be slow. Most other options are ignored.")
  {
    registerOptions();
//    optionXml.hide();
    optionLow.hide();
    optionHigh.hide();
  }

  const char *name()
  {
    return "_tropicalintersection";
  }
  int main()
  {
    FileParser P(Stdin);

    PolynomialSet theInput=P.parsePolynomialSetWithRing();
    int n=theInput.numberOfVariablesInRing();

    if(optionStableIntersection.getValue())
      {
        PolyhedralFan f=PolyhedralFan::fullSpace(n);

        for(PolynomialSet::const_iterator i=theInput.begin();i!=theInput.end();i++)
        {
            PolyhedralFan f2=PolyhedralFan::normalFanOfNewtonPolytope(*i);
            if(f.size()==0)break;
            f=weilDivisor(f,f2);
        }
        f.printWithIndices(&pout,
                    FPF_multiplicities|
                    (optionParseSymmetry.getValue()?FPF_group|FPF_conesCompressed:0)|
                    (optionIgnoreCones.getValue()?0:FPF_conesExpanded)|
                           (optionTPlane.getValue()?FPF_boundedInfo|FPF_tPlaneSort:0)|
                    FPF_maximalCones|FPF_cones,0);
        return 0;
      }


    SymmetryGroup sym(n);
    if(optionParseSymmetry.getValue()||optionExploitSymmetry.getValue())sym.computeClosure(P.parseIntegerVectorList());

    PolyhedralFan F(n);
    /*    if(optionMinkowskiRefinement.getValue())
      {
	F=PolyhedralFan::fullSpace(n);
	for(PolynomialSet::const_iterator i=theInput.begin();i!=theInput.end();i++)
	  F=refinement(F,PolyhedralFan::normalFanOfNewtonPolytope(*i),n-1,false);
      }
      else*/

    if(optionTestIfTropicalBasis.getValue())
    {
    	HalfOpenConeProcessorTropicalBasisTester myProcessor(theInput);
    	tropicalHyperSurfaceIntersectionWithProcessor(n,theInput, myProcessor);
    }

    if(0)
      {
	F=tropicalPrincipalIntersection(n, theInput); // dimension of lineality space could be computed to speed up computations
      }
    else
      {
//	log1 fprintf(Stderr,"WARINING USING EXPERIMENTAL TROPICAL HYPERSURFACE INTERSECTION ROUTINE!!\n");

	if(optionRestrict.getValue()||optionExploitSymmetry.getValue())
	  {
	    IntegerVectorList inequalities;
       	    IntegerVectorList equations;

	    if(optionRestrict.getValue())
	      {
		PolynomialSet theConeAsPolys=P.parsePolynomialSet(theInput.getRing());
		inequalities=wallInequalities(theConeAsPolys);
	      }
	    else
	      {
		inequalities=sym.fundamentalDomainInequalities();
		equations=commonHomogeneitySpaceGenerators(theInput);
	      }
	    PolyhedralCone c(inequalities,equations,n);
	    c.canonicalize();

	    AsciiPrinter P(Stderr);
	    c.print(&P);

	    F=tropicalHyperSurfaceIntersectionClosed(n, theInput,&c,true,/*true*/false,optionLow.getValue(),optionHigh.getValue());//saveresult==false
	  }

	else
	  F=tropicalHyperSurfaceIntersectionClosed(n, theInput);
      }

    if(optionTPlane.getValue())
      {
	PolyhedralFan temp=PolyhedralFan::halfSpace(n,0);
	F=refinement(F,temp);
      }
    //    if(optionIncidencePrinting.getValue())
      {
	AsciiPrinter p(Stdout);
	PolyhedralFan a=F;
	//a.makePure();
	/////////a.printWithIndices(&p,false,&sym,false,false,optionXml.getValue(),optionTPlane.getValue());
	a.printWithIndices(&p,
		    FPF_multiplicities|
		    (optionParseSymmetry.getValue()?FPF_group|FPF_conesCompressed:0)|
		    (optionIgnoreCones.getValue()?0:FPF_conesExpanded)|
			   (optionTPlane.getValue()?FPF_boundedInfo|FPF_tPlaneSort:0)|
		    FPF_maximalCones|FPF_cones,
//			   FPF_default|
//			   (optionParseSymmetry.getValue()?FPF_group|FPF_conesCompressed:0) |
//			   (optionXml.getValue()?FPF_xml:0) |
//			   (optionTPlane.getValue()?FPF_boundedInfo|FPF_tPlaneSort:0),
			    &sym);
	//	a.printWithIndices(&p,false,&sym,optionParseSymmetry.getValue(),false,optionXml.getValue(),optionTPlane.getValue());
      }

    //AsciiPrinter(Stdout).printPolyhedralFan(F);

    //    AsciiPrinter Temp(Stdout);
    //    F.printWithIndices(&Temp,false,0);





/*    if(optionTestIfTropicalBasis.getValue())
      {
	fprintf(Stdout,"\nA list of relative interior points:\n");
	AsciiPrinter(Stdout).printVectorList(F.getRelativeInteriorPoints());

	PolynomialSet I=theInput;
	IntegerVector grading=IntegerVector::allOnes(I.numberOfVariablesInRing());

	PolynomialSet h=I.homogenization(I.getRing().withVariablesAppended("H"));

	IntegerVectorList r=F.getRelativeInteriorPoints();


	IntegerVectorList trueRays;
	IntegerVectorList falseRays;

	for(IntegerVectorList::const_iterator i=r.begin();i!=r.end();i++)
	  {
	    int n=h.numberOfVariablesInRing();
	    IntegerVector weight(n);
	    for(int j=0;j<n-1;j++)weight[j]=(*i)[j];
	    weight[n-1]=0;
	    PolynomialSet h2=h;
	    {
	      WeightReverseLexicographicTermOrder t(weight);
	      fprintf(Stdout,"Computing the initial ideal with respect to:");
	      AsciiPrinter(Stdout).printVector(weight);
	      fprintf(Stdout,"\n");
	      buchberger(&h2,t);
	      fprintf(Stdout,"Done computing the initial ideal.\n");
	    }

	    PolynomialSet wall=initialFormsAssumeMarked(h2,weight);

	    if(containsMonomial(wall))
	      {
		fprintf(Stdout,"The (homogenized) initial ideal contains a monomial: ");
		AsciiPrinter(Stdout).printPolynomial(Polynomial(computeTermInIdeal(wall)));
		fprintf(Stdout,"\n");
		falseRays.push_back(*i);
	      }
	    else
	      {
		fprintf(Stdout,"The initial ideal contains no monomial.\n");
		trueRays.push_back(*i);
	      }
	    fprintf(Stdout,"\n");
	  }
	fprintf(Stdout,"The set of tested interior points that are in the tropical variety of the ideal generated by the input:\n");
	AsciiPrinter(Stdout).printVectorList(trueRays);
	fprintf(Stdout,"The set of tested interior points that are not in the tropical variety of the ideal generated by the input:\n");
	AsciiPrinter(Stdout).printVectorList(falseRays);
      }*/
    return 0;
  }
};

static TropicalIntersectionApplication theApplication;
