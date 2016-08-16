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
#include "polyhedralcone.h"
#include "gfanapplication.h"
#include "saturation.h"
#include "symmetrictraversal.h"
#include "traverser_tropical.h"

#define FILENAME "_templatex"

static void printVertex(Printer &p, const char *comment, PolynomialSet const &s)
{
  p.printString("{\\bf ");
  p.printString(comment);
  p.printString("}\n");
  p.printPolynomialSet(s,true);
  p.printNewLine();
}


class InteractiveApplication : public GFanApplication
{
  SimpleOption optionLatex;
  SimpleOption optionExitImmediately;
  SimpleOption optionPrintFlipped;
  SimpleOption optionPrintWallIdeal;
  SimpleOption optionPrintInequalities;
  SimpleOption optionPrintWeightVector;
  SimpleOption optionTropicalVariety;
public:
  const char *helpText()
  {
    return "This is a program for doing interactive walks in the Groebner fan of an ideal. "
      "The input is a Groebner basis defining the starting Groebner cone of the walk. "
      "The program will list all flippable facets of the Groebner cone and ask the user to choose one. "
      "The user types in the index (number) of the facet in the list. "
      "The program will walk through the selected facet and display the new Groebner basis and a list of new facet normals for the user to choose from. "
      "Since the program reads the user's choices through the the standard input it is recommended not to redirect the standard input for this program.\n";
  }
  InteractiveApplication():
    optionLatex("-L","Latex mode. The program will try to show the current Groebner basis in a readable form by invoking LaTeX and xdvi.\n"),
    optionExitImmediately("-x","Exit immediately.\n"),
    optionPrintWallIdeal("-w","Tell the program to list (a Groebner basis with respect to the current term order for) the initial ideal for each flippable wall in the current Groebner cone.\n"),
    optionPrintFlipped("-f","Tell the program to list the flipped reduced Groebner basis of the initial ideal for each flippable wall in the current Groebner cone.\n"),
    optionPrintInequalities("-i","Tell the program to list the defining set of inequalities of the non-restricted Groebner cone as a set of vectors after having listed the current Groebner basis.\n"),
    optionPrintWeightVector("-W","Print weight vector. This will make the program print an interior vector of the current Groebner cone and a relative interior point for each flippable facet of the current Groebner cone.\n"),
    optionTropicalVariety("--tropical","Traverse a tropical variety interactively.")
    {
    registerOptions();
  }

  const char *name()
  {
    return "_interactive";
  }


  int tropical()
  {
	  AsciiPrinter P(Stdout);
	  PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();
	  PolynomialSet g2=FileParser(Stdin).parsePolynomialSet(g.getRing());
	  debug<<"test1\n";
	  TropicalTraverser traverser(g,g2);
	  while(1)
	  {
debug<<"test2\n";
		  P<<"test3\n";
		  P<< traverser.refToPolyhedralCone();
		  IntegerVector ridge=FileParser(Stdin).parseIntegerVector();
		  P<<"RAYS:"<< traverser.link(ridge);
		  IntegerVector ray=FileParser(Stdin).parseIntegerVector();
		  traverser.changeCone(ridge,ray);
	  }
	  return 0;
  }

  int main()
  {
	  if(optionTropicalVariety.getValue())return tropical();
	  LexicographicTermOrder myOrder;

    bool outputLatex=optionLatex.getValue();

    Printer *P;
    LatexPrinter *Q=0;
    FILE *latexFile=0;

    PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();
    PolynomialRing theRing=g.getRing();
    PolynomialSet gOld(theRing);


    //    buchberger(&g,myOrder);
    g.scaleMarkedCoefficientsToOne();
    fprintf(Stderr,"Minimizing and autoreducing input...\n");
    minimize(&g);
    autoReduce(&g, LexicographicTermOrder());

    if(outputLatex)system("xdvi " FILENAME ".dvi&");

    gOld=g;

    while(1)
      {
	if(outputLatex)
	  {
	    latexFile=fopen(FILENAME ".tex","w");
	    Q = new LatexPrinter(latexFile);
	    Q->printLatexStart();
	    P=Q;
	  }
	else
	  P= new AsciiPrinter(Stdout);


	//	P->printString("-------------------------------------------------\n");
	printVertex(*P,"Old Vertex",gOld);
	printVertex(*P,"Current Vertex",g);

	//	P->printPolynomialSet(g.markedTermIdeal());

	P->printNewLine();

	PolynomialSet wall(theRing);
	IntegerVectorList normals=wallInequalities(g);
	IntegerVectorList flipableNormals=wallFlipableNormals(g,false);


	{
	  wall=wallIdeal(g,*flipableNormals.begin());
	  PolynomialSet sat=nonHomogeneousSaturation(wall);
	  P->printString("Saturated wall ideal: ");
	  P->printPolynomialSet(sat);
	}

	if(optionPrintInequalities.getValue())
	  {
	    P->printString("Inequalities:");
	    P->printNewLine();
	    P->printVectorList(normals);
	  }

	if(optionPrintWeightVector.getValue())
	  {
	    P->printString("An interior point: ");
	    PolyhedralCone C=PolyhedralCone(normals,IntegerVectorList());
	    P->printVector(intersection(PolyhedralCone::positiveOrthant(C.ambientDimension()),C).getRelativeInteriorPoint());
	    P->printNewLine();
	    //interiorPoint(normals,Stdout,true);
	  }

	if(outputLatex)
	  P->printString("\\begin{enumerate}");

	IntegerVectorList facets;

	int itemIndex=0;
	for(IntegerVectorList::const_iterator i=flipableNormals.begin();i!=flipableNormals.end();i++)
	  {
	    //	    if(isFacet(normals,i))
	    //  if(wallContainsPositiveVector(*i))

		{
		  itemIndex++;
		  facets.push_back(*i);
		  if(outputLatex)
		    P->printString("\\item ");
		  else
		    {
		      P->printString("--");
		      P->printInteger(itemIndex);
		    }
		  wall=wallIdeal(g,*i);


		  WeightTermOrder flipOrder(*i);

		  wall.markAndScale(flipOrder);
		  PolynomialSet oldWall=wall;

		  P->printNewLine();
		  P->printString("Facet normal: ");
		  P->printVector(*i);

		  if(optionPrintWeightVector.getValue())
		    {
		      P->printNewLine();
		      P->printString("Relative interior point: ");

		      IntegerVectorList l;
		      l.push_back(*i);
		      //		    AsciiPrinter(Stdout).printVector(PolyhedralCone(normals,l).relativeInteriorPoint(true));
		      P->printVector(intersection(PolyhedralCone(normals,l),PolyhedralCone::positiveOrthant(i->size())).getRelativeInteriorPoint());
		    }

		  //int numberOfPolynomials=wall.size();
		  //P->printString("Number of polynomials\n");
		  //P->printInteger(numberOfPolynomials);

		  //int numberOfMonomials=0;
		  //for(PolynomialSet::const_iterator j=wall.begin();j!=wall.end();j++)
		    //  if(j->isMonomial())numberOfMonomials++;

		  //		  P->printString("Number of monomials\n");
		  //P->printInteger(numberOfMonomials);

		  P->printNewLine();

		  if(optionPrintWallIdeal.getValue())
		    {
		      //P->printString("Gr\\\"obner basis of wall");
		      P->printString("Wall ideal: ");
		      P->printPolynomialSet(wall);
		      P->printNewLine();
		      PolynomialSet sat=nonHomogeneousSaturation(wall);
		      P->printString("Saturated wall ideal: ");
		      P->printPolynomialSet(sat);
		    }

		  if(outputLatex||optionPrintFlipped.getValue())
		    {
		      flipOrder=WeightTermOrder(-*i);

		      buchberger(&wall,flipOrder);

		      P->printNewLine();
		      P->printString("New Gr\\\"obner basis of wall");
		      P->printPolynomialSet(wall);
		    }
		}
	  }
	if(outputLatex)
	  P->printString("\\end{enumerate}");

	if(outputLatex)
	  {
	    Q->printLatexEnd();
	    fclose(latexFile);

	    system("latex "FILENAME".tex >/dev/null");
	  }
	{
	  if(optionExitImmediately.getValue())
	    {
	      fprintf(Stderr,"exiting.\n");
	      break;
	    }
	  int t;
	  fscanf(Stdin,"%i",&t);
	  t--;
	  if(t<0)break;
	  IntegerVectorList::const_iterator i=facets.begin();
	  while(t>0)
	    {
	      i++;
	      t--;
	    }
	  gOld=g;
	  g=flip(g,*i);
	}
      }
    return 0;
  }
};

static InteractiveApplication theApplication;

