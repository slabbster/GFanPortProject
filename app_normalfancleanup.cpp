#include <iostream>
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
#include "log.h"


class NormalFanCleanUpApplication : public GFanApplication
{
  SimpleOption optionPreprocess;

public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program takes a list of polynomials and outputs a subcomplex of its normalfan. After the polynomial list follows a generating set for the symmetry group and then a list of vectors. The cones which are output are the ones picked by the list of vectors. The output is symmetric with respect to the symmetry group.\n";
  }
  NormalFanCleanUpApplication():
    optionPreprocess("-pre","Preprocess cones before building polyhedral fan - that is find representatives for each orbit first.\n")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_normalfancleanup";
  }
  int main()
  {
    FileParser P(Stdin);

    PolynomialSet g=P.parsePolynomialSetWithRing();
    int n=g.numberOfVariablesInRing();

    SymmetryGroup sym(n);
    sym.computeClosure(P.parseIntegerVectorList());

    IntegerVectorList l=P.parseIntegerVectorList();

    cerr<<"Number of input rays: " << l.size() <<"\n";

    if(optionPreprocess.getValue())
      {
	set<IntegerVector> l2;
	for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
	  {
	    {
	      static int a;
	      if(!((++a)&255))
		{
		  cerr<<a<<"\n";

		  FILE *f=fopen("iteraTIon","w");
		  if(f)
		    {
		      fprintf(f,"%i:%i\n",a,l2.size());
		      fclose(f);
		    }
		}
	    }
	    {
	      static int a;
	      if(!((++a)&2047))
		{
		  FILE *f=fopen("parTialOutPUt","w");
		  if(f)
		    {
		      IntegerVectorList l3;
		      for(set<IntegerVector>::const_iterator i=l2.begin();i!=l2.end();i++)l3.push_back(*i);
		      AsciiPrinter(f)<<l3;
		      fclose(f);
		    }
		}
	    }
	    PolyhedralCone c=normalConeOfMinkowskiSum(g,*i);
	    c.canonicalize();
	    cerr << "Dim " <<c.dimension()<<endl;
	    l2.insert(sym.orbitRepresentative(c.getUniquePoint()));
	  }

	cerr<<"Number of input cones up to symmetry: " << l2.size() <<"\n";


	IntegerVectorList l3;
	for(set<IntegerVector>::const_iterator i=l2.begin();i!=l2.end();i++)l3.push_back(*i);
	AsciiPrinter(Stdout)<<l3;

	l=l3;
      }


    PolyhedralFan F(n);
    for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
      {
	{
	  static int a;
	  if(!((++a)&255))
	    {
	      cerr<<a<<"\n";

	      FILE *f=fopen("iteraTIon","w");
	      if(f)
		{
		  fprintf(f,"%i\n",a);
		  fclose(f);
		}
	    }
	}
	PolyhedralCone c=normalConeOfMinkowskiSum(g,*i);
	c.canonicalize();
	F.insert(c);
      }
    AsciiPrinter p(Stdout);
    F.printWithIndices(&p,
		       FPF_group|FPF_conesCompressed|
		       FPF_maximalCones|FPF_cones,
		       &sym);


    return 0;
  }
};

static NormalFanCleanUpApplication theApplication;
