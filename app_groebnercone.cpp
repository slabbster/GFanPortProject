#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "wallideal.h"
#include "lp.h"
#include "polyhedralcone.h"
#include "gfanapplication.h"
#include "polyhedralfan.h"
#include "halfopencone.h"
#include "linalg.h"
#include "log.h"

class GroebnerConeApplication : public GFanApplication
{
  SimpleOption optionRestrict;
  SimpleOption optionPair;
  SimpleOption optionAsFan;
  SimpleOption optionVectorInput;
public:
  const char *helpText()
  {
    return "This program computes a Groebner cone. Three different cases are handled. The input may be a marked reduced Groebner basis in which case its Groebner cone is computed. The input may be just a marked minimal basis in which case the cone computed is not a Groebner cone in the usual sense but smaller. (These cones are described in [Fukuda, Jensen, Lauritzen, Thomas]). The third possible case is that the Groebner cone is possibly lower dimensional and given by a pair of Groebner bases as it is useful to do for tropical varieties, see option --pair. The facets of the cone can be read off in section FACETS and the equations in section IMPLIED_EQUATIONS.\n";
  }
  GroebnerConeApplication():
    optionRestrict("--restrict","Add an inequality for each coordinate, so that the the cone is restricted to the non-negative orthant."),
    optionPair("--pair","The Groebner cone is given by a pair of compatible Groebner bases. The first basis is for the initial ideal and the second for the ideal itself. See the tropical section of the manual."),
    optionAsFan("--asfan","Writes the cone as a polyhedral fan with all its faces instead. In this way the extreme rays of the cone are also computed."),
    optionVectorInput("--vectorinput","Compute a cone given list of inequalities rather than a Groebner cone. The input is an integer which specifies the dimension of the ambient space, a list of inequalities given as vectors and a list of equations.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_groebnercone";
  }

  int main()
  {
    IntegerVectorList equalities;
    IntegerVectorList normals;
    int n;

    if(optionVectorInput.getValue())
      {
    	FileParser P(Stdin);
    	n=P.parseInt();
    	normals=P.parseIntegerVectorList();
    	equalities=P.parseIntegerVectorList();
      }
    else
    {
      PolynomialSet m=FileParser(Stdin).parsePolynomialSetWithRing();
      PolynomialSet g(m.getRing());
      if(optionPair.getValue())
      {
    	  g=FileParser(Stdin).parsePolynomialSet(m.getRing());
      }
      else
      {
    	  g=m;
    	  m=g.markedTermIdeal();
      }

      {
    	  //Check that markings are consistent
    	  PolynomialSet M1=g.markedTermIdeal();
    	  PolynomialSet M2=m.markedTermIdeal();
    	  assert(M1.size()==M2.size());
    	  PolynomialSet::const_iterator i1=M1.begin();
    	  for(PolynomialSet::const_iterator i2=M2.begin();i2!=M2.end();i1++,i2++)
    	  {
    		  assert((*i1-*i2).isZero());
    	  }
      }

      n=g.getRing().getNumberOfVariables();

      equalities=wallInequalities(m);
      normals=(wallInequalities(g));
      if(optionRestrict.getValue())
      {
    	  for(int i=0;i<n;i++)
    		  normals.push_back(IntegerVector::standardVector(n,i));
      }
    }



    {
      FieldMatrix A=integerMatrixToFieldMatrix(rowsToIntegerMatrix(equalities,n),Q);
      A.reduce();
      A.removeZeroRows();
      equalities=IntegerVectorList();
      for(int i=0;i<A.getHeight();i++)
	equalities.push_back(A[i].primitive());

      set<IntegerVector> newInequalities;
      for(IntegerVectorList::const_iterator i=normals.begin();i!=normals.end();i++)
	newInequalities.insert(A.canonicalize(integerVectorToFieldVector(*i,Q)).primitive());

      normals=IntegerVectorList();
      for(set<IntegerVector>::const_iterator i=newInequalities.begin();i!=newInequalities.end();i++)
	normals.push_back(*i);
    }


    log1 fprintf(Stderr,"Inequalities");
    log1 AsciiPrinter(Stderr).printVectorList(normals);
    log1 fprintf(Stderr,"Equations");
    log1 AsciiPrinter(Stderr).printVectorList(equalities);

    PolyhedralCone c(normals,equalities,n);
    c.canonicalize();
    if(!optionAsFan.getValue())
      AsciiPrinter(Stdout).printPolyhedralCone(c);
    else
      {
	AsciiPrinter P(Stdout);
	IntegerVectorList empty;
	HalfOpenCone C(n,c.getEquations(),c.getHalfSpaces(),empty,true);
	//	PolyhedralFan F=faceComplexOfCone(C);
	PolyhedralCone C2=C.closure();
	C2.canonicalize();
	PolyhedralFan F(n);F.insert(C2);
	F.printWithIndices(&P,FPF_default);
	//	F.printWithIndices(&P,false,0,false,false,optionXml.getValue());
      }

    return 0;
  }
};

static GroebnerConeApplication theApplication;

