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
#include "polyhedralfan.h"
#include "linalg.h"
#include "polymakefile.h"
#include "tropical_weildivisor.h"
#include "log.h"
#include "symmetry.h"

class TropicalWeilDivisorApplication : public GFanApplication
{
  StringOption inputOption1;
  StringOption inputOption2;
public:
  const char *helpText()
  {
    return "This program computes the tropical Weil divisor of piecewise linear (or tropical rational) function on a tropical k-cycle. See the Gfan manual for more information.\n";
  }
  TropicalWeilDivisorApplication():
    inputOption1("-i1","Specify the name of the Polymake input file containing the k-cycle.","polymake1.out"),
    inputOption2("-i2","Specify the name of the Polymake input file containing the piecewise linear function.","polymake2.out")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_tropicalweildivisor";
  }

  IntegerVectorList readLinealitySpace(const char *filename)
  {
    PolymakeFile inFile;
    inFile.open(filename);
    int n=inFile.readCardinalProperty("AMBIENT_DIM");
    int homog=inFile.readCardinalProperty("LINEALITY_DIM");
    IntegerMatrix rays=inFile.readMatrixProperty("LINEALITY_SPACE",homog,n);
    return rays.getRows();
  }
  IntegerVectorList readRays(const char *filename)
  {
    PolymakeFile inFile;
    inFile.open(filename);
    int n=inFile.readCardinalProperty("AMBIENT_DIM");
    int nRays=inFile.readCardinalProperty("N_RAYS");
    IntegerMatrix rays=inFile.readMatrixProperty("RAYS",nRays,n);
    return rays.getRows();
  }
  PolyhedralFan readFan(const char *filename, bool readRayValues, bool readLinearForms, bool readMultiplicities)
  {
    PolymakeFile inFile;
    inFile.open(filename);

    int n=inFile.readCardinalProperty("AMBIENT_DIM");
    int nRays=inFile.readCardinalProperty("N_RAYS");
    int lindim=inFile.readCardinalProperty("LINEALITY_DIM");
    IntegerMatrix rays=inFile.readMatrixProperty("RAYS",nRays,n);
    vector<list<int> > cones=inFile.readMatrixIncidenceProperty("MAXIMAL_CONES");

    IntegerMatrix multiplicities;
    if(readMultiplicities)
      multiplicities=inFile.readMatrixProperty("MULTIPLICITIES",cones.size(),1);

    PolyhedralFan ret(n);

    IntegerVectorList equations=inFile.readMatrixProperty("LINEALITY_SPACE",lindim,n).getRows();

    IntegerMatrix linearForms(cones.size(),n);
    if(readLinearForms)linearForms=inFile.readMatrixProperty("LINEAR_FORMS",cones.size(),n);

    IntegerMatrix rayValues(nRays,1);
    if(readRayValues)rayValues=inFile.readMatrixProperty("RAY_VALUES",nRays,1);

    IntegerMatrix linealityValues(lindim,1);
    if(readRayValues)linealityValues=inFile.readMatrixProperty("LINEALITY_VALUES",lindim,1);


    int I=0;
    for(vector<list<int> >::const_iterator i=cones.begin();i!=cones.end();i++,I++)
      {
	IntegerVectorList coneGenerators;
	FieldMatrix LFA(Q,lindim+i->size(),n);
	FieldVector LFB(Q,lindim+i->size()+n);
	int J=0;
	for(list<int>::const_iterator j=i->begin();j!=i->end();j++,J++)
	  {
	    coneGenerators.push_back(rays[*j]);
	    LFA[J]=integerVectorToFieldVector(rays[*j],Q);
	    LFB[J]=Q.zHomomorphism(rayValues[*j][0]);
	  }
	int K=0;
	for(IntegerVectorList::const_iterator j=equations.begin();j!=equations.end();j++,J++,K++)
	  {
	    LFA[J]=integerVectorToFieldVector(*j,Q);
	    LFB[J]=Q.zHomomorphism(linealityValues[K][0]);
	  }
	if(readRayValues)
	  {
	    AsciiPrinter P(Stderr);//log0 P<<LFA<<LFB;

	    FieldVector LFX=LFA.solver().canonicalize(LFB);
	    if(LFX.subvector(0,LFX.size()-n).isZero())
	      {
		linearForms[I]=fieldVectorToIntegerVector(LFX.subvector(LFX.size()-n,LFX.size()));
	      }
	    else
	      {
		cerr<<"Values on cone are not linear" <<endl;
		assert(0);
	      }
	  }
	//AsciiPrinter P(Stderr);
	PolyhedralCone d(coneGenerators,equations,n);
	d.canonicalize();
	//d.print(P);
	IntegerVectorList inequalities=d.extremeRays();
	IntegerVectorList equations2;
	IntegerMatrix A=rowsToIntegerMatrix(coneGenerators,n);
	IntegerMatrix B=rowsToIntegerMatrix(equations,n);///!!!!
	FieldMatrix M=combineOnTop(integerMatrixToFieldMatrix(A,Q),integerMatrixToFieldMatrix(B,Q));
	FieldMatrix kerM=M.reduceAndComputeKernel();
	for(int i=0;i<kerM.getHeight();i++)
	  {
	    equations2.push_back(kerM[i].primitive());
	  }
	PolyhedralCone c(inequalities,equations2,n);
	c.canonicalize();
	if(readLinearForms || readRayValues)c.setLinearForm(linearForms[I]);

	if(readMultiplicities)
	  c.setMultiplicity(multiplicities[I][0]);

	ret.insert(c);
      }
    return ret;
  }

  int main()
  {
    //    LpSolver::printList(Stderr);
    //    lpSetSolver("cddgmp");


    PolyhedralFan f=readFan(inputOption1.getValue(),false,false,true);
    PolyhedralFan g=readFan(inputOption2.getValue(),true,false,false);

    IntegerVectorList originalRays=readRays(inputOption1.getValue());
    IntegerVectorList linealitySpace=readLinealitySpace(inputOption1.getValue());
    /*    FileParser P(Stdin);

    PolynomialRing R=P.parsePolynomialRing();

    Polynomial p=P.parsePolynomial(R);
    */

    PolyhedralFan divisor=weilDivisor(f,g);

    AsciiPrinter Q(stdout);

    SymmetryGroup sym(f.getAmbientDimension());
    vector<string> comments=f.renamingStrings(divisor.getRaysInPrintingOrder(&sym),originalRays,linealitySpace,&sym);
    //    divisor.printWithIndices(&Q,true,0,false,false,false,false,&comments);
    divisor.printWithIndices(&Q,FPF_default|FPF_multiplicities,0,&comments);

    /*    FileParser P(Stdin);

    PolynomialRing R=P.parsePolynomialRing();
    PolynomialSet s=P.parsePolynomialSet(R);

    Polynomial sum(R);
    for(PolynomialSet::const_iterator i=s.begin();i!=s.end();i++)
      sum+=*i;

    AsciiPrinter(Stdout).printPolynomial(sum);
    */


    return 0;
  }
};

static TropicalWeilDivisorApplication theApplication;
