#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "lp.h"
#include "gfanapplication.h"
#include "polyhedralcone.h"

#include "polymakefile.h"
#include "determinant.h"
#include "subspace.h"
#include "triangulation.h"

class TriangulateApplication : public GFanApplication
{
  FieldOption theFieldOption;
  StringOption inputOption;
  StringOption outputOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program ........\n";
  }
  TriangulateApplication():
    inputOption("-i","Specify the name of the input file.","polymake.out"),
    outputOption("-o","Specify the name of the output file.","polymake.out2")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_triangulate";
  }

  void printIntList(list<int> const &v)
  {
    FILE *f=Stderr;
    fprintf(f,"{");
    for(list<int>::const_iterator i=v.begin();i!=v.end();i++)
      {
	if(i!=v.begin())fprintf(f," ");
	fprintf(f,"%i",*i);
      }
    fprintf(f,"}\n");
  }
  void printIntListList(list<list<int> > const &l)
  {
    for(list<list<int> >::const_iterator i=l.begin();i!=l.end();i++)
      printIntList(*i);
  }


  int main()
  {
    LpSolver::printList(Stderr);
    lpSetSolver("cddgmp");

    PolymakeFile inFile;
    fprintf(Stderr,"Test\n");
    inFile.open(inputOption.getValue());
    fprintf(Stderr,"Test\n");


    int n=inFile.readCardinalProperty("AMBIENT_DIM");
    //int d=inFile.readCardinalProperty("DIM");
    int nRays=inFile.readCardinalProperty("N_RAYS");

    fprintf(Stderr,"%i %i\n",n,nRays);

    IntegerMatrix rays=inFile.readMatrixProperty("RAYS",nRays,n);

    {/* Make sure that representatives of rays are in the same subspace. */
      IntegerVectorList rays2;
      IntegerMatrix linealitySpace=inFile.readMatrixProperty("LINEALITY_SPACE",nRays,n);
      Subspace l(linealitySpace.getRows(),linealitySpace.getWidth());
      for(int i=0;i<rays.getHeight();i++)
	rays2.push_back(l.canonicalizeVector(rays[i]));
      rays=rowsToIntegerMatrix(rays2,n);
    }


    /*    if(1)//just for printing - no triangulation here.
      {
	//	int m=inFile.readCardinalProperty("NMAXCONES");
	vector<list<int> > cones=inFile.readMatrixIncidenceProperty("MAXIMAL_CONES_COMPRESSED");
	IntegerVectorList r;

	for(int i=0;i<cones.size();i++)
	  {
	    IntegerVector v(n);
	    for(list<int>::const_iterator j=cones[i].begin();j!=cones[i].end();j++)
	      {

		//AsciiPrinter(Stderr).printVector(v);fprintf(Stderr,"%i\n",*j);
		v+=rays[*j];
		//AsciiPrinter(Stderr).printVector(rays[*j]);fprintf(Stderr,"%i\n",*j);
	      }
	    r.push_back(v);
	  }


	PolymakeFile outFile;
	outFile.create(outputOption.getValue(),"topaz","SimplicialComplex");

	outFile.writeMatrixProperty("INTERIOR_POINTS",rowsToIntegerMatrix(r));
	outFile.close();

	}*/


    //    vector<Triangulation::Cone > cones=inFile.readMatrixIncidenceProperty("MAXIMAL_CONES");
    vector<list<int> > cones=inFile.readMatrixIncidenceProperty("MAXIMAL_CONES");

    for(vector<list<int> >::const_iterator i=cones.begin();i!=cones.end();i++)
      printIntList(*i);

    AsciiPrinter(Stderr).printVectorList(rays.getRows());

    vector<Triangulation::Cone> simplicialComplex;

    //    for(vector<Triangulation::Cone>::const_iterator i=cones.begin();i!=cones.end();i++)
    for(vector<list<int> >::const_iterator i=cones.begin();i!=cones.end();i++)
      {
	//	fprintf(Stderr,"Triangulating:\n");
	//	printIntList(*i);
	//	list<Triangulation::Cone> coneList=Triangulation::triangulate(*i,rays);
	list<Triangulation::Cone> coneList=Triangulation::triangulate(*i,rays);
	// fprintf(Stderr,"Result:\n");
	for(list<Triangulation::Cone>::const_iterator j=coneList.begin();j!=coneList.end();j++)
	  {
	    // printIntList(*j);
	    simplicialComplex.push_back(*j);
	  }
      }

    PolymakeFile outFile;
    outFile.create(outputOption.getValue(),"topaz","SimplicialComplex");

    outFile.writeCardinalProperty("N_VERTICES",nRays);
    //    outFile.writeCardinalProperty("DIM",d);

    outFile.writeIncidenceMatrixProperty("INPUT_FACES",Triangulation::removeOrientation(simplicialComplex),nRays);

    outFile.close();

    return 0;
  }
};

static TriangulateApplication theApplication;

