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

#include "symmetry.h"

class CombineRaysApplication : public GFanApplication
{
  StringOption inputOption;
  StringOption sectionOption;
  SimpleOption dumpOption;
  IntegerOption dumpHeightOption;
  IntegerOption dumpWidthOption;
public:
  const char *helpText()
  {
    return "This program combines rays from the specified polymake file by adding according to a list of vectors of indices given on the standard input.\n";
  }
  CombineRaysApplication():
    inputOption("-i","Specify the name of the input file.","examples/grassmann3_7.out"),
	sectionOption("--section","Specify a section of the polymake file to use as input, rather than standard input.",0),
	dumpOption("--dump","Dump specified section as a matrix rather than combining the rays"),
        dumpHeightOption("--dheight","Specify height of matrix to be dumped.",1),
        dumpWidthOption("--dwidth","Specify width of matrix to be dumped.",1)
  {
    dumpOption.hide();
    dumpHeightOption.hide();
    dumpWidthOption.hide();
    registerOptions();
  }

  const char *name()
  {
    return "_combinerays";
  }

  int main()
  {
    PolymakeFile inFile;
    inFile.open(inputOption.getValue());

    if(dumpOption.getValue())
      {
        pout<<inFile.readMatrixProperty(sectionOption.getValue(),dumpHeightOption.getValue(),dumpWidthOption.getValue()).getRows();
        return 0;
      }

    int N=inFile.readCardinalProperty("AMBIENT_DIM");

    int nRays=inFile.readCardinalProperty("N_RAYS");
    IntegerMatrix rays=inFile.readMatrixProperty("RAYS",nRays,N);


    FileParser P(Stdin);

    IntegerVectorList comb;
    if(sectionOption.getValue())
    {
		vector<list<int> > l=inFile.readMatrixIncidenceProperty(sectionOption.getValue());
		for(vector<list<int> >::const_iterator i=l.begin();i!=l.end();i++)
		{
			IntegerVector temp(i->size());
			int J=0;
			for(list<int>::const_iterator j=i->begin();j!=i->end();j++,J++)temp[J]=*j;
			comb.push_back(temp);
		}
    }
    else
    	comb=P.parseIntegerVectorList();

    IntegerVectorList result;
    for(IntegerVectorList::const_iterator j=comb.begin();j!=comb.end();j++)
      {
	IntegerVector interiorPoint(N);
	for(int i=0;i<j->size();i++)
	  {
	    interiorPoint+=rays[(*j)[i]];
	  }
	result.push_back(interiorPoint);
      }

    AsciiPrinter(Stdout).printVectorList(result);
    fprintf(Stdout,"\n");

    return 0;
  }
};

static CombineRaysApplication theApplication;
