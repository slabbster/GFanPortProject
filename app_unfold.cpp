#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "lp.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"

using namespace std;

class UnfoldApplication : public GFanApplication
{
	class Edge{
	public:
		int i,j;
	};

	class Facet{
	  vector<Edge> edges;
	public:

	};

	class Surface{
		vector<Edge> edges;
		vector<Facet> facets;
	public:

	};


public:
  bool includeInDefaultInstallation() // Not included since the program has no relation to the main programs
  {
    return false;
  }
  UnfoldApplication():
	  input1Option("-i1","Specify the name of the first input file.","polymake.out")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_unfold";
  }
  int main()
  {
    FileParser P(Stdin);

    PolyhedralFan f1=PolyhedralFan::readFan(input1Option.getValue());

    assert(f1.getAmbientDimension()==4);

    int boxSize=2;
    IntegerVectorList equations;
	IntegerVectorList inequalities;
	inequalities.push_back(boxSize*IntegerVector::standardVector(4,0)+IntegerVector::standardVector(n,1));
	inequalities.push_back(boxSize*IntegerVector::standardVector(4,0)+IntegerVector::standardVector(n,2));
	inequalities.push_back(boxSize*IntegerVector::standardVector(4,0)+IntegerVector::standardVector(n,3));

	PolyhedralCone C(inequalities,equalities,4);
	C.canonicalize();
	PolyhedralFan F(4);
	F.insert(C);

	PolyhedralFan f2=refinement(f1,F);

	IntegerVectorList rays=f2.getRays();

	Surface s;

	for(PolyhedralFan::coneIterator i=f2.conesBegin();i!=f2.conesEnd();i++)
	{
		if(i->dimension()==3)
		{
			PolyhedralFan f3=PolyhedralFan::facetsOfCone(*i);
			for(PolyhedralFan::coneIterator i=f3.conesBegin();i!=f3.conesEnd();i++)
			{
				Facet F;
				int J=0;
				for(IntegerVectorList::const_iterator j=rays.begin();j!=rays.end();j++,J++)
					if(i->contains(*j))
						F.edges.push_back(J);
			}
		}
		s.facets.push_back(F);
	}


    IntegerVectorList v=P.parseIntegerVectorList();
    fprintf(Stderr,"Rank:%i\n",rankOfMatrix(v));
    AsciiPrinter(Stdout).printVectorList(transposeIntegerVectorList(v));
    return 0;
  }
  const char *helpText()
  {
    return "\n";
  }
};

static UnfoldApplication theApplication;
