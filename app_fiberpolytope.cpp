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
#include "gfanapplication.h"
#include "timer.h"
#include "log.h"
#include "matrix.h"
#include "lll.h"
#include "polyhedralfan.h"
#include "linalg.h"
#include "determinant.h"
#include "triangulation.h"
#include "intsinpolytope.h"
#include "graph.h"

#include "triangulation2.h"

#include "traverser_secondaryfan.h"
#include "symmetrictraversal.h"

#include <iostream>
#include <algorithm>

class FiberPolytopeApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program computes the normal fan of a fiber polytope. The input is a list of vertices of a polytope, followed by a list linear map given by a matrix.\n";
  }
  FiberPolytopeApplication()
  {
     registerOptions();
  }

  const char *name()
  {
    return "_fiberpolytope";
  }

  int main()
  {
	    IntegerMatrix A=rowsToIntegerMatrix(FileParser(Stdin).parseIntegerVectorList());
	    IntegerMatrix map=rowsToIntegerMatrix(FileParser(Stdin).parseIntegerVectorList());

	    IntegerMatrix transformedPoints=(A*map).transposed();
	    IntegerVectorList temp=transformedPoints.getRows();
	    temp.push_back(IntegerVector::allOnes(transformedPoints.getWidth()));
	    IntegerMatrix B=rowsToIntegerMatrix(temp);

	    /* If the vector configuration B does not have full rank then
	       change coordinates. */
	    if(rank(B)!=B.getHeight())
	      {
	        FieldMatrix M=integerMatrixToFieldMatrix(B,Q);
	        M.reduce(false,true);//force integer operations - preserving volume
	        M.removeZeroRows();
	        B=fieldMatrixToIntegerMatrix(M);
	      }



	    AsciiPrinter P(Stdout);
	    P<<B.transposed().getRows();

	    IntegerVectorList empty;

	    log1 debug<<"Computing restricting cone\n";
	    PolyhedralCone C=PolyhedralCone::givenByRays(empty,A.transposed().getRows(),A.getHeight());
	    log1 debug<<"Done computing restricting cone\n";

	    debug<<B.getRows();
	    debug<<int(rank(B));

	    Triangulation2 t(B);
log1 debug<<"Computing initial triangulation\n";
	    // Convert a Triangulation to a Triangulation2
	    {
	      list<Triangulation::Cone> T=Triangulation::triangulate(B.transposed());
	      for(list<Triangulation::Cone>::const_iterator i=T.begin();i!=T.end();i++)
		{
		  IntegerVector v=i->size();
		  int J=0;
		  for(Triangulation::Cone::const_iterator j=i->begin();j!=i->end();j++,J++)
		    v[J]=*j;
		    t.bases.insert(v);
		}
	    }
	    log1 debug<<"Done computing initial triangulation\n";
	    AsciiPrinter p(Stdout);
	    int n=B.getWidth();
	    SymmetryGroup s(n);
	    SymmetricTargetFanBuilder target(n,s);

	    log1 debug<<"Initializing\n";
	    SecondaryFanTraverser traverser(triangulationWithFullDimensionalIntersection(t,C),C);
	    log1 debug<<"Done initializing\n";
	    log1 debug<<"Traversing\n";
		symmetricTraverse(traverser,target,&s);
	    log1 debug<<"Done traversing\n";

		target.getFanRef().printWithIndices(&p,
											FPF_default/*|
											(symmetryOption.getValue()?FPF_group|FPF_conesCompressed:0)*/,
											&s);

    return 0;
/*    if(symmetryOption.getValue())
      {
	IntegerVectorList generators=FileParser(Stdin).parseIntegerVectorList();
	for(IntegerVectorList::const_iterator i=generators.begin();i!=generators.end();i++)
	  {
	    assert(i->size()==n);
	    FieldMatrix M1=integerMatrixToFieldMatrix(A,Q);
	    FieldMatrix M2=integerMatrixToFieldMatrix(rowsToIntegerMatrix(SymmetryGroup::permuteIntegerVectorList(A.getRows(),*i)),Q);
	    M1.reduce();
	    M1.REformToRREform(true);
	    M2.reduce();
	    M2.REformToRREform(true);

	    if(!(M1==M2))
	      {
		AsciiPrinter(Stderr) << "Permutation "<< *i <<
		  " does not keep the configuration fixed.\n";
		assert(0);
	      }
	  }
	s.computeClosure(generators);
      }
*/
    /* If the vector configuration A does not have full rank then
       change coordinates. */
/*    if(rank(A)!=A.getHeight())
      {
	FieldMatrix M=integerMatrixToFieldMatrix(A,Q);
	M.reduce(false,true);//force integer operations - preserving volume
	M.removeZeroRows();
	A=fieldMatrixToIntegerMatrix(M);
      }


    Triangulation2 t;
    t.A=A;

    // Convert a Triangulation to a Triangulation2
    {
      list<Triangulation::Cone> T=Triangulation::triangulate(A.transposed());
      for(list<Triangulation::Cone>::const_iterator i=T.begin();i!=T.end();i++)
	{
	  IntegerVector v=i->size();
	  int J=0;
	  for(Triangulation::Cone::const_iterator j=i->begin();j!=i->end();j++,J++)
	    v[J]=*j;
	    t.bases.insert(v);
	}
    }

{
    	AsciiPrinter p(Stdout);
		SymmetricTargetFanBuilder target(n,s);

		if(!optionRestrictingFan.getValue())
		{
			SecondaryFanTraverser traverser(t);
			symmetricTraverse(traverser,target,&s);
		}
		else
		{
		    PolyhedralFan f1=PolyhedralFan::readFan(optionRestrictingFan.getValue(),true,0,0,0,true);

		    for(PolyhedralFan::coneIterator i=f1.conesBegin();i!=f1.conesEnd();i++)
			{
		    	static int a;
		    	log2 cerr<<"Processing Cone "<<a++<<" which has dimension "<<i->dimension()<<endl;
		    	SecondaryFanTraverser traverser(triangulationWithFullDimensionalIntersection(t,*i),*i);
				symmetricTraverse(traverser,target,&s);
			}
		}

		target.getFanRef().printWithIndices(&p,
											FPF_default|
											(symmetryOption.getValue()?FPF_group|FPF_conesCompressed:0),
											&s);

		}*/
    return 0;
  }
};

static FiberPolytopeApplication theApplication;

