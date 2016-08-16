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
#include "matrix.h"
#include "latticeideal.h"
#include "subspace.h"
#include "scarf.h"

class ScarfComplexApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program .\n";
  }
  ScarfComplexApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_scarf_complex";
  }
  int main()
  {
    LpSolver::printList(Stdout);
    lpSetSolver("cddgmp");

    FileParser P(Stdin);

    IntegerVectorList ivl=P.parseIntegerVectorList();
    IntegerMatrix A=rowsToIntegerMatrix(ivl);

    bool s1=satisfiesA1(A);
    fprintf(Stdout,"A1 satisfied:%i\n",s1);
    if(s1)
      {
	IntegerVectorList N=neighbours(A);
	AsciiPrinter Q(Stdout);
	Q.printVectorList(N);
	fprintf(Stdout,"A2 satisfied:%i\n",satisfiesA2(A));
	fprintf(Stdout,"A3 satisfied:%i\n",satisfiesA3(A,&N));

	N=orientedNeighbours(N,-A[0]);
	IntegerMatrix N2(1,A.getWidth());
	N2.append(rowsToIntegerMatrix(N));
	N2.append((-1)*rowsToIntegerMatrix(N));

	Q.printVectorList(N2.getRows());

	/*	simplex=StringParser("(0,5,6)").parseIntegerVector();
	Q.printVector(simplex);
	simplex=kFlip(A,N2,simplex,0);
	Q.printVector(simplex);
	*/
	IntegerVector simplex=computeMaximalScarfSimplex(A,N2);
	traverseScarfComplex(A,N2,simplex);
      }
    else
      fprintf(Stdout,"Neighbours for this matrix cannot be computed with the implemented method.\n");

    return 0;
  }
};

static ScarfComplexApplication theApplication;
