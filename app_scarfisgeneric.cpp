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

class ScarfIsGenericApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program takes a matrix as input and checks if it satisfies Scarf's generality conditions. The rows of the matrix are listed on the input. The A1 condition is that there exists a strictly poistive vector in the co-kernel of the matrix. The A2 condition is that te maximal minors of the matrix are non-zero. A3\n";
  }
  ScarfIsGenericApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_scarf_isgeneric";
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
	N=orientedNeighbours(N,-A[0]);
	Q.printVectorList(N);
	fprintf(Stdout,"A2 satisfied:%i\n",satisfiesA2(A));
	for(int i=0;i<A.getHeight();i++)
	  fprintf(Stdout,"A3(%i) satisfied:%i\n",i,satisfiesA3i(A,i,&N));
	fprintf(Stdout,"A3 satisfied:%i\n",satisfiesA3(A,&N));
      }
    else
      fprintf(Stdout,"Neighbours for this matrix cannot be computed with the implemented method.\n");

    return 0;
  }
};

static ScarfIsGenericApplication theApplication;
