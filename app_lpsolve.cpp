#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "fieldlp.h"

class LPSolveApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  LPSolveApplication()
  {
    registerOptions();
  }
  const char *name()
  {
    return "_lpsolve";
  }
  int main()
  {
    FileParser P(Stdin);
    FieldMatrix A=integerMatrixToFieldMatrix(rowsToIntegerMatrix(P.parseIntegerVectorList()),Q);
    IntegerVector b0=P.parseIntegerVector();
    FieldVector b=integerVectorToFieldVector(b0,Q);

    FieldLP lp(A,b);
    FieldVector w=integerVectorToFieldVector(P.parseIntegerVector(),Q);
    lp.setObjectiveFunction(w);

    /*    set<int> B;
    for(int i=0;i<3;i++)B.insert(2*i+1);
    lp.setCurrentBasis(B);
    */

    if(!lp.findFeasibleBasis())
      {
	fprintf(Stderr,"LP is infeasible\n");
      }
    else
      {
	int status;
	do
	  {
	    AsciiPrinter Q(Stdout);
	    lp.print(Q);
	    status=lp.step();
	  }
	while(status==1);
	fprintf(Stderr,status?"LP is unbounded.\n":"Optimal solution found.\n");
      }

    return 0;
  }
  const char *helpText()
  {
    return "Test program for FieldLP.\n";
  }
};

static LPSolveApplication theApplication;
