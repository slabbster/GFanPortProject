#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "lp.h"
#include "gfanapplication.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "symmetry.h"

#include "polymakefile.h"

class IsConnectedApplication : public GFanApplication
{
  StringOption inputOption;
  SimpleOption optionSymmetry;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program checks if a PURE d-dimensional fan is connected in codimension one (in the sense of [Bogart, Jensen, Speyer, Sturmfels, Thomas]).\n";
  }
  IsConnectedApplication():
    inputOption("-i","Specify the name of the input file.","polymake.out"),
    optionSymmetry("--symmetry",
		   "Reads in the cone stored with symmetry. The generators of the symmetry group are given on the input. IF SYMMETRY IS SPECIFIED THEN "
		   "THE CHECK IS NOT COMPLETE AND IT WILL ONLY BE CHECKED IF THAT GIVEN"
		   "TWO CONES THAT EXIST ELEMENTS IN THE TWO RESPECTIVE ORBITS"
		   "WHICH ARE CONNECTED BY A RIDGE PATH\n")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_isconnected";
  }

  int main()
  {
    PolyhedralFan f=PolyhedralFan::readFan(inputOption.getValue(),true,0,0,0,true);
    int n=f.getAmbientDimension();
    SymmetryGroup s(n);

    if(optionSymmetry.getValue())
      {
	IntegerVectorList generators=FileParser(Stdin).parseIntegerVectorList();
	s.computeClosure(generators);
      }

    cout << f.isConnected(&s) << endl;

    return 0;
  }
};

static IsConnectedApplication theApplication;

