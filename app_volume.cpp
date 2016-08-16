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
#include "polyhedralfan.h"
#include "symmetry.h"

class VolumeApplication : public GFanApplication
{
  StringOption inputOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program computes the volume of a polyhedral fan intersected with the standard cube with edge length 2 centered at the origin. This is done for each dimension\n";
  }
  VolumeApplication():
    inputOption("-i","Specify the name of the input file.","polymake.out")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_volume";
  }

  int main()
  {
    PolyhedralFan F=PolyhedralFan::readFan(inputOption.getValue(),true,0,0,0,true);

    IntegerVectorList generators=FileParser(Stdin).parseIntegerVectorList();

    SymmetryGroup s(F.getAmbientDimension());
    s.computeClosure(generators);

    cout << F.volume(21,&s).toString();

    return 0;
  }
};

static VolumeApplication theApplication;
