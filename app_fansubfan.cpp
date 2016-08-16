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

class SubfanApplication : public GFanApplication
{
  StringOption inputOption;
  SimpleOption optionSymmetry;
public:
  const char *helpText()
  {
    return "This program takes a polyhedral fan and a list of vectors and computes the smallest subfan of the fan having the list of vectors in its support.\n";
  }
  SubfanApplication():
    inputOption("-i","Specify the name of the input file.","polymake.out"),
    optionSymmetry("--symmetry",
                   "Reads in the cone stored with symmetry. The generators of the symmetry group must be given on the standard input.\n")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_fansubfan";
  }

  int main()
  {
    IntegerVectorList theVectors=FileParser(Stdin).parseIntegerVectorList();
    assert(theVectors.size());
    int n=theVectors.begin()->size();
    SymmetryGroup s(n);

    if(optionSymmetry.getValue())
      {
        IntegerVectorList generators=FileParser(Stdin).parseIntegerVectorList();
        s.computeClosure(generators);
      }
    PolyhedralFan f1=PolyhedralFan::readFan(inputOption.getValue(),true,0,0,(optionSymmetry.getValue())?&s:0);

    PolyhedralFan f2(n);

    for(IntegerVectorList::const_iterator i=theVectors.begin();i!=theVectors.end();i++)
      f2.insert(f1.coneContaining(*i));

    AsciiPrinter P(Stdout);
    f2.printWithIndices(&P,FPF_default | FPF_primitiveRays | FPF_multiplicities);

    return 0;
  }
};

static SubfanApplication theApplication;

