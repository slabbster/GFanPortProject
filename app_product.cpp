#include "parser.h"
#include "printer.h"
#include "lp.h"
#include "gfanapplication.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "polymakefile.h"

class ProductApplication : public GFanApplication
{
  StringOption input1Option;
  StringOption input2Option;
public:
  const char *helpText()
  {
    return "This program takes two polyhedral fans and computes their product.\n";
  }
  ProductApplication():
    input1Option("-i1","Specify the name of the first input file.","polymake.out"),
    input2Option("-i2","Specify the name of the second input file.","polymake.out")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_fanproduct";
  }

  int main()
  {
    PolyhedralFan f1=PolyhedralFan::readFan(input1Option.getValue());
    PolyhedralFan f2=PolyhedralFan::readFan(input2Option.getValue());

    PolyhedralFan f=product(f1,f2);

    AsciiPrinter P(Stdout);

    f.printWithIndices(&P,FPF_default|FPF_multiplicities/*|FPF_values*/);

    return 0;
  }
};

static ProductApplication theApplication;

