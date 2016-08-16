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

class FanConesApplication : public GFanApplication
{
  StringOption inputOption;
  SimpleOption resultantOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program lists the cones of a polyhedral fan.\n";
  }
  FanConesApplication():
    inputOption("-i","Specify the name of the input file.","polymake.out"),
    resultantOption("--resultant","Take codim 1 skeleton and wipe out bad cones.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_fancones";
  }

  bool zeroOrTwo(int v)
  {
    return (v==0) || (v==2);
  }
  int main()
  {
    PolyhedralFan f1=PolyhedralFan::readFan(inputOption.getValue(),true,0,0,0);
    AsciiPrinter P(Stdout);

    if(resultantOption.getValue())
      {
        PolyhedralFan f2=f1.facetComplex();
        PolyhedralFan f3(f2.getAmbientDimension());
        for(PolyhedralFan::coneIterator i=f2.conesBegin();i!=f2.conesEnd();i++)
          {
            IntegerVector v=i->getEquations().front();
            IntegerVector u=v.supportAsZeroOneVector();

            if(zeroOrTwo(u[0]+u[1]+u[2]+u[3]) &&
                zeroOrTwo(u[4]+u[5]+u[6]) &&
                zeroOrTwo(u[7]+u[8]/*+u[5]*/))
              f3.insert(*i);
          }
        f3.printWithIndices(&pout);
        //        pout << f2;
      }
    else
    for(PolyhedralFan::coneIterator i=f1.conesBegin();i!=f1.conesEnd();i++)
      {
        P<<*i;
        P<<"-------------------------------------\n";
      }



    return 0;
  }
};

static FanConesApplication theApplication;

