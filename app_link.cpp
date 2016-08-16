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

class LinkApplication : public GFanApplication
{
  //  FieldOption theFieldOption;
  StringOption inputOption;
  SimpleOption optionSymmetry;
  SimpleOption optionConeIndices;
  SimpleOption optionNonMaximal;
  SimpleOption optionStar;
public:
  const char *helpText()
  {
    return "This program takes a polyhedral fan and a vector and computes the link of the polyhedral fan around that vertex. The link will have lineality space dimension equal to the dimension of the relative open polyhedral cone of the original fan containing the vector.\n";
  }
  LinkApplication():
    inputOption("-i","Specify the name of the input file.","polymake.out"),
    optionSymmetry("--symmetry",
		   "Reads in a fan stored with symmetry. The generators of the symmetry group must be given on the standard input.\n"),
    optionConeIndices("--coneIndices","Restrict the fan to a set of cones. The indices of the cones are given on the input."),
    optionNonMaximal("--nonMaximal","Tells the program that the cone indices are with respect to the non maximal cone list."),
    optionStar("--star","Computes the star instead. The star is defined as the smallest polyhedral fan containing all cones of the original fan containing the vector.")
  {
    registerOptions();
    optionConeIndices.hide();
    optionNonMaximal.hide();
  }

  const char *name()
  {
    return "_fanlink";
  }

  int main()
  {
    IntegerVector v=FileParser(Stdin).parseIntegerVector();
    int n=v.size();
    SymmetryGroup s(n);



    if(optionSymmetry.getValue())
      {
	IntegerVectorList generators=FileParser(Stdin).parseIntegerVectorList();
	s.computeClosure(generators);
      }
    PolyhedralFan f(n);
    if(optionConeIndices.getValue())
      {
	IntegerVector temp=FileParser(Stdin).parseIntegerVector();
	set<int> coneIndices;
	for(int i=0;i<temp.size();i++)coneIndices.insert(temp[i]);
	f=PolyhedralFan::readFan(inputOption.getValue(),!optionNonMaximal.getValue(),&v,&coneIndices,(optionSymmetry.getValue())?&s:0);
      }
    else
      f=PolyhedralFan::readFan(inputOption.getValue(),true,&v,0,(optionSymmetry.getValue())?&s:0);

    if(!optionStar.getValue())f=f.link(v);
    AsciiPrinter P(Stdout);
    f.printWithIndices(&P,FPF_default | FPF_primitiveRays | FPF_multiplicities);
    //    f.printWithIndices(&P,false,0,0);
    //f.printWithIndices(&P,false,&s,true);

    return 0;
  }
};

static LinkApplication theApplication;

