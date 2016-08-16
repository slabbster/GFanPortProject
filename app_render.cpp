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
#include "renderer.h"
#include "xfig.h"
#include <math.h>

class RenderApplication : public GFanApplication
{
  SimpleOption optionLargerTriangle;
  IntegerOption optionShiftVariablesWhenDrawing;
public:
  const char *helpText()
  {
    return "This program renders a Groebner fan as an xfig file. To be more precise, the input is the list of all reduced Groebner bases of an ideal. The output is a drawing of the Groebner fan intersected with a triangle. The corners of the triangle are (1,0,0) to the right, (0,1,0) to the left and (0,0,1) at the top. If there are more than three variables in the ring these coordinates are extended with zeros. It is possible to shift the 1 entry cyclic with the option --shiftVariables.\n";
  }
  RenderApplication():
    optionLargerTriangle("-L",
			 "Make the triangle larger so that the shape of the Groebner region appears."),
    optionShiftVariablesWhenDrawing("--shiftVariables",
				    "Shift the positions of the variables in the drawing. For example with the value equal to 1 the corners will be right: (0,1,0,0,...), left: (0,0,1,0,...) and top: (0,0,0,1,...). The shifting is done modulo the number of variables in the polynomial ring. The default value is 0.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_render";
  }

  int main()
  {
    FileParser P(Stdin);

    PolynomialSetList l;
    l=P.parsePolynomialSetListWithRing();

    XfigEnumerationPrinter ep(optionLargerTriangle.getValue());
    ep.setVariableShift(optionShiftVariablesWhenDrawing.getValue());

    ep.open(Stdout);

    if(l.size())
      {
	ep.beginEnumeration(*l.begin());
	for(PolynomialSetList::const_iterator i=l.begin();i!=l.end();i++)ep.basis(*i);
	ep.endEnumeration();
      }
    ep.close();

    return 0;
  }
};

static RenderApplication theApplication;
