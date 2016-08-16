/*
  TO DO:
  Find a good way to input the name of the new variable and adjust this to helpText.
*/
#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "termorder.h"
#include "gfanapplication.h"
#include "wallideal.h"

class HomogenizeApplication : public GFanApplication
{
  SimpleOption optionIdeal;
  SimpleOption optionHomogenizationVector;
  SimpleOption optionH;
public:
  const char *helpText()
  {
    return "This program homogenises a list of polynomials by introducing an extra variable. The name of the variable to be introduced is read from the input after the list of polynomials. Without the -w option the homogenisation is done with respect to total degree.\n"
      "Example:\n"
      "Input:\n"
      "Q[x,y]{y-1}\n"
    "z\n"
    "Output:\n"
      "Q[x,y,z]{y-z}\n";
;
  }
  HomogenizeApplication():
    optionIdeal("-i","Treat input as an ideal. This will make the program compute the homogenisation of the input ideal. This is done by computing a degree Groebner basis and homogenising it."),
    optionHomogenizationVector("-w","Specify a homogenisation vector. The length of the vector must be the same as the number of variables in the ring. The vector is read from the input after the list of polynomials.\n"),
    optionH("-H","Let the name of the new variable be H rather than reading in a name from the input.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_homogenize";
  }

  int main()
  {
    FileParser P(Stdin);

    PolynomialSet g=P.parsePolynomialSetWithRing();

    string homogenizingVariableName=(optionH.getValue())?"H":P.parseVariableName();

    IntegerVector grading;
    if(optionHomogenizationVector.getValue())
      grading=P.parseIntegerVector();
    else
      grading=IntegerVector::allOnes(g.numberOfVariablesInRing());

    if(optionIdeal.getValue())
      {
	WeightReverseLexicographicTermOrder t(grading);
	buchberger(&g,t);
      }

    PolynomialSet h(g.getRing().withVariablesAppended(homogenizingVariableName));

    AsciiPrinter(Stdout).printPolynomialRing(h.getRing());
    AsciiPrinter(Stdout).printNewLine();

    for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
      {
	h.push_back(i->homogenization(h.getRing(),&grading));
      }

    AsciiPrinter(Stdout).printPolynomialSet(h);

    return 0;
  }
};

static HomogenizeApplication theApplication;
