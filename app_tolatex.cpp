#include "parser.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"


class ToLatexApplication : public GFanApplication
{
  SimpleOption optionAddHeader;
  SimpleOption optionPolynomialSet;
  SimpleOption optionPolynomialSetList;
public:
  const char *helpText()
  {
    return "This program converts ASCII math to TeX math. The data-type is specified by the options.\n";
  }
  ToLatexApplication():
    optionPolynomialSet("--polynomialset_","The data to be converted is a list of polynomials."),
    optionPolynomialSetList("--polynomialsetlist_","The data to be converted is a list of lists of polynomials."),
    optionAddHeader("-h","Add a header to the output. Using this option the output will be LaTeXable right away.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_tolatex";
  }

  int main()
  {
    LatexPrinter out(Stdout);
    FileParser in(Stdin);


    if(optionAddHeader.getValue())
      fprintf(Stdout,"\\documentclass{article}\n"
	    "\\usepackage[english]{babel}\n"
	    "\\begin{document}\n");

    if(optionPolynomialSet.getValue())
      {
	out.printPolynomialSet(in.parsePolynomialSetWithRing());
      }
    else if(optionPolynomialSetList.getValue())
      {
	out.printPolynomialSetList(in.parsePolynomialSetListWithRing());
      }
    else

      {
	fprintf(Stderr,"No type specified\n");
	assert(0);
      }

    if(optionAddHeader.getValue())
      fprintf(Stdout,"\\end{document}\n");
    return 0;
  }
};

static ToLatexApplication theApplication;
