#include "monomial.h"

#include "printer.h"
#include <sstream>

Monomial::Monomial(PolynomialRing const &r,const IntegerVector &v):exponent(v),theRing(r)
{
  if(v.size()!=r.getNumberOfVariables())
    {
      AsciiPrinter(Stderr).printPolynomialRing(r);
      AsciiPrinter(Stderr).printVector(v); 
      assert(v.size()==r.getNumberOfVariables());
    }
}


string Monomial::toString(bool alwaysWriteSign, bool writeIfOne, bool latex/*, bool mathMode*/)const
{
  stringstream s;

  /*  if(latex & !mathMode)
    s << "$";
  */
  const int sign=1;

  bool variablePrinted=false;
  for(int i=0;i<exponent.size();i++)if(exponent[i]*sign>0)
    {
      s << getRing().getVariableName(i);
      if(int(exponent[i]*sign)!=1)
	{
	  s << "^";
	  if(latex)
	    s << "{";
	  s << int(exponent[i]*sign);
	  if(latex)
	    s << "}";
	}
      variablePrinted=true;
    }
  if(!variablePrinted && writeIfOne)
    s<< "1";

  /*  if(latex & !mathMode)
    s << "$";
  */
  return s.str();
}
