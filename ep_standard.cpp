#include "ep_standard.h"
#include "printer.h"


StandardEnumerationPrinter::StandardEnumerationPrinter():
  insertComma(false)
{
}


void StandardEnumerationPrinter::onOpened()
{
  fprintf(file,"{");
}


void StandardEnumerationPrinter::onClose()
{
  fprintf(file,"}");
}


void StandardEnumerationPrinter::onClosed()
{
}


void StandardEnumerationPrinter::beginEnumeration(const PolynomialSet &groebnerBasis)
{
  basisCounter=0;
}


void StandardEnumerationPrinter::endEnumeration()
{
}


bool StandardEnumerationPrinter::basis(const PolynomialSet &groebnerBasis)
{
  if(insertComma)fprintf(file,",\n");
  AsciiPrinter(file).printPolynomialSet(groebnerBasis);
  insertComma=true;

  basisCounter++;

  return true;
}


string StandardEnumerationPrinter::extension()
{
  return "";
}
