#ifndef EP_XFIG_H_INCLUDED
#define EP_XFIG_H_INCLUDED

#include "enumeration.h"

#include <list>

class XfigEnumerationPrinter: public EnumerationFilePrinter
{
  int basisCounter;
  bool largerTriangle;
  int variableShift;
  class XFig *xfig;
public:
  XfigEnumerationPrinter(bool _largerTriangle=false);
  bool basis(const PolynomialSet &groebnerBasis);
  string extension();
  void beginEnumeration(const PolynomialSet &groebnerBasis);
  void endEnumeration();
  void onOpened();
  void onClose();
  void onClosed();
  void setVariableShift(int shift);
};

#endif
