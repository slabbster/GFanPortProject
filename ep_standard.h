#ifndef EP_STANDARD_H_INCLUDED
#define EP_STANDARD_H_INCLUDED

#include "enumeration.h"

class StandardEnumerationPrinter: public EnumerationFilePrinter
{
  int basisCounter;
  bool insertComma;
public:
  StandardEnumerationPrinter();
  bool basis(const PolynomialSet &groebnerBasis);
  string extension();
  void beginEnumeration(const PolynomialSet &groebnerBasis);
  void endEnumeration();
  void onOpened();
  void onClose();
  void onClosed();
};

#endif
