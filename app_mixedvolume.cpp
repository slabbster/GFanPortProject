#include <iostream>
#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "mixedvolume.h"
#include "gfanapplication.h"

class MixedVolumeApplication : public GFanApplication
{
public:
  const char *helpText()
  {
    return "This program computes the mixed volume of the Newton polytopes of a list of polynomials. The ring is specified on the input. After this follows the list of polynomials.\n";
  }
  MixedVolumeApplication()
  {
    registerOptions();
  }
  const char *name()
  {
    return "_mixedvolume";
  }
  int main()
  {
    FileParser P(Stdin);
    PolynomialSet s=P.parsePolynomialSetWithRing();
    cout << mixedVolume(s) << endl;
    return 0;
  }
};

static MixedVolumeApplication theApplication;
