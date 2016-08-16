#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"

class SupportIndicesApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  SupportIndicesApplication()
  {
    registerOptions();
  }
  const char *name()
  {
    return "_supportindices";
  }
  int main()
  {
    FileParser P(Stdin);
    IntegerVector v=P.parseIntegerVector();

    for(int i=0;i<v.size();i++)
      if(v[i]!=0)printf("%i\n",i);

    return 0;
  }
  const char *helpText()
  {
    return "Takes a vektor and prints the indices for which the vector is nonzero.\n"
;
  }
};

static SupportIndicesApplication theApplication;
