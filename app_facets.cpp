#include "parser.h"
#include "printer.h"
#include "wallideal.h"
#include "lp.h"
#include "gfanapplication.h"
#include "log.h"

class FacetsApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program computes the facet normals of a full-dimensional cone specified on the input. The input is a list of vectors each defining a half space. The cone is the intersection of these half spaces. The output is a minimal list of normals defining the cone.\n";
  }
  FacetsApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_facets";
  }
  int main()
  {
    FileParser P(Stdin);

    IntegerVectorList normals=P.parseIntegerVectorList();

    normals=wallRemoveScaledInequalities(normals);
    IntegerVectorList facets;
    for(IntegerVectorList::const_iterator i=normals.begin();i!=normals.end();i++)
      if(isFacet(normals,i))
        {
    	  facets.push_back(*i);
        }

    AsciiPrinter(Stdout).printVectorList(facets);

    log1 fprintf(Stderr,"Rank:%i\n",rankOfMatrix(facets));

    return 0;
  }
};

static FacetsApplication theApplication;
