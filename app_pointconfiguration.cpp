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

class PointConfigurationApplication : public GFanApplication
{
  IntegerOption dimension1Option;
  IntegerOption dimension2Option;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program produces the point configuration of the product of two standard simplices together with their symmetries.\n";
  }
  PointConfigurationApplication():
    dimension1Option("-d1","Number of vertices in the first simplex.",2),
    dimension2Option("-d2","Number of vertices in the second simplex.",2)
  {
    registerOptions();
  }

  const char *name()
  {
    return "_pointconfiguration";
  }

  int main()
  {
    int d1=dimension1Option.getValue();
    int d2=dimension2Option.getValue();
    assert(d1>=2);
    assert(d2>=2);
    IntegerMatrix A(d1*d2,d1+d2);
    IntegerVector p1(d1*d2);
    IntegerVector p2(d1*d2);
    IntegerVector p3(d1*d2);
    IntegerVector p4(d1*d2);
    for(int b=0;b<d2;b++)
      for(int a=0;a<d1;a++)
        {
          A[b*d1+a][b]=1;
          A[b*d1+a][a+d2]=1;

          p1[b*d1+a]=b*d1+((a==0 || a==1)?1-a:a);
          p2[b*d1+a]=b*d1+((a+1)%d1);
          p3[b*d1+a]=((b==0 || b==1)?1-b:b)*d1+a;
          p4[b*d1+a]=((b+1)%d2)*d1+a;
        }
    pout<<A.getRows();

    IntegerVectorList p;
    p.push_back(p1);
    p.push_back(p2);
    p.push_back(p3);
    p.push_back(p4);
    pout<<p;

    return 0;
  }
};

static PointConfigurationApplication theApplication;
