#include "ep_xfig.h"

#include "printer.h"
#include "wallideal.h"
#include "lp.h"
#include "xfig.h"
#include "log.h"

XfigEnumerationPrinter::XfigEnumerationPrinter(bool _largerTriangle):
  largerTriangle(_largerTriangle),
  variableShift(0),
  xfig(0)
{
}


void XfigEnumerationPrinter::onOpened()
{
  xfig= new XFig(file);
}


void XfigEnumerationPrinter::onClose()
{
}


void XfigEnumerationPrinter::onClosed()
{
  delete xfig;
}


void XfigEnumerationPrinter::beginEnumeration(const PolynomialSet &groebnerBasis)
{
  basisCounter=0;
}


void XfigEnumerationPrinter::endEnumeration()
{
}



bool XfigEnumerationPrinter::basis(const PolynomialSet &groebnerBasis)
{
  assert(xfig);
  XFig::Polygon p;

  if(largerTriangle)
    {
      p.push_back(XFig::Point(3,-1,-1));
      p.push_back(XFig::Point(-1,3,-1));
      p.push_back(XFig::Point(-1,-1,3));
    }
  else
    {
      p.push_back(XFig::Point(1,0,0));
      p.push_back(XFig::Point(0,1,0));
      p.push_back(XFig::Point(0,0,1));
    }

  IntegerVectorList normals=wallInequalities(groebnerBasis);
  for(IntegerVectorList::const_iterator i=normals.begin();i!=normals.end();i++)
    //if(wallContainsPositiveVector(*i))
    if(isFacet(normals,i))
        {
          if(i->size()>=3)
            {
	      //  AsciiPrinter(Stderr).printVector(*i);
	      XFig::Point n((*i)[variableShift%i->size()],(*i)[(variableShift+1)%i->size()],(*i)[(variableShift+2)%i->size()]);
              p=xfig->intersect(p,n);
	      //  fprintf(Stderr,"%i\n",variableShift);
            }
        }

  log2 xfig->printPolygon(p);
  xfig->drawPolygon(p,0);

  basisCounter++;

  log2 fprintf(Stderr,"basisCounter:%i\n",basisCounter);

  return true;
}


string XfigEnumerationPrinter::extension()
{
  return "";
}

void XfigEnumerationPrinter::setVariableShift(int shift)
{
  variableShift=shift;
}
