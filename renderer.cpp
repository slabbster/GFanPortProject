#include "renderer.h"

#include <math.h>

int StandardMonomialRenderer::getOffsetX()
{
  return boxSize*maxEntry*(1+2*(position%numberOfDrawingsPerLine));
}


int StandardMonomialRenderer::getOffsetY()
{
  return boxSize*maxEntry*(1+2*(position/numberOfDrawingsPerLine));
}


bool StandardMonomialRenderer::isInInitialIdeal(const PolynomialSet &s,int x, int y, int z)
{
  IntegerVector v(3);
  v[0]=x;
  v[1]=y;
  v[2]=z;
  PolynomialSet::const_iterator i;
  for(i=s.begin();i!=s.end();i++)
    if(i->getMarked().m.exponent.divides(v))break;

  return i!=s.end();
}


int StandardMonomialRenderer::trace(const PolynomialSet &s, int x, int y, int z, int d)
{
  int pos=0;

  while(1)
    {
      bool wasIn=isInInitialIdeal(s,x,y,z);
      if(!wasIn)return -1;

      switch(d)
	{
	case 0:
	  pos=x;
	  x--;
	  break;
	case 1:
	  pos=y;
	  y--;
	  break;
	case 2:
	  pos=z;
	  z--;
	  break;
	}

      if(x<0 || y<0 || z<0)break;

      bool isIn=isInInitialIdeal(s,x,y,z);

      if(!isIn && wasIn)
	return pos;
    }

  return 0;
}


void StandardMonomialRenderer::putPoint(float x, float y, int rotation, float size)
{ //rotation is clockwise
  fprintf(f," %i %i",getOffsetX()+int(x+size*sin(rotation*3.1415927f/3)),getOffsetY()+int(y-size*cos(rotation*3.1415927f/3)));
}


void StandardMonomialRenderer::drawQuad(int x, int y, int z, int rotation, int color, int intensity)
{
  fprintf(f,"2 3 0 1 0 %i 50 0 %i 0.000 0 0 -1 0 0 5\n        ",color, intensity);
  putPoint((x-y)*boxSize*sqrt(0.75),boxSize*(0.5*(x+y)-z),0,0);
  putPoint((x-y)*boxSize*sqrt(0.75),boxSize*(0.5*(x+y)-z),rotation,boxSize);
  putPoint((x-y)*boxSize*sqrt(0.75),boxSize*(0.5*(x+y)-z),rotation+1,boxSize);
  putPoint((x-y)*boxSize*sqrt(0.75),boxSize*(0.5*(x+y)-z),rotation+2,boxSize);
  putPoint((x-y)*boxSize*sqrt(0.75),boxSize*(0.5*(x+y)-z),0,0);
  fprintf(f,"\n");
}


StandardMonomialRenderer::StandardMonomialRenderer(FILE *f):
  numberOfDrawingsPerLine(5),
  maxEntry(8),
  position(0),
  boxSize(120)
{
  this->f=f;
  fprintf(f,"#FIG 3.2\n"
          "Landscape\n"
          "Center\n"
          "Inches\n"
          "Letter\n"
          "100.00\n"
          "Single\n"
          "-2\n"
          "1200 2\n");
}


void StandardMonomialRenderer::setMaxEntry(int maxEntry)
{
  this->maxEntry=maxEntry;
}


void StandardMonomialRenderer::setBoxSize(int boxSize)
{
  this->boxSize=boxSize;
}


void StandardMonomialRenderer::setNumberOfDrawingsPerLine(int numberOfDrawingsPerLine)
{
  this->numberOfDrawingsPerLine=numberOfDrawingsPerLine;
}


void StandardMonomialRenderer::render(const PolynomialSet &s)
{
  int x,y,z;
  bool color=false;

  x=maxEntry;
  for(y=0;y<maxEntry;y++)
    for(z=0;z<maxEntry;z++)
      {
	int pos=trace(s,x,y,z,0);
	if(pos>=0)
	  {
	    if(color)
	      drawQuad(pos,y,z,4,pos?15:7);
	    else
	      drawQuad(pos,y,z,4,7,pos?10:20);//dark
	  }
      }

  y=maxEntry;
  for(x=0;x<maxEntry;x++)
    for(z=0;z<maxEntry;z++)
      {
	int pos=trace(s,x,y,z,1);
	if(pos>=0)
	  {
	    if(color)
	      drawQuad(x,pos,z,0,pos?17:7);
	    else
	      drawQuad(x,pos,z,0,7,pos?12:20);
	  }
      }

  z=maxEntry;
  for(x=0;x<maxEntry;x++)
    for(y=0;y<maxEntry;y++)
      {
	int pos=trace(s,x,y,z,2);
	if(pos>=0)
	  {
	    if(color)
	      drawQuad(x,y,pos,2,pos?3:7);
	    else
	      drawQuad(x,y,pos,2,7,pos?14:20);//light
     	  }
      }

  position++;
}
