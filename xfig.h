#ifndef XFIG_H_INCLUDED
#define XFIG_H_INCLUDED

#include <stdio.h>
#include "vektor.h"
#include <string>

// In time this class should also be used in renderer.cpp and in ep_xfig.cpp

class XFig
{
  FILE *f;
  IntegerVectorList p;
  bool arrowOrigin;
  bool arrowTarget;
  int offsetX,offsetY;
 public:
  XFig(FILE *f_):
    f(f_),
    offsetX(0),
    offsetY(0)
    {
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

  void setOffset(int x, int y)
    {
      offsetX=x;
      offsetY=y;
    }
  void beginDrawLine(bool arrowOrigin_=false, bool arrowTarget_=false, int thickness=1)
    {
      arrowOrigin=arrowOrigin_;
      arrowTarget=arrowTarget_;
      fprintf(f,"2 1 0 %i 0 7 50 -1 -1 0.000 0 0 -1 %i %i",thickness,int(arrowTarget),int(arrowOrigin));
    }
  void addPoint(int x, int y)
    {
      IntegerVector v(2);
      v[0]=x+offsetX;
      v[1]=y+offsetY;
      p.push_back(v);
    }
  void endDrawLine()
    {
      assert(p.size()>0);
      fprintf(f," %i\n",p.size());
      if(arrowOrigin)
	fprintf(f,"\t 2 1 1.00 60.00 120.00\n");
      if(arrowTarget)
	fprintf(f,"\t 2 1 1.00 60.00 120.00\n");
      for(IntegerVectorList::const_iterator i=p.begin();i!=p.end();i++)
	{
	  fprintf(f," %i %i",(*i)[0]+offsetX,(*i)[1]+offsetY);
	}
      p=IntegerVectorList();
      fprintf(f,"\n");
    }
  void drawString(int x, int y, string const &s, int size=12)
    {
      fprintf(f,"4 0 0 48 -1 0 %i 0.0000 4 135 270 %i %i ",size,x+offsetX,y+offsetY);
      fprintf(f,"%s",s.c_str());
      fprintf(f,"\\001\n");
    }   
  void drawDot(int x, int y, int color=0, int radius=75)
    {
      fprintf(f,"1 3 0 1 %i %i 49 -1 20 0.000 1 0.0000 %i %i %i %i %i %i %i %i\n",0,color,x+offsetX,y+offsetY,radius,radius,x+offsetX,y+offsetY,x+offsetX+radius,y+offsetY);
    }

  // code for drawing intersection of halfspaces 
  struct Point
  {
    float x,y,z;
    float dot(const Point &p)const{return x*p.x+y*p.y+z*p.z;}
    bool isInside(const Point &p)const{return dot(p)>-0.05;}
    //  bool isInside(const Point &p)const{return dot(p)>0;}
    Point(float x_, float y_, float z_):x(x_),y(y_),z(z_){}
    void print()const;
    Point intersect(Point a, Point b)const
    {
      float A=dot(a);
      float B=dot(b);
      float cA=B/(B-A);
      float cB=-A/(B-A);

      return Point(a.x*cA+b.x*cB,a.y*cA+b.y*cB,a.z*cA+b.z*cB);
    }
  };
  typedef list<Point> Polygon;
  void kickPoint(const Point &p, int mode);
  void printPolygon(const Polygon &p);
  void drawPolygon(const Polygon &vertices,int mode);
  Polygon intersect(const Polygon &polygon, const Point &normal);
};


#endif
