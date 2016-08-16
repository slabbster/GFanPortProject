#include "xfig.h"
#include "printer.h"

void XFig::Point::print()const
{
  fprintf(Stderr,"(%f,%f,%f)\n",x,y,z);
}

void XFig::printPolygon(const Polygon &p)
{
  fprintf(Stderr,"{\n");
  for(Polygon::const_iterator i=p.begin();i!=p.end();i++)
    i->print();
  fprintf(Stderr,"}\n");
}


XFig::Polygon XFig::intersect(const Polygon &polygon, const Point &normal)
{
  //  fprintf(Stderr,"intersect input:\n");
  //  printPolygon(polygon);

  Polygon ret;
  Polygon temp;

  Polygon::const_iterator i=polygon.begin();
  while(i!=polygon.end())
    {
      if(normal.isInside(*i))break;
      i++;
    }

  if(i==polygon.end())return ret;
  
  for(Polygon::const_iterator j=i;j!=polygon.end();j++)
    temp.push_back(*j);
  for(Polygon::const_iterator j=polygon.begin();j!=i;j++)
    temp.push_back(*j);
  temp.push_back(*i);

//  fprintf(Stderr,"intersect temp:\n");
//  printPolygon(temp);

  Polygon::const_iterator j=temp.begin();


  while(j!=temp.end())
    {
      Point last(0,0,0);
      while(j!=temp.end()&&normal.isInside(*j))
        {
          //fprintf(Stderr,"inside:");
          //j->print();
          last=*j;
          ret.push_back(*j);
          j++;
        }
      if(j==temp.end())break;
      Point leaving=normal.intersect(last,*j);
      while(!normal.isInside(*j))
        {
          //fprintf(Stderr,"outside:");
          //j->print();
          last=*j;
          j++;
        }        
      Point entering=normal.intersect(last,*j);
      ret.push_back(leaving);
      ret.push_back(entering);
      ret.push_back(*j);
      j++;
    }
  ret.pop_back();

  //  fprintf(Stderr,"intersect ret:\n");
  //  printPolygon(ret);

  return ret;
}


void XFig::kickPoint(const Point &p, int mode)
{
  //  fprintf(file," %i %i",(int)(p.x*10000),(int)(p.y*10000));
  switch(mode)
    {
    case 0:
      {
	float invSqrt2=0.707106781;
	float invSqrt6=0.40824829;
	fprintf(f," %i %i",(int)((invSqrt2+invSqrt2*p.x-invSqrt2*p.y)*10000+offsetX),(int)((2*invSqrt6+invSqrt6*p.x+invSqrt6*p.y-2*invSqrt6*p.z)*10000+offsetY));
      }
      break;
    case 1:
      {
	fprintf(f," %i %i",(int)((p.x)*10000+offsetX),(int)((p.y)*10000+offsetY));
      }
    }
}


void XFig::drawPolygon(const Polygon &vertices, int mode)
{
  if(vertices.size())
    {
      fprintf(f,"2 3 0 1 0 %i 50 0 25 0.000 0 0 -1 0 0 %i\n        ",7,vertices.size()+1);
      
      for(Polygon::const_iterator i=vertices.begin();i!=vertices.end();i++)
        kickPoint(*i,mode);
      kickPoint(*vertices.begin(),mode);
      fprintf(f,"\n");
    }
  fflush(f);
}

