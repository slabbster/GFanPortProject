#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "termorder.h"
#include "ep_standard.h"
#include "ep_xfig.h"
#include "gfanapplication.h"
#include "matrix.h"
#include "latticeideal.h"
#include "subspace.h"
#include "scarf.h"
#include "xfig.h"

class ScarfVisualizeApplication : public GFanApplication
{
public:
  IntegerOption optionNumberOfSteps;
  int radius;
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program program produces a projective drawing of a slice of the Scarf fan. The Scarf fan is not a fan in the strict polyhedral sense. The input for the program is either a\n";
  }
  ScarfVisualizeApplication():
    optionNumberOfSteps("-n","Number of steps",10)
  {
    registerOptions();
  }

  const char *name()
  {
    return "_scarf_visualize";
  }

  int colorMap(bool A1, bool A2, bool A3)
  {
    if(A3)
      {
	return 7;
      }
    else
      {
	if(A2)
	  {
	    if(A1)
	      return 4;
	    else
	      return 3;
	  }
	else
	  {
	    if(A1)
	      return 1;
	    else
	      return 0;
	  }
      }
    return 0;
  }

  void box(XFig &xfig)
  {
    xfig.beginDrawLine(false,false,4);
    xfig.addPoint(-10000,-10000);
    xfig.addPoint(10000,-10000);
    xfig.addPoint(10000,10000);
    xfig.addPoint(-10000,10000);
    xfig.addPoint(-10000,-10000);
    xfig.endDrawLine();
  }
  void printVector(XFig &xfig,int x, int y, IntegerVector const &v)
  {
    if(v.size()==3)
      {
	char s[50];
	sprintf(s,"%i %i %i",v[0],v[1],v[2]);
	xfig.drawString(x,y,string(s),60);
      }
    else
    if(v.size()==2)
      {
	char s[50];
	sprintf(s,"%i  %i",v[0],v[1]);
	xfig.drawString(x,y,string(s),60);
      }
    else
      {
	assert(0);
      }
  }
  void renderFace(IntegerMatrix A, IntegerMatrix const &mRot2, XFig &xfig,int steps)
  {
    IntegerVector v(3);

    v[2]=steps;
    for(int x=-steps;x<=steps;x++)
      {
	v[0]=x;
	for(int y=-steps;y<=steps;y++)
	  {
	    v[1]=y;

	    A[0]=mRot2.vectormultiply(v);

	    if(x==-steps && y==-steps)
	      printVector(xfig,-10000,-10000+800,A[0]);
	    if(x==0 && y==-steps)
	      printVector(xfig,0,-10000+800,A[0]);
	    if(x==-steps && y==0)
	      printVector(xfig,-10000,800,A[0]);
	    if(x==0 && y==0)
	      printVector(xfig,0,0+800,A[0]);


	    //	    AsciiPrinter(Stderr).printVectorList(A.getRows());

	    bool A3=false;
	    bool A1=satisfiesA1(A);
	    bool A2=satisfiesA2(A);
	    if(A1 && A2)
	      {
		IntegerVectorList N=neighbours(A);
		A3=satisfiesA3(A,&N);

		if(A1 && A2 && A3)
		  {

		    XFig::Polygon p;
		    p.push_back(XFig::Point(-1,-1,1));
		    p.push_back(XFig::Point(1,-1,1));
		    p.push_back(XFig::Point(1,1,1));
		    p.push_back(XFig::Point(-1,1,1));

		    IntegerVectorList N2=orientedNeighbours(N,A[0]);
		    for(IntegerVectorList::const_iterator i=N2.begin();i!=N2.end();i++)
		      {
			IntegerVector N=mRot2.transposed().vectormultiply(*i);
			XFig::Point n(N[0],N[1],N[2]);
			fprintf(Stderr,"(%f,%f,%f)\n",n.x,n.y,n.z);

			p=xfig.intersect(p,n);
		      }
		    xfig.drawPolygon(p,1);
		  }
	      }
	    if(!A3)
	      xfig.drawDot((10000*x)/steps,(10000*y)/steps,colorMap(A1,A2,A3),radius);

	    //	    fprintf(Stderr,"test\n");
	  }
      }
  }

  void signature(XFig &xfig, int x, int y)
  {
    for(int i=0;i<4;i++)
      {
	bool A1=i&1;
	bool A2=i&2;
	bool A3=i&4;

	char s[100];
	sprintf(s,"%s,%s",(A1)?"A1":"a1",(A2)?"A2":"a2");

	xfig.drawString(x+500, 950*i +y,string(s),60);
	xfig.drawDot(x, 950*i +y-300,colorMap(A1,A2,A3),250);
      }
  }
  int main()
  {
    LpSolver::printList(Stderr);
    lpSetSolver("cddgmp");

    FileParser P(Stdin);

    IntegerVectorList ivl=P.parseIntegerVectorList();
    IntegerMatrix A=rowsToIntegerMatrix(ivl);

    XFig xfig(Stdout);

    int steps=optionNumberOfSteps.getValue();
    assert(steps>0);
    radius= 2500/steps;

    if(A.getWidth()==2)
      {
	signature(xfig,12000,-15000);

	xfig.drawString(16000, -15000,string(" x1   x2"),60);
	xfig.drawString(16000, -15000+900,string(" y1   y2"),60);
	for(int i=2;i<A.getHeight();i++)
	  printVector(xfig,16000, 950*i-15000, A[i]);

	//    int rotX=2;
	//    int rotY=3;


	for(int rotX=0;rotX<4;rotX++)
	  for(int rotY=0;rotY<4;rotY++)
	    {

	      xfig.setOffset(-rotX*2*10000,-rotY*2*10000);
	      box(xfig);

	      IntegerMatrix mRot=rowsToIntegerMatrix(StringParser("{(0,-1),(1,0)}").parseIntegerVectorList());
	      IntegerMatrix mRotX=IntegerMatrix::identity(2);
	      IntegerMatrix mRotY=IntegerMatrix::identity(2);
	      AsciiPrinter(Stderr).printVectorList(mRot.getRows());
	      AsciiPrinter(Stderr).printVectorList(mRotY.getRows());
	      AsciiPrinter(Stderr).printVectorList(mRotX.getRows());
	      for(int i=0;i<rotX;i++)
		{
		  IntegerMatrix A=mRotX;
		  IntegerMatrix B=mRot;
		  mRotX=A*B;
		}
	      for(int i=0;i<rotY;i++)
		mRotY=mRotY*mRot;

	      for(int x=-steps;x<=steps;x++)
		{
		  IntegerVector v(2);
		  v[1]=steps;
		  v[0]=x;
		  A[0]=mRotX.vectormultiply(v);

		  for(int y=-steps;y<=steps;y++)
		    {
		      IntegerVector v(2);
		      v[1]=steps;
		      v[0]=y;
		      A[1]=mRotY.vectormultiply(v);

		      if(x==-steps && y==-steps && rotY==0)
			{
			  printVector(xfig,-11000,11400,A[0]);
			}
		      if(x==-steps && y==-steps && rotX==0)
			{
			  printVector(xfig,11000,-10000+300,A[1]);
			}

		      /*	    if(x==-steps && y==-steps)
				    {
				    printVector(xfig,-10000,-10000,A[0]);
				    printVector(xfig,-10000,-10000+800,A[1]);
				    }
		      */



		      bool A3=false;
		      bool A1=satisfiesA1(A);
		      bool A2=satisfiesA2(A);
		      if(A1 && A2)A3=satisfiesA3(A);

		      if(A1 && A2 && A3)
			{
			  IntegerVectorList N=neighbours(A);

			  XFig::Polygon p;
			  p.push_back(XFig::Point(-1,-1,1));
			  p.push_back(XFig::Point(1,-1,1));
			  p.push_back(XFig::Point(1,1,1));
			  p.push_back(XFig::Point(-1,1,1));


			  for(int j=0;j<2;j++)
			    {
			      IntegerVectorList N2=orientedNeighbours(N,A[j]);
			      for(IntegerVectorList::const_iterator i=N2.begin();i!=N2.end();i++)
				{
				  IntegerVector N;
				  if(j==0)
				    N=mRotX.transposed().vectormultiply(*i);
				  else
				    N=mRotY.transposed().vectormultiply(*i);
				  XFig::Point n(N[0],N[0],N[1]);
				  if(j==0)
				    {
				      n.y=0;
				    }
				  else
				    {
				      n.x=0;
				    }
				  fprintf(Stderr,"(%f,%f,%f)\n",n.x,n.y,n.z);

				  p=xfig.intersect(p,n);
				}
			    }
			  xfig.drawPolygon(p,1);
			}
		      else
			xfig.drawDot((10000*x)/steps,(10000*y)/steps,colorMap(A1,A2,A3),radius);
		    }

		}
	    }
      }
    else if(A.getWidth()==3)
      {
	signature(xfig,-41000,-15000);

	xfig.drawString(-36000, -15000,string("   x   y   z"),60);
	for(int i=1;i<A.getHeight();i++)
	  printVector(xfig,-36000, 950*i-15000, A[i]);


	IntegerMatrix mRot=rowsToIntegerMatrix(StringParser("{(0,0,-1),(0,1,0),(1,0,0)}").parseIntegerVectorList());
	IntegerMatrix mRot2=IntegerMatrix::identity(3);

	for(int i=0;i<4;i++)
	  {
	    xfig.setOffset(-i*2*10000,0);
	    box(xfig);
	    renderFace(A,mRot2,xfig,steps);
	    mRot2=mRot*mRot2;
	  }

	IntegerMatrix mRot3=rowsToIntegerMatrix(StringParser("{(1,0,0),(0,0,1),(0,-1,0)}").parseIntegerVectorList());
	mRot2=mRot*mRot3;
	//	mRot2=mRot3*mRot2;
	//	mRot2=mRot3*mRot2;
	xfig.setOffset(-20000,20000);
	box(xfig);
	renderFace(A,mRot2,xfig,steps);
	mRot2=mRot2*mRot3;
	mRot2=mRot2*mRot3;
	xfig.setOffset(-20000,-20000);
	box(xfig);
	renderFace(A,mRot2,xfig,steps);

      }

    /*    IntegerVectorList N=neighbours(A);
    AsciiPrinter Q(Stdout);
    Q.printVectorList(N);
    fprintf(Stdout,"A1 satisfied:%i\n",satisfiesA1(A));
    fprintf(Stdout,"A2 satisfied:%i\n",satisfiesA2(A));
    fprintf(Stdout,"A3 satisfied:%i\n",satisfiesA3(A));
    */
    return 0;
  }
};

static ScarfVisualizeApplication theApplication;
