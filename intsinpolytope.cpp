#include "intsinpolytope.h"

#include "latticeideal.h"
#include "printer.h"
#include "linalg.h"
#include "field_rationals.h"

static IntegerVectorList::const_iterator findImproving(IntegerVectorList const &b, IntegerVector const &x)
{
  IntegerVectorList::const_iterator i=b.begin();

  while(i!=b.end())
    {
      if(i->divides(x))break;
      i++;
    }
  return i;
}

static int rek(IntegerVectorList const &b, IntegerVector &x, IntegerVectorList &output)
{
  //  fprintf(stdout,"rek\n");

  //  AsciiPrinter(Stdout).printVector(x);
  int ret=1;
  output.push_back(x);

  for(IntegerVectorList::const_iterator i=b.begin();i!=b.end();i++)
    {
      x+=*i;
      if(x.isNonNegative())
	if(findImproving(b,x)==i)ret+=rek(b,x,output);
      x-=*i;
    }
  return ret;
}


bool solveIntegerProgramIneq(IntegerMatrix const &M, IntegerVector const &rightHandSide, IntegerVector &solution)
{
  int d=M.getHeight();
  int n=M.getWidth();

  IntegerMatrix M2(d,n+d);
  for(int i=0;i<d;i++)
    {
      for(int j=0;j<n;j++)
	M2[i][j]=M[i][j];
      M2[i][n+i]=-1;
    }
  

  AsciiPrinter P(Stderr);
  P.printVectorList(M2.getRows());

  IntegerVectorList b=latticeIdealRevLex(M2);

  IntegerVector ret=IntegerVector(d+n);

  {
    IntegerVectorList::const_iterator i;
    while((i=findImproving(b,ret))!=b.end())
      {
	ret-=*i;
      }
  }

  solution=ret.subvector(n,n+d);

  return ret.subvector(0,n).isZero();
}


IntegerVectorList intsInPolytopeGivenIneqAndPt(IntegerMatrix const &M, IntegerVector const &rightHandSide, IntegerVector const &p)
{
  IntegerVectorList b=latticeIdealRevLex(M);

  //  AsciiPrinter(Stdout).printVectorList(p);

  IntegerVector p2=rightHandSide-M.vectormultiply(p);

  {
    IntegerVectorList::const_iterator i;
    while((i=findImproving(b,p2))!=b.end())
      {
	p2-=*i;
      }
  }

  IntegerVectorList points;
  rek(b,p2,points);


  FieldMatrix Mf=integerMatrixToFieldMatrix(M,Q);

  AsciiPrinter PP(Stderr);
  FieldMatrix solver=Mf.solver();
  solver.printMatrix(PP);
  IntegerVectorList ret;
  for(IntegerVectorList::const_iterator i=points.begin();i!=points.end();i++)
    {

      FieldVector temp=solver.canonicalize(concatenation(integerVectorToFieldVector(rightHandSide-*i,Q),
										 FieldVector(Q,
											     Mf.getWidth()
											     )));
      ret.push_back(fieldVectorToIntegerVector(temp.subvector(Mf.getHeight(),Mf.getHeight()+Mf.getWidth())));
    }



  //  AsciiPrinter(Stdout).printVectorList(points);

  
  

  return ret;
}

IntegerVectorList intsInPolytopeGivenIneq(IntegerMatrix const &M, IntegerVector const &rightHandSide)
{
  IntegerVectorList ret;
  IntegerVector solution;
  if(solveIntegerProgramIneq(M,rightHandSide,solution))
    intsInPolytopeGivenIneqAndPt(M,rightHandSide,solution);

  return ret;
}
