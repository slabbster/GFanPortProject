#include "lp_cdd.h"
#include "setoper.h"
#include "cdd.h"
#include "cdd_f.h"
#include "termorder.h"
#include "printer.h"
#include "log.h"

//--------------------------------------------------
// LpSolverCdd (double precision)
//--------------------------------------------------

static ddf_MatrixPtr vectorList2Matrix(int n, const IntegerVectorList &g, ddf_ErrorType *Error)
{
  ddf_MatrixPtr M=NULL;
  ddf_rowrange m_input,i;
  ddf_colrange d_input,j;
  ddf_RepresentationType rep=ddf_Inequality;
  ddf_boolean found=ddf_FALSE, newformat=ddf_FALSE, successful=ddf_FALSE;
  char command[ddf_linelenmax], comsave[ddf_linelenmax];
  ddf_NumberType NT;

  (*Error)=ddf_NoError;

  rep=ddf_Inequality; newformat=ddf_TRUE;

  m_input=g.size();
  d_input=n+1;//g.begin()->size()+1;

  NT=ddf_Rational;

  M=ddf_CreateMatrix(m_input, d_input);
  M->representation=rep;
  M->numbtype=NT;

  IntegerVectorList::const_iterator I=g.begin();
  for (i = 0; i < m_input; i++) {
    ddf_set_si(M->matrix[i][0],0);
    for (j = 1; j < d_input; j++) {
      ddf_set_si(M->matrix[i][j],(*I)[j-1]);
    }
    I++;
  }

  successful=ddf_TRUE;

  return M;
}

static void cddinit()
{
  static bool initialized;
  if(!initialized)
    {
      ddf_set_global_constants();  /* First, this must be called. */
      initialized=true;
    }
}

bool LpSolverCdd::isFacet(const IntegerVectorList &g, IntegerVectorList::const_iterator i)
{
  bool ret;
  int index;
  ddf_MatrixPtr M=NULL,M2=NULL,M3=NULL;
  ddf_colrange d;
  ddf_ErrorType err=ddf_NoError;
  ddf_rowset redrows,linrows,ignoredrows, basisrows;
  ddf_colset ignoredcols, basiscols;
  //  long rank;
  mytype val;
  ddf_DataFileType inputfile;
  FILE *reading=NULL;


  cddinit();

  //  ddf_init(val);

  M=vectorList2Matrix(i->size(), g, &err);

  if (err!=ddf_NoError) goto _L99;

  d=M->colsize;

  //  redrows=ddf_RedundantRows(M, &err);
  //  set_fwrite(Stderr, redrows);


  index=0;
  for(IntegerVectorList::const_iterator J=g.begin();J!=g.end()&&J!=i;J++){index++;}
  if(index==g.size())assert(0);


  static ddf_Arow temp;

  ddf_InitializeArow(i->size()+1,&temp);

  ret= !ddf_Redundant(M,index+1,temp,&err);

  ddf_FreeMatrix(M);
  ddf_FreeArow(i->size()+1,temp);
  //  ddf_FreeArow(i->size(),temp);

  if (err!=ddf_NoError) goto _L99;
  return ret;
 _L99:
  assert(0);
  return false;

}

static LpSolverCdd theLpSolverCdd;


//--------------------------------------------------
// LpSolverCddGmp (exact arithmetics)
//--------------------------------------------------


static void cddinitGmp()
{
  static bool initialized;
  if(!initialized)
    {
      dd_set_global_constants();  /* First, this must be called. */
      initialized=true;
    }
}


static dd_MatrixPtr vectorList2MatrixGmp(int n, const IntegerVectorList &g, dd_ErrorType *Error)
{
  dd_MatrixPtr M=NULL;
  dd_rowrange m_input,i;
  dd_colrange d_input,j;
  dd_RepresentationType rep=dd_Inequality;
  dd_boolean found=dd_FALSE, newformat=dd_FALSE, successful=dd_FALSE;
  char command[dd_linelenmax], comsave[dd_linelenmax];
  dd_NumberType NT;

  (*Error)=dd_NoError;

  rep=dd_Inequality; newformat=dd_TRUE;

  if(n==-1)
    {
      assert(g.size());
      n=g.begin()->size();
    }
  m_input=g.size();
  //  d_input=g.begin()->size()+1;
  d_input=n+1;
  if(m_input)
    {
      if(g.begin()->size()!=n)
	{
	  AsciiPrinter(Stderr).printVectorList(g);
	}


      assert(g.begin()->size()==n);
    }

  NT=dd_Rational;

  M=dd_CreateMatrix(m_input, d_input);
  M->representation=rep;
  M->numbtype=NT;

  IntegerVectorList::const_iterator I=g.begin();
  for (i = 0; i < m_input; i++) {
    dd_set_si(M->matrix[i][0],0);
    for (j = 1; j < d_input; j++) {
      dd_set_si(M->matrix[i][j],(*I)[j-1]);
    }
    I++;
  }

  successful=dd_TRUE;

  return M;
}


static dd_MatrixPtr vectorList2MatrixIncludingFirstColumnGmp(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations, dd_ErrorType *Error)
{
  dd_MatrixPtr M=NULL;
  dd_rowrange m_input,i;
  dd_colrange d_input,j;
  dd_RepresentationType rep=dd_Inequality;
  dd_boolean found=dd_FALSE, newformat=dd_FALSE, successful=dd_FALSE;
  //  char command[dd_linelenmax], comsave[dd_linelenmax];
  dd_NumberType NT;

  (*Error)=dd_NoError;

  int numberOfEquations=equations.size();
  int numberOfInequalities=inequalities.size();
  m_input=numberOfEquations+numberOfInequalities;
  assert(m_input>0);

  rep=dd_Inequality; newformat=dd_TRUE;

  d_input=n;
  assert(d_input>0);

  NT=dd_Rational;

  M=dd_CreateMatrix(m_input, d_input);
  M->representation=rep;
  M->numbtype=NT;

  IntegerVectorList::const_iterator I=inequalities.begin();
  for (i = 0; i < numberOfInequalities; i++) {
    for (j = 0; j < d_input; j++)dd_set_si(M->matrix[i][j],(*I)[j]);
    I++;
  }
  I=equations.begin();
  for (; i < m_input; i++) {
    set_addelem(M->linset, i+1);
    for (j = 0; j < d_input; j++)dd_set_si(M->matrix[i][j],(*I)[j]);
    I++;
  }

  successful=dd_TRUE;

  return M;
}

bool LpSolverCddGmp::hasHomogeneousSolution(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations)
{
  if(n==0)return false;
  int nrows=inequalities.size()+equations.size();
  if(nrows==0)return true;

  dd_LPSolverType solver=dd_DualSimplex;
  dd_MatrixPtr A=NULL;
  dd_LPSolutionPtr lps1;
  dd_ErrorType err=dd_NoError;
  int ret=0;
  cddinitGmp();

  dd_LPPtr lp,lp1;

  A=vectorList2MatrixIncludingFirstColumnGmp(n, inequalities, equations, &err);

  if (err!=dd_NoError) goto _L99;

  A->objective=dd_LPmax;
  lp=dd_Matrix2LP(A, &err);
  if (err!=dd_NoError) goto _L99;


  //  dd_WriteMatrix(Stderr,A);

  /* Find an interior point with cdd LP library. */

  lp1=dd_MakeLPforInteriorFinding(lp);
  dd_LPSolve(lp1,solver,&err);
  if (err!=dd_NoError) goto _L99;

  //  dd_WriteMatrix(Stderr,A);
  //  dd_WriteLP(Stderr,lp);
  //  dd_WriteLP(Stderr,lp1);

  /* Write an interior point. */
  lps1=dd_CopyLPSolution(lp1);
  //  dd_WriteLPSolution(Stderr,lps1);
  /*  {
    fprintf(Stderr,"(");
    for (int j=1; j <(lps1->d); j++) {
      if(j!=1)fprintf(Stderr,", ");
      dd_WriteNumber(Stderr,lps1->sol[j]);
    }
    fprintf(Stderr,")\n");
    }*/

  //  dd_WriteNumber(Stderr,lps1->optvalue);

  if(dd_Positive(lps1->optvalue))ret=1;
  if(dd_Negative(lps1->optvalue))ret=-1;

  //  fprintf(Stderr,"ret=%i",ret);

  dd_FreeLPData(lp);
  dd_FreeLPSolution(lps1);
  dd_FreeLPData(lp1);
  dd_FreeMatrix(A);

  return ret>=0;
 _L99:
  assert(0);
  return 0;
}


bool LpSolverCddGmp::isFacet(const IntegerVectorList &g, IntegerVectorList::const_iterator i)
{
  bool ret;
  int index;
  dd_MatrixPtr M=NULL,M2=NULL,M3=NULL;
  dd_colrange d;
  dd_ErrorType err=dd_NoError;
  dd_rowset redrows,linrows,ignoredrows, basisrows;
  dd_colset ignoredcols, basiscols;
  //  long rank;
  mytype val;
  dd_DataFileType inputfile;
  FILE *reading=NULL;


  cddinitGmp();


  //  dd_init(val);

  M=vectorList2MatrixGmp(i->size(),g, &err);

  if (err!=dd_NoError) goto _L99;

  d=M->colsize;

  //  redrows=dd_RedundantRows(M, &err);
  //  set_fwrite(Stderr, redrows);


  index=0;
  for(IntegerVectorList::const_iterator J=g.begin();J!=g.end()&&J!=i;J++){index++;}
  if(index==g.size())assert(0);


  static dd_Arow temp;

  dd_InitializeArow(i->size()+1,&temp);

  ret= !dd_Redundant(M,index+1,temp,&err);


  dd_FreeMatrix(M);
  dd_FreeArow(i->size()+1,temp);

  if (err!=dd_NoError) goto _L99;
  return ret;
 _L99:
  assert(0);
  return false;

}

int staticInteriorPoint(int n, mpq_t *point, const IntegerVectorList &g, bool strictlyPositive, IntegerVector const *equalitySet=0)
{// copy-paste from testshoot.c in cdd
  //  AsciiPrinter(Stderr).printVectorList(g);
  //  fprintf(Stderr,"strictlyPositive:%i\n",strictlyPositive);

  //if(equalitySet)  AsciiPrinter(Stderr).printVector(*equalitySet);
  dd_LPSolverType solver=dd_DualSimplex;
  dd_MatrixPtr A=NULL;
  dd_LPSolutionPtr lps1;
  dd_ErrorType err=dd_NoError;
  int ret=0;
  cddinitGmp();

  dd_LPPtr lp,lp1;

  assert(g.begin()!=g.end());
  if(strictlyPositive)
    {
      int n=g.begin()->size();
      IntegerVectorList G=g;
      for(int i=0;i<n;i++)
	G.push_back(IntegerVector::standardVector(n,i));
      A=vectorList2MatrixGmp(n,G, &err);
    }
  else
    A=vectorList2MatrixGmp(n, g, &err);
  if (err!=dd_NoError) goto _L99;

  if(equalitySet)
    {
      for(int i=0;i<g.size();i++)
	if(!(*equalitySet)[i])
	  dd_set_si(A->matrix[i][0],-1);

      assert(g.size()>=equalitySet->size());

      for(int i=0;i<equalitySet->size();i++)
	if((*equalitySet)[i])set_addelem(A->linset, i+1);
    }

  A->objective=dd_LPmax;
  lp=dd_Matrix2LP(A, &err);
  if (err!=dd_NoError) goto _L99;


  //  dd_WriteMatrix(Stderr,A);

  /* Find an interior point with cdd LP library. */

  lp1=dd_MakeLPforInteriorFinding(lp);
  dd_LPSolve(lp1,solver,&err);
  if (err!=dd_NoError) goto _L99;

  //  dd_WriteMatrix(Stderr,A);
  //  dd_WriteLP(Stderr,lp1);

  /* Write an interior point. */
  lps1=dd_CopyLPSolution(lp1);

  if(dd_Positive(lps1->optvalue))ret=1;
  if(dd_Negative(lps1->optvalue))ret=-1;

  //  fprintf(Stderr,"ret=%i",ret);

  if (ret>=0)
    //  if (ret>0)
    for (int j=1; j <(lps1->d)-1; j++)
      mpq_set(point[j-1],lps1->sol[j]);

  dd_FreeLPData(lp);
  dd_FreeLPSolution(lps1);
  dd_FreeLPData(lp1);
  dd_FreeMatrix(A);
  return ret;
 _L99:
  assert(0);
  return 0;
}

int staticRelativeInteriorPoint(int n, mpq_t *point, const IntegerVectorList &g, bool strictlyPositive, IntegerVector const *equalitySet=0)
{// copy-paste from testshoot.c in cdd
  //  AsciiPrinter(Stderr).printVectorList(g);
  //  fprintf(Stderr,"strictlyPositive:%i\n",strictlyPositive);

  //if(equalitySet)  AsciiPrinter(Stderr).printVector(*equalitySet);
  dd_LPSolverType solver=dd_DualSimplex;
  dd_MatrixPtr A=NULL;
  dd_LPSolutionPtr lps1;
  dd_ErrorType err=dd_NoError;
  int ret=0;
  cddinitGmp();

  dd_LPPtr lp,lp1;

  //  assert(g.begin()!=g.end());
  if(strictlyPositive)
    {
      //  int n=g.begin()->size();
      IntegerVectorList G=g;
      for(int i=0;i<n;i++)
	G.push_back(IntegerVector::standardVector(n,i));
      A=vectorList2MatrixGmp(n, G, &err);
    }
  else
    A=vectorList2MatrixGmp(n, g, &err);
  if (err!=dd_NoError) goto _L99;

  if(equalitySet)
    {
      for(int i=0;i<g.size();i++)
	if(!(*equalitySet)[i])
	  dd_set_si(A->matrix[i][0],-1);

      assert(g.size()>=equalitySet->size());

      for(int i=0;i<equalitySet->size();i++)
	if((*equalitySet)[i])set_addelem(A->linset, i+1);
    }

  A->objective=dd_LPmax;
  lp=dd_Matrix2LP(A, &err);
  if (err!=dd_NoError) goto _L99;


  //  dd_WriteMatrix(Stderr,A);

  /* Find an interior point with cdd LP library. */

  lp1=dd_MakeLPforInteriorFinding(lp);
  dd_LPSolve(lp1,solver,&err);
  if (err!=dd_NoError) goto _L99;

  //  dd_WriteMatrix(Stderr,A);
  //  dd_WriteLP(Stderr,lp1);

  /* Write an interior point. */
  lps1=dd_CopyLPSolution(lp1);

  if(dd_Positive(lps1->optvalue))ret=1;
  if(dd_Negative(lps1->optvalue))ret=-1;

  //  fprintf(Stderr,"ret=%i",ret);

  if (ret>=0)
    //  if (ret>0)
    for (int j=1; j <(lps1->d)-1; j++)
      mpq_set(point[j-1],lps1->sol[j]);

  dd_FreeLPData(lp);
  dd_FreeLPSolution(lps1);
  dd_FreeLPData(lp1);
  dd_FreeMatrix(A);
  return ret;
 _L99:
  assert(0);
  return 0;
}

bool LpSolverCddGmp::hasInteriorPoint(const IntegerVectorList &g, bool strictlyPositive, IntegerVector const *equalitySet)
{
  assert(!g.empty());
  int n=g.begin()->size();
  mpq_t *point = new mpq_t [n];
  for(int i=0;i<n;i++)mpq_init(point[i]);

  int ret=staticInteriorPoint(n, point,g,strictlyPositive,equalitySet);

  // fprintf(Stderr,"%i\n",ret);

  for(int i=0;i<n;i++)mpq_clear(point[i]);
  delete [] point;

  if(equalitySet)return ret>=0; //THIS NEEDS TO BE FIXED
  return ret>0;
}


IntegerVector arrayToIntegerVector(mpq_t *point, int n)
{
  IntegerVector ret(n);

  for (int j=0; j <n; j++)
    {
      int den;
      int num;

      if((!mpz_fits_sint_p(mpq_denref(point[j])))||(!mpz_fits_sint_p(mpq_numref(point[j]))))
	{
	  fprintf(stderr,"INTEGER OVERFLOW IN POLYHEDRAL COMPUTATION\n");
	  assert(0);
	}
      den=mpz_get_si(mpq_denref(point[j]));
      num=mpz_get_si(mpq_numref(point[j]));

      assert(den==1);
      ret[j]=num;
    }

  return ret;
}

void scaleToIntegerVector(mpq_t *point, int n)
{
  mpz_t lcm;
  mpz_t gcd;
  mpz_init_set_ui(lcm, 1);
  mpz_init_set_ui(gcd, 0);

  for(int j=0;j<n;j++)
    {
      mpz_lcm(lcm,lcm,mpq_denref(point[j]));
      mpz_gcd(gcd,gcd,mpq_numref(point[j]));
    }

  if(mpz_sgn(gcd)!=0)
if(mpz_cmp(gcd,lcm)!=0)
    {
      assert(mpz_sgn(gcd)>0);
      mpq_t scale;
      mpq_init(scale);
      mpq_set_den(scale,gcd);
      mpq_set_num(scale,lcm);
      mpq_canonicalize(scale); //according to the gmp manual this call is needed, although it slows things down a bit
      for(int j=0;j<n;j++)
	{
	  mpq_mul(point[j],point[j],scale);
	}
      mpq_clear(scale);
    }
  mpz_clear(lcm);
  mpz_clear(gcd);
}


IntegerVector LpSolverCddGmp::relativeInteriorPoint(int n, const IntegerVectorList &g, IntegerVector const *equalitySet)
{
  //  assert(!g.empty());
  //  int n=g.begin()->size();
  mpq_t *point = new mpq_t [n];
  for(int i=0;i<n;i++)mpq_init(point[i]);

  int ret=staticRelativeInteriorPoint(n, point,g,false,equalitySet);

  assert(ret>=0);//-- any cone has a relative interior point
  //  if (ret>0){
  if (ret>=0)
  {
    //    fprintf(stderr,"TEST1\n");
    scaleToIntegerVector(point,n);
    //   fprintf(stderr,"TEST2\n");
  }


  IntegerVector result=arrayToIntegerVector(point,n);


  for(int i=0;i<n;i++)mpq_clear(point[i]);
  delete [] point;

  return result;
}


bool LpSolverCddGmp::interiorPoint(const IntegerVectorList &g, IntegerVector &result, bool strictlyPositive, IntegerVector const *equalitySet)
{
  int n=g.begin()->size();
  mpq_t *point = new mpq_t [n];
  for(int i=0;i<n;i++)mpq_init(point[i]);

  int ret=staticInteriorPoint(n, point,g,strictlyPositive,equalitySet);

  //  if (ret>0){
  if (ret>=0)
  {
    scaleToIntegerVector(point,n);
  }


  result=arrayToIntegerVector(point,n);
  //  if (ret>0){
  /*  if (ret>=0){
    fprintf(f,"(");
    for (int j=0; j <n; j++) {
      if(j!=0)fprintf(f,", ");
      dd_WriteNumber(f,point[j]);
    }
    fprintf(f,")\n");
  }
  if (ret<0)
    fprintf(f,"The feasible region is empty.\n");
  if (ret==0)
    fprintf(f,"The feasible region is nonempty but has no interior point.\n");
  */


  for(int i=0;i<n;i++)mpq_clear(point[i]);
  delete [] point;

  return (ret>=0);
}

// the following two routines are static to avoid including gmp.h in the header file. Maybe that should be changed...

static bool lexicographicShootCompare(IntegerVector const &a, IntegerVector const &b, mpq_t const &aDot, mpq_t const &bDot, mpq_t &aTemp, mpq_t &bTemp)
{
  int n=a.size();
  assert(b.size()==n);
  for(int i=0;i<n;i++)
    {
      mpq_set_si(aTemp,a[i],1);
      mpq_set_si(bTemp,b[i],1);
      mpq_mul(aTemp,aTemp,bDot);
      mpq_mul(bTemp,bTemp,aDot);
      int cmp=mpq_cmp(aTemp,bTemp);
      if(cmp>0)return true;
      if(cmp<0)return false;
    }
  assert(0);
  return false;
}

static void computeDotProduct(mpq_t &ret, IntegerVector const &iv, mpq_t const *qv, mpq_t &temp)
{
  mpq_set_si(ret,0,1);
  int n=iv.size();
  for(int i=0;i<n;i++)
    {
      mpq_set_si(temp,iv[i],1);
      mpq_mul(temp,temp,qv[i]);
      mpq_add(ret,ret,temp);
    }
}

IntegerVectorList::const_iterator LpSolverCddGmp::shoot(const IntegerVectorList &g)
{
  //  fprintf(Stderr,"shoot begin\n");

  //  AsciiPrinter(Stderr).printVectorList(g);

  int n=g.begin()->size();
  mpq_t *point = new mpq_t [n];
  for(int i=0;i<n;i++)mpq_init(point[i]);

  //  fprintf(Stderr,"\nVektor list with no interior point:\n");
  //  AsciiPrinter(Stderr).printVectorList(g);


  int ret=staticInteriorPoint(n, point,g,true);

  //fprintf(Stderr,"Interior point\n");
  //  AsciiPrinter(Stderr).printIntegerVector(point);

  assert(ret>0);

  IntegerVectorList::const_iterator bestIterator=g.end();
  mpq_t tempA;
  mpq_init(tempA);
  mpq_t tempB;
  mpq_init(tempB);
  mpq_t bestDotProduct;
  mpq_init(bestDotProduct);
  mpq_t currentDotProduct;
  mpq_init(currentDotProduct);

  IntegerVectorList::const_iterator i=g.begin();
  for(;i!=g.end();i++)
    if(LexicographicTermOrder()(*i,*i-*i))
      {
	computeDotProduct(bestDotProduct,*i,point,tempA);
	bestIterator=i;
	break;
      }
  //  fprintf(Stderr,"A\n");
  //if(bestIterator!=g.end())
  //  AsciiPrinter(Stderr).printVector(*bestIterator);
  if(i!=g.end())
    for(i++;i!=g.end();i++)
      {
	computeDotProduct(currentDotProduct,*i,point,tempA);
	if(lexicographicShootCompare(*bestIterator,*i,bestDotProduct,currentDotProduct,tempA,tempB))
	  {
	    bestIterator=i;
	    mpq_set(bestDotProduct,currentDotProduct);
	  }
      }
  mpq_clear(currentDotProduct);
  mpq_clear(bestDotProduct);
  mpq_clear(tempB);
  mpq_clear(tempA);

  for(int i=0;i<n;i++)mpq_clear(point[i]);
  delete [] point;

  //  fprintf(Stderr,"shoot end\n");
  //if(bestIterator!=g.end())
  //  AsciiPrinter(Stderr).printVector(*bestIterator);
  return bestIterator;
}


//-----------------------------------------
// Positive vector in kernel
//-----------------------------------------

bool LpSolverCddGmp::positiveVectorInKernel(const IntegerVectorList &g_, IntegerVector *result)
{// copy-paste from testshoot.c in cdd
  //  AsciiPrinter(Stderr).printVectorList(g);
  //  fprintf(Stderr,"strictlyPositive:%i\n",strictlyPositive);


  IntegerVectorList g=g_;


  //  AsciiPrinter(Stderr).printVectorList(g);

  int dim2=g.size();

  assert(g.size()!=0);
  int dimension=g.begin()->size();

  for(int i=0;i<dimension;i++)
    g.push_back(IntegerVector::standardVector(dimension,i));

  dd_LPSolverType solver=dd_DualSimplex;
  dd_MatrixPtr A=NULL;
  dd_LPSolutionPtr lps1;
  dd_ErrorType err=dd_NoError;
  //  int ret=0;
  cddinitGmp();

  dd_LPPtr lp;

  A=vectorList2MatrixGmp(dimension, g, &err);

  for(int i=0;i<dim2;i++)set_addelem(A->linset, i+1);//???
  // set objective function
  for (int i = 0; i < dimension; i++)
    dd_set_si(A->rowvec[i+1],-1);


  // set right hand side
  for (int i = 0; i < dimension; i++)
    dd_set_si(A->matrix[i+dim2][0],-1);


  if (err!=dd_NoError) goto _L99;
  A->objective=dd_LPmax;


  /*  for(int i=0;i<dimension;i++)
    {
      for(int j=0;j<dimension;j++)
	dd_WriteNumber(f,point[j]);

      fprintf(Stderr,"\n");
    }
  */

  //  dd_WriteMatrix(Stderr,A);
  //  assert(0);

  lp=dd_Matrix2LP(A, &err);
  if (err!=dd_NoError) goto _L99;

//  dd_WriteLP(Stderr,lp);
  /* Find an interior point with cdd LP library. */

  dd_LPSolve(lp,solver,&err);


  if (err!=dd_NoError) goto _L99;

  //  IF INFEASIBLE RETURN FALSE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(lp->LPS==dd_Inconsistent)
    {
      dd_FreeLPData(lp);
      dd_FreeMatrix(A);

      return false;
    }
//  dd_WriteLPResult(stdout,lp,err);

  /* Write an interior point. */
  lps1=dd_CopyLPSolution(lp);

  if (dd_Positive(lps1->optvalue))
    {
//      fprintf(Stderr,"Returning false\n");
      return false;
    }




/*  {
    fprintf(Stderr,"(");
    for (int j=1; j <(lps1->d); j++) {
      if(j!=1)fprintf(Stderr,", ");
      dd_WriteNumber(Stderr,lps1->sol[j]);
    }
    fprintf(Stderr,")\n");
  }
*/

  //transform into integer vectors
  {
    int n=dimension;
    dd_Arow point=lps1->sol+1;

    mpz_t lcm;
    mpz_t gcd;
    mpz_init_set_ui(lcm, 1);
    mpz_init_set_ui(gcd, 0);

    for(int j=0;j<n;j++)
      {
	mpz_lcm(lcm,lcm,mpq_denref(point[j]));
	mpz_gcd(gcd,gcd,mpq_numref(point[j]));
      }

    mpq_t scale;
    mpq_init(scale);
    mpq_set_den(scale,gcd);
    mpq_set_num(scale,lcm);
    for(int j=0;j<n;j++)
      {
	mpq_mul(point[j],point[j],scale);
      }

    mpz_clear(lcm);
    mpz_clear(gcd);
  }



/*  {
    fprintf(Stderr,"(");
    for (int j=1; j <(lps1->d); j++) {
      if(j!=1)fprintf(Stderr,", ");
      dd_WriteNumber(Stderr,lps1->sol[j]);
    }
    fprintf(Stderr,")\n");
  }
*/
  for (int j=1; j <(lps1->d); j++)
    {
      int den;
      int num;

      den=mpz_get_si(mpq_denref(lps1->sol[j]));
      num=mpz_get_si(mpq_numref(lps1->sol[j]));

      assert(den==1);
      if(result)
	{
	  assert(j-1<result->size());
	  (*result)[j-1]=num;
	}
    }

  //  dd_WriteLPSolution(lps1);

  dd_FreeLPData(lp);
  dd_FreeLPSolution(lps1);
  dd_FreeMatrix(A);
  return true;
 _L99:
  assert(0);
  return 0;
}

//-----------------------------------------
// Rank of matrix
//-----------------------------------------

#include "subspace.h"
int LpSolverCddGmp::rankOfMatrix(const IntegerVectorList &g_)
{
  if(g_.size()==0)return 0;
  FieldMatrix temp=integerMatrixToFieldMatrix(rowsToIntegerMatrix(g_),Q);
  return temp.reduceAndComputeRank();
  return Subspace(g_).dimension();
  dd_rowset r=NULL;
  int result;
  IntegerVectorList g=g_;

  int dim2=g.size();

  assert(g.size()!=0);
  int dimension=g.begin()->size();

  /*  for(int i=0;i<dimension;i++)
    g.push_back(IntegerVector::standardVector(dimension,i));
  */
  dd_LPSolverType solver=dd_DualSimplex;
  dd_MatrixPtr A=NULL;
  dd_ErrorType err=dd_NoError;
  //  int ret=0;
  cddinitGmp();


  A=vectorList2MatrixGmp(dimension, g, &err);

  //  dd_WriteMatrix(Stderr,A);

  for(int i=0;i<g.size();i++)
  set_addelem(A->linset,i+1);
  //  for(int i=0;i<dim2;i++)set_addelem(A->linset, i+1);//???
  // set objective function
  for (int i = 0; i < dimension; i++)
    dd_set_si(A->rowvec[i+1],-1);


  // set right hand side
  /*for (int i = 0; i < dimension; i++)
    dd_set_si(A->matrix[i+dim2][0],-1);
  */

  if (err!=dd_NoError) goto _L99;
  A->objective=dd_LPmax;

  //fprintf(Stderr,"rank of matrix matrix:\n");

  //  dd_WriteMatrix(Stderr,A);

  r=dd_RedundantRows(A,&err);
  if (err!=dd_NoError) goto _L99;

  result=dim2-set_card(r);

//  fprintf(Stderr,"dim2==%i set_card(r)==%i\n",dim2,set_card(r));

  set_free(r);

  dd_FreeMatrix(A);
  return result;
 _L99:
  assert(0);
  return 0;
}

//-----------------------------------------
// Extreme Rays' Inequality Indices
//-----------------------------------------

// this procedure is take from cddio.c.
static void dd_ComputeAinc(dd_PolyhedraPtr poly)
{
/* This generates the input incidence array poly->Ainc, and
   two sets: poly->Ared, poly->Adom.
*/
  dd_bigrange k;
  dd_rowrange i,m1;
  dd_colrange j;
  dd_boolean redundant;
  dd_MatrixPtr M=NULL;
  mytype sum,temp;

  dd_init(sum); dd_init(temp);
  if (poly->AincGenerated==dd_TRUE) goto _L99;

  M=dd_CopyOutput(poly);
  poly->n=M->rowsize;
  m1=poly->m1;


  /*  fprintf(Stderr,"m1:%i\n",m1);
  fprintf(Stderr,"n:%i\n",poly->n);
  fprintf(Stderr,"m:%i\n",poly->m);
  */
   /* this number is same as poly->m, except when
      poly is given by nonhomogeneous inequalty:
      !(poly->homogeneous) && poly->representation==Inequality,
      it is poly->m+1.   See dd_ConeDataLoad.
   */
  poly->Ainc=(set_type*)calloc(m1, sizeof(set_type));
  for(i=1; i<=m1; i++) set_initialize(&(poly->Ainc[i-1]),poly->n);
  set_initialize(&(poly->Ared), m1);
  set_initialize(&(poly->Adom), m1);

  for (k=1; k<=poly->n; k++){
    for (i=1; i<=poly->m; i++){
      dd_set(sum,dd_purezero);
      for (j=1; j<=poly->d; j++){
        dd_mul(temp,poly->A[i-1][j-1],M->matrix[k-1][j-1]);
        dd_add(sum,sum,temp);
      }
      if (dd_EqualToZero(sum)) {
        set_addelem(poly->Ainc[i-1], k);
      }
    }
    if (!(poly->homogeneous) && poly->representation==dd_Inequality){
      if (dd_EqualToZero(M->matrix[k-1][0])) {
        set_addelem(poly->Ainc[m1-1], k);  /* added infinity inequality (1,0,0,...,0) */
      }
    }
  }

  for (i=1; i<=m1; i++){
    if (set_card(poly->Ainc[i-1])==M->rowsize){
      set_addelem(poly->Adom, i);
    }
  }
  for (i=m1; i>=1; i--){
    if (set_card(poly->Ainc[i-1])==0){
      redundant=dd_TRUE;
      set_addelem(poly->Ared, i);
    }else {
      redundant=dd_FALSE;
      for (k=1; k<=m1; k++) {
        if (k!=i && !set_member(k, poly->Ared)  && !set_member(k, poly->Adom) &&
            set_subset(poly->Ainc[i-1], poly->Ainc[k-1])){
          if (!redundant){
            redundant=dd_TRUE;
          }
          set_addelem(poly->Ared, i);
        }
      }
    }
  }
  dd_FreeMatrix(M);
  poly->AincGenerated=dd_TRUE;
_L99:;
  dd_clear(sum);  dd_clear(temp);
}


IntegerVectorList LpSolverCddGmp::extremeRaysInequalityIndices(const IntegerVectorList &inequalityList)
{
  int result;
  IntegerVectorList g=inequalityList;

  int dim2=g.size();

  //  assert(g.size()!=0);
  if(g.size()==0)return IntegerVectorList();

  int dimension=g.begin()->size();

  dd_MatrixPtr A=NULL;
  dd_ErrorType err=dd_NoError;

  cddinitGmp();

  A=vectorList2MatrixGmp(dimension, g, &err);

  //  dd_WriteMatrix(Stderr,A);

  dd_PolyhedraPtr poly;
  poly=dd_DDMatrix2Poly2(A, dd_LexMin, &err);


  //  dd_WritePolyFile(Stderr,poly);

  //  dd_WriteIncidence(Stderr,poly);
  if (poly->child==NULL || poly->child->CompStatus!=dd_AllFound) assert(0);
  if (poly->AincGenerated==dd_FALSE) dd_ComputeAinc(poly);


  //    dd_WriteIncidence(Stderr,poly);

  IntegerVectorList ret;

  //  fprintf(Stderr,"%i %i\n",poly->n,poly->m1);

  /*
    How do we interpret the cddlib output?  For a long ting gfan has
    been using poly->n as the number of rays of the cone and thus
    returned sets of indices that actually gave the lineality
    space. The mistake was then caught later in PolyhedralCone. On Feb
    17 2009 gfan was changed to check the length of each set to make
    sure that it does not specify the lineality space and only return
    those sets giving rise to rays.  This does not seem to be the best
    strategy and might even be wrong.
   */


  for (int k=1; k<=poly->n; k++)
    //for (int k=1; k<=poly->m1; k++)
    {
      int length=0;
      for (int i=1; i<=poly->m1; i++)
	if(set_member(k,poly->Ainc[i-1]))length++;
      if(length!=dim2)
	{
	  IntegerVector v(length);
	  int j=0;
	  for (int i=1; i<=poly->m1; i++)
	    if(set_member(k,poly->Ainc[i-1]))v[j++]=i-1;
	  ret.push_back(v);
	}
    }

  //  AsciiPrinter(Stderr).printVectorList(ret);
  //  dd_WriteIncidence(Stderr,poly);

  dd_FreeMatrix(A);
  dd_FreePolyhedra(poly);

  return ret;
 _L99:
  assert(0);
  return IntegerVectorList();
}




//-----------------------------------------
// Remove Redundant Rows
//-----------------------------------------

void LpSolverCddGmp::removeRedundantRows(IntegerVectorList *inequalities, IntegerVectorList *equalities, bool removeInequalityRedundancies)
{
  int numberOfEqualities=equalities->size();
  int numberOfInequalities=inequalities->size();
  int numberOfRows=numberOfEqualities+numberOfInequalities;


  //  AsciiPrinter(Stderr).printVectorList(*inequalities);
  //  AsciiPrinter(Stderr).printVectorList(*equalities);

  dd_rowset r=NULL;
  IntegerVectorList g=*inequalities;
  for(IntegerVectorList::const_iterator i=equalities->begin();i!=equalities->end();i++)
    g.push_back(*i);

  if(numberOfRows==0)return;//the full space, so description is alredy irredundant
  assert(numberOfRows>0);

  dd_LPSolverType solver=dd_DualSimplex;
  dd_MatrixPtr A=NULL;
  dd_ErrorType err=dd_NoError;

  cddinitGmp();

  A=vectorList2MatrixGmp(g.begin()->size(),g, &err);

  for(int i=numberOfInequalities;i<numberOfRows;i++)
    set_addelem(A->linset,i+1);

  //  dd_WriteMatrix(Stderr,A);

  //  for(int i=0;i<dim2;i++)set_addelem(A->linset, i+1);//???
  // set objective function
  //  for (int i = 0; i < dimension; i++)
  //    dd_set_si(A->rowvec[i+1],-1);


  // set right hand side
  /*for (int i = 0; i < dimension; i++)
    dd_set_si(A->matrix[i+dim2][0],-1);
  */


  IntegerVectorList newLin;
  IntegerVectorList newIn;
  int index=0;
  if (err!=dd_NoError) goto _L99;
  A->objective=dd_LPmax;

  //  dd_WriteMatrix(Stderr,A);

  //  r=dd_RedundantRows(A,&err);

  dd_rowset impl_linset;
  dd_rowset redset;
  dd_rowindex newpos;

  if(removeInequalityRedundancies)
    dd_MatrixCanonicalize(&A, &impl_linset, &redset, &newpos, &err);
  else
    dd_MatrixCanonicalizeLinearity(&A, &impl_linset, &newpos, &err);

  if (err!=dd_NoError) goto _L99;

  //  set_fwrite(stderr,redset);
  //  set_fwrite(stderr,impl_linset);

  // Maybe the following should be changed... what if cononicalize generates new rows?

  if(1)
    {
      //  dd_WriteMatrix(Stderr,A);
      int rowsize=A->rowsize;
      int n=A->colsize-1;
      // fprintf(Stderr,"rowsize: %i ColSize:%i\n",rowsize,n);
      for(int i=0;i<rowsize;i++)
	{
	  mpq_t *point = new mpq_t [n];
	  for(int j=0;j<n;j++)mpq_init(point[j]);

	  for(int j=0;j<n;j++)mpq_set(point[j],A->matrix[i][j+1]);

	  IntegerVector v=arrayToIntegerVector(point, n);
	  //  AsciiPrinter(Stderr).printVector(v);

	  for(int j=0;j<n;j++)mpq_clear(point[j]);
	  delete [] point;
	  if(set_member(i+1,A->linset))
	    newLin.push_back(v);
	  else
	    newIn.push_back(v);
	}
    }
    else
    {
      for(IntegerVectorList::iterator i=inequalities->begin();i!=inequalities->end();index++)
	{
	  int i2=newpos[index+1];
	  if(i2)
	    {
	      if(set_member(i2,A->linset))
		newLin.push_back(*i);
	      else
		newIn.push_back(*i);
	    }
	  i++;
	}

      for(IntegerVectorList::iterator i=equalities->begin();i!=equalities->end();index++)
	{
	  int i2=newpos[index+1];
	  if(i2)
	    {
	      if(set_member(i2,A->linset))
		newLin.push_back(*i);
	      else
		newIn.push_back(*i);
	    }
	  i++;
	}
      assert(index==numberOfRows);
    }

  assert(set_card(A->linset)==newLin.size());

  if(A->rowsize!=newLin.size()+newIn.size())
    {
      fprintf(stderr,"A->rowsize: %i\n",(int)A->rowsize);
      fprintf(stderr,"newLin.size(): %i\n",newLin.size());
      fprintf(stderr,"newIn.size(): %i\n",newIn.size());

      dd_WriteMatrix(Stderr,A);

      AsciiPrinter(Stderr).printVectorList(newLin);
      AsciiPrinter(Stderr).printVectorList(newIn);

    }
  assert(A->rowsize==newLin.size()+newIn.size());


  set_free(impl_linset);
  if(removeInequalityRedundancies)
    set_free(redset);
  free(newpos);
  //  set_free();//HOW DO WE FREE newpos?

  *equalities=newLin;
  *inequalities=newIn;

  dd_FreeMatrix(A);
  return;
 _L99:
  assert(0);
}

static dd_MatrixPtr vectorLists2MatrixGmp(int n, IntegerVectorList const &inequalities, IntegerVectorList const &equations, dd_ErrorType *err)
{
  IntegerVectorList g=inequalities;
  int numberOfInequalities=inequalities.size();
  for(IntegerVectorList::const_iterator i=equations.begin();i!=equations.end();i++)
    g.push_back(*i);

  //  assert(g.size());  // If this restriction turns out to be a problem it should be fixed in vectorList2MatrixGmp()
  int numberOfRows=g.size();

  dd_MatrixPtr A=NULL;

  cddinitGmp();

  A=vectorList2MatrixGmp(n, g, err);

  for(int i=numberOfInequalities;i<numberOfRows;i++)
    set_addelem(A->linset,i+1);
  return A;
}


static IntegerVectorList getConstraints(dd_MatrixPtr A, bool returnEquations)
{
  IntegerVectorList ret;

  int rowsize=A->rowsize;
  int n=A->colsize-1;
  // fprintf(Stderr,"rowsize: %i ColSize:%i\n",rowsize,n);

  mpq_t *point = new mpq_t [n];
  for(int j=0;j<n;j++)mpq_init(point[j]);

  for(int i=0;i<rowsize;i++)
    {
      bool isEquation=set_member(i+1,A->linset);
      if(isEquation==returnEquations)
	{
	  for(int j=0;j<n;j++)mpq_set(point[j],A->matrix[i][j+1]);

	  scaleToIntegerVector(point,n);
	  IntegerVector v=arrayToIntegerVector(point, n);
      //  AsciiPrinter(Stderr).printVector(v);

	  ret.push_back(v);
	}
    }

  for(int j=0;j<n;j++)mpq_clear(point[j]);
  delete [] point;

  return ret;
}


void LpSolverCddGmp::dual(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations, IntegerVectorList *dualInequalities, IntegerVectorList *dualEquations)
{
  IntegerVectorList dummy1;
  IntegerVectorList dummy2;
  if(!dualInequalities)dualInequalities=&dummy1;
  if(!dualEquations)dualEquations=&dummy2;

  int result;

  dd_MatrixPtr A=NULL;
  dd_ErrorType err=dd_NoError;

  cddinitGmp();

  A=vectorLists2MatrixGmp(n, inequalities, equations, &err);
  //  A=vectorList2MatrixGmp(g, &err);

  //    dd_WriteMatrix(Stderr,A);

  dd_PolyhedraPtr poly;
  poly=dd_DDMatrix2Poly2(A, dd_LexMin, &err);

  //  dd_WriteIncidence(Stderr,poly);
  if (poly->child==NULL || poly->child->CompStatus!=dd_AllFound) assert(0);

  //  dd_WriteIncidence(Stderr,poly);
  //  dd_WritePolyFile(Stderr,poly);
  //  dd_WriteMatrix(Stderr,poly->child->A);
  dd_MatrixPtr      A2=dd_CopyGenerators(poly);
  //  dd_WriteMatrix(Stderr,A2);

  *dualInequalities=getConstraints(A2,false);
  *dualEquations=getConstraints(A2,true);

  dd_FreeMatrix(A2);

  //  AsciiPrinter(Stderr).printVectorList(*dualInequalities);
  //  AsciiPrinter(Stderr).printVectorList(*dualEquations);


  //  assert(0);

  //  AsciiPrinter(Stderr).printVectorList(ret);
  //  dd_WriteIncidence(Stderr,poly);

  dd_FreeMatrix(A);
  dd_FreePolyhedra(poly);

  return;
 _L99:
  assert(0);
}

static LpSolverCddGmp theLpSolverCddGmp;
