#include "lp.h"

#include <assert.h>
#include <math.h>
#include <string>

#include "linalg.h"
#include "timer.h"

static Timer lpTimer("LP",100);
static Timer lpTimer2("LP2",100);
static Timer lpTimer3("LP3 - rank of matrix",100);

LpSolver *LpSolver::list;


LpSolver::LpSolver()
{
   next=list;
   list=this;
}


LpSolver *LpSolver::find(const char *name)
{
   LpSolver *l=list;
   while(l)
      {
         if(std::string(l->name())==std::string(name))break;
         l=l->next;
      }
   return l;
}


void LpSolver::printList(FILE *f)
{
   fprintf(f,"List of linked LP solvers:\n");
   LpSolver *l=list;
   while(l)
      {
         fprintf(f," %s\n",l->name());
         l=l->next;
      }
}


bool LpSolver::interiorPoint(const IntegerVectorList &g, IntegerVector &result, bool strictlyPositive, IntegerVector const *equalitySet)
{
  fprintf(stderr,"interiorPoint method not supported in \"%s\" LP class\n",name());
  assert(0);

  return false;
}


bool LpSolver::hasInteriorPoint(const IntegerVectorList &g, bool strictlyPositive, IntegerVector const *equalitySet)
{
  fprintf(stderr,"hasInteriorPoint method not supported in \"%s\" LP class\n",name());
  assert(0);
}


IntegerVectorList::const_iterator LpSolver::shoot(const IntegerVectorList &g)
{
  fprintf(stderr,"shoot method not supported in \"%s\" LP class\n",name());
  assert(0);
  return g.begin();
}



bool LpSolver::positiveVectorInKernel(const IntegerVectorList &g, IntegerVector *result)
{
  fprintf(stderr,"positiveVectorInKernel method not supported in \"%s\" LP class\n",name());
  assert(0);
  return false;
}


int LpSolver::rankOfMatrix(const IntegerVectorList &g)
{
  fprintf(stderr,"rankOfMatrix method not supported in \"%s\" LP class\n",name());
  assert(0);
  return 0;
}


IntegerVectorList LpSolver::extremeRaysInequalityIndices(const IntegerVectorList &inequalityList)
{
  fprintf(stderr,"extremeRaysInequalityIndices not supported in \"%s\" LP class\n",name());
  assert(0);
  return IntegerVectorList();
}


void LpSolver::removeRedundantRows(IntegerVectorList *inequalities, IntegerVectorList *equalities, bool removeInequalityRedundancies)
{
  fprintf(stderr,"removeRedundantRows method not supported in \"%s\" LP class\n",name());
  assert(0);
}

IntegerVector LpSolver::relativeInteriorPoint(int n, const IntegerVectorList &g, IntegerVector const *equalitySet)
{
  fprintf(stderr,"relativeInteriorPoint method not supported in \"%s\" LP class\n",name());
  assert(0);
  return IntegerVector();
}

void LpSolver::dual(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations, IntegerVectorList *dualInequalities, IntegerVectorList *dualEquations)
{
  fprintf(stderr,"dual method not supported in \"%s\" LP class\n",name());
  assert(0);
}


bool LpSolver::hasHomogeneousSolution(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations)
{
  fprintf(stderr,"hasHomogeneousSolution method not supported in \"%s\" LP class\n",name());
  assert(0);
}

static LpSolver *soplex,*soplexCddGmp,*huber,*cdd,*cddgmp,*default_;
static bool initialized;


bool lpSetSolver(const char *name)
{
  soplexCddGmp=LpSolver::find("SoPlexCddGmp");
  soplex=LpSolver::find("SoPlex");
  huber=LpSolver::find("Huber's");
  cdd=LpSolver::find("cdd");
  cddgmp=LpSolver::find("cddgmp");
  LpSolver *selected=LpSolver::find(name);
  default_=huber;
  if(soplex)default_=soplex;
  if(cdd)default_=cdd;
  if(cddgmp)default_=cddgmp;
  if(soplexCddGmp)default_=soplexCddGmp;
  if(selected)default_=selected;
  initialized=true;
  assert(default_);
  fprintf(stderr,"LP algorithm being used: \"%s\".\n",default_->name()); //TODO: change to debug printer
  //if(default_==soplexCddGmp)fprintf(stderr,"USING SoPlexCddGmp\n");

  return selected;
}


static void LpInit()
{
   if(!initialized)
      {
	lpSetSolver("");
      }
}


bool isFacet(const IntegerVectorList &g, IntegerVectorList::const_iterator i)
{
  TimerScope ts(&lpTimer);
  LpInit();

  return default_->isFacet(g,i);
}


bool interiorPoint(const IntegerVectorList &g, IntegerVector &result, bool strictlyPositive, IntegerVector const *equalitySet)
{
  LpInit();
  return default_->interiorPoint(g,result,strictlyPositive,equalitySet);
}


bool hasInteriorPoint(const IntegerVectorList &g, bool strictlyPositive, IntegerVector const *equalitySet)
{
  LpInit();
  return default_->hasInteriorPoint(g,strictlyPositive, equalitySet);
}


IntegerVectorList::const_iterator shootRay(const IntegerVectorList &g)
{
  TimerScope ts(&lpTimer2);
  LpInit();

  if(g.empty())return g.end();

  return default_->shoot(g);
}


bool positiveVectorInKernel(const IntegerVectorList &g, IntegerVector *result)
{
  LpInit();

  return default_->positiveVectorInKernel(g,result);
}


int rankOfMatrix(const IntegerVectorList &g)
{
  TimerScope ts(&lpTimer3);
  LpInit();

  return default_->rankOfMatrix(g);
}


IntegerVectorList extremeRaysInequalityIndices(const IntegerVectorList &inequalityList)
{
  /* If cone is simplicial, then the rays are easy to find...*/
  if(rankOfMatrix(inequalityList)==inequalityList.size())
    {
      IntegerVectorList ret;
      int m=inequalityList.size();
      for(int i=0;i<m;i++)
	{
	  IntegerVector v(m-1);
	  for(int j=0;j<i;j++)
	    v[j]=j;
	  for(int j=i+1;j<m;j++)
	    v[j-1]=j;
	  ret.push_back(v);
	}
      return ret;
    }


  LpInit();

  return default_->extremeRaysInequalityIndices(inequalityList);
}


void removeRedundantRows(IntegerVectorList *inequalities, IntegerVectorList *equalities, bool removeInequalityRedundancies)
{
  LpInit();

  return default_->removeRedundantRows(inequalities,equalities,removeInequalityRedundancies);
}


IntegerVector relativeInteriorPoint(int n, const IntegerVectorList &g, IntegerVector const *equalitySet)
{
  LpInit();

  return default_->relativeInteriorPoint(n,g,equalitySet);
}


void dual(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations, IntegerVectorList *dualInequalities, IntegerVectorList *dualEquations)
{
  LpInit();

  return default_->dual(n,inequalities,equations,dualInequalities,dualEquations);
}


bool hasHomogeneousSolution(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations)
{
  LpInit();

  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
    if(i->size()!=n)
      {
	AsciiPrinter(Stderr) << "Inequality length does not match. n="<<n<<" *i="<<*i<<"\n";
	assert(0);
      }
  for(IntegerVectorList::const_iterator i=equations.begin();i!=equations.end();i++)
    if(i->size()!=n)
      {
	AsciiPrinter(Stderr) << "Equation length does not match. n="<<n<<" *i="<<*i<<"\n";
	assert(0);
      }

  return default_->hasHomogeneousSolution(n, inequalities, equations);
}


bool isInNonNegativeSpan(IntegerVector const &v, IntegerVectorList const &rays, IntegerVectorList const &linealitySpace)
{
  int n=v.size();
  /*
  Solve Ax>=0 with A being:
    0| 1  | 000
    0|  1 | 000
    0|   1| 000
   -v| rrr| lll\
   -v| aaa| iii \ Added as
   -v| yyy| nnn / equations
   -v| sss| eee/
   */

  FieldMatrix A1=combineLeftRight(combineLeftRight(FieldMatrix(Q,rays.size(),1),FieldMatrix::identity(Q,rays.size())),FieldMatrix(Q,rays.size(),linealitySpace.size()));
  FieldMatrix temp=FieldMatrix(Q,1,n);
  temp[0]=integerVectorToFieldVector(-v,Q);
  FieldMatrix A2=combineLeftRight(combineLeftRight(temp.transposed(),integerMatrixToFieldMatrix(rowsToIntegerMatrix(rays,n),Q).transposed()),
				  integerMatrixToFieldMatrix(rowsToIntegerMatrix(linealitySpace,n),Q).transposed());

  return hasHomogeneousSolution(1+rays.size()+linealitySpace.size(),fieldMatrixToIntegerMatrixPrimitive(A1).getRows(),fieldMatrixToIntegerMatrixPrimitive(A2).getRows());
}
