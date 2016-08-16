#include "lp_soplexcdd.h"

#include "printer.h"

#include "spxdefines.h"
#include "spxsolver.h"

#include "timer.h"
#include "spxpricer.h"
//#include "spxdefaultpr.h"
#include "spxparmultpr.h"
#include "spxdevexpr.h"
#include "spxhybridpr.h"
#include "spxsteeppr.h"
#include "spxweightpr.h"
#include "spxratiotester.h"
#include "spxharrisrt.h"
#include "spxdefaultrt.h"
#include "spxfastrt.h"
#include "spxsimplifier.h"
//#include "spxaggregatesm.h"
//#include "spxredundantsm.h"
//#include "spxrem1sm.h"
//#include "spxgeneralsm.h"
#include "spxscaler.h"
#include "spxequilisc.h"
#include "spxsumst.h"
#include "spxweightst.h"
#include "spxvectorst.h"
#include "slufactor.h"
#include "soplex.h"
#include "continuedfractions.h"
#include "matrix.h"

#include "log.h"
using namespace soplex;

/** Here comes a simple derived class from #SoPlex, which uses #terminate() as
 *  callback method for outputting statistics.
 */
class MySoPlex : public SoPlex
{
private:
   SLUFactor m_slu;

public:
   /// default constructor
   MySoPlex(SPxSolver::Type p_type = SPxSolver::LEAVE, SPxSolver::Representation p_rep = SPxSolver::COLUMN)
      : SoPlex(p_type, p_rep)
   {



   bool                   print_solution = false;
   bool                   print_quality  = false;
   NameSet                rownames;
   NameSet                colnames;
   SPxStarter*            starter        = 0;
   SPxSolver::Type           type           = SPxSolver::LEAVE;
   SPxSolver::Representation representation = SPxSolver::COLUMN;
   int                    precision;
   Real                   delta          = DEFAULT_BND_VIOL;
   Real                   epsilon        = DEFAULT_EPS_ZERO;
   int                    verbose        = 0;
   SLUFactor::UpdateType  update         = SLUFactor::FOREST_TOMLIN;
   Real                   timelimit      = 1.0;-1.0;
   SPxPricer*             pricer         = 0;
   SPxRatioTester*        ratiotester    = 0;
   SPxScaler*             scaler         = 0;
   SPxSimplifier*         simplifier     = 0;

   precision = int(-log10(delta)) + 1;

   Param::setEpsilon(epsilon);
   Param::setVerbose(verbose);



   //options -p4 -t2 -g1 -s0 -c0
   setUtype(update);
   setTerminationTime(timelimit);
   setDelta(delta);

   assert(isConsistent());

   pricer = new SPxSteepPR;
   setPricer(pricer);
   assert(isConsistent());

   ratiotester = new SPxFastRT;
   setTester(ratiotester);
   assert(isConsistent());

   /*   scaler = new SPxEquili(representation == SoPlex::COLUMN, true);
   setScaler(scaler);
   assert(isConsistent());
   */
   setSimplifier(simplifier);
   assert(isConsistent());

   setStarter(starter);
   assert(isConsistent());
   }

   virtual bool terminate()
   {
      /*      if (iteration() % 100 == 0)
         std::cout << iteration() << ":\t" << value() << std::endl;
      */
     return SoPlex::terminate();
   }

   void setUtype(SLUFactor::UpdateType tp)
   {
      m_slu.setUtype(tp);
   }

  void build(const IntegerVectorList &g, IntegerVectorList::const_iterator i)
   {
     int width=g.size()-1;
     int height=i->v.size();
     
     LPColSet cols(width,width*height);
     DSVector c1(height);
     
     //Build Cols
     for(IntegerVectorList::const_iterator j=g.begin();j!=g.end();j++)
       {
	 c1.clear();
	 for(int t=0;t<height;t++)
	   if(j->v[t]!=0)
	     c1.add(t,(double)(j->v[t]));
	 
	 if(j!=i)
	   {
	     Real obj=0;
	     Real upper=infinity;
	     Real lower=0;
	     
	     cols.add(obj,lower,c1,upper);
	   }
       }
     
     LPRowSet rows(height,width*height);
     DSVector r1(width);
     
     //Change Rows
     for(int t=0;t<height;t++)
       rows.add(i->v[t],r1,i->v[t]);
     
     addRows(rows);
     addCols(cols);
     
     changeSense(SPxLP::MINIMIZE);
     
     assert(isConsistent());
   }
  void buildType2(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations)
   {
     int width=n;
     int nInequalities=inequalities.size();
     int height=inequalities.size()+equations.size();

     IntegerMatrix m=rowsToIntegerMatrix(inequalities,n);
     m.append(rowsToIntegerMatrix(equations,n));
     
     LPColSet cols(width,width*height);
     DSVector c1(height);
     
     //Build Cols
     for(int i=0;i<width;i++)
       {
	 c1.clear();
	 for(int t=0;t<height;t++)
	   if(m[t][i]!=0)
	     c1.add(t,(double)(m[t][i]));

	 Real obj=0;
	 Real upper=infinity;
	 //Real lower=0;
	 Real lower=-infinity;//is this correct?

	 if(i==0)
	   {
	     upper=1;
	     lower=1;
	   }

	 cols.add(obj,lower,c1,upper);
       }
     
     LPRowSet rows(height,width*height);
     DSVector r1(width);
     
     //Change Rows
     for(int t=0;t<nInequalities;t++)
       rows.add(0,r1,infinity);
     for(int t=nInequalities;t<height;t++)
       rows.add(0,r1,0);
     
     addRows(rows);
     addCols(cols);
     
     changeSense(SPxLP::MINIMIZE);
     
     assert(isConsistent());
   }
};

static int toint(float r)
{
   return *((int*)&r);
}

static void printLP(SPxLP &w)
{
      std::cout << "LP has " 
             << w.nRows() 
             << "\trows and\n       "
             << w.nCols() 
             << "\tcolumns" 
             << std::endl;
      int nr=w.nRows();
      int nc=w.nCols();

      for(int i=0;i<nr;i++)
         {      
            for(int j=0;j<nc;j++)
               {
                  LPRow R;
                  w.getRow(i,R);
                  //  if(j<R.rowVector().size())
                     std::cout<<R.rowVector()[j]<<" ";
                  //  else
                     //                     std::cout<<(0.0)<<" ";
               }
            std::cout<<std::endl;
         }
      std::cout<<std::endl;

      for(int i=0;i<nr;i++)
         {      
            for(int j=0;j<nc;j++)
               {
                  LPCol C;
                  w.getCol(j,C);
                  //                  if(i<C.colVector().size())
                     std::cout<<C.colVector()[i]<<" ";
                     // else
                     // std::cout<<(0.0)<<" ";
               }
            std::cout<<std::endl;
         }

      std::cout<<"cols:"<<std::endl;

      for(int j=0;j<nc;j++)
         {
            LPCol C;
            w.getCol(j,C);
            std::cout<<C.lower()<<" "<<C.upper()<<" "<<C.obj()<<std::endl;
            std::cout<<toint(C.lower())<<" "<<toint(C.upper())<<" "<<toint(C.obj())<<std::endl;
         }

      std::cout<<"rows:"<<std::endl;

      for(int i=0;i<nr;i++)
         {
            LPRow R;
            w.getRow(i,R);
            std::cout<<toint(R.lhs())<<" "<<toint(R.rhs())<<" "<<R.type()<<std::endl;
         }
}


MySoPlex work(SPxSolver::LEAVE, SPxSolver::COLUMN);


static bool isFeasibleSolution(IntegerVector const &solution, int denominator, IntegerVectorList const &g, IntegerVectorList::const_iterator i)
{
  if(denominator<=0)return false;
  // To do: Truncate  
  IntegerVector sum=denominator*(*i);

  int t=0;
  for(IntegerVectorList::const_iterator j=g.begin();j!=g.end();j++)
    if(j!=i)
      {
	if(solution[t]<0)return false;
	sum-=solution[t]*(*j);
	t++;
      }
  return sum.isZero();
}

static bool isInfeasibilityCertificate(IntegerVector const &certificate, IntegerVectorList const &g, IntegerVectorList::const_iterator i)
{
  IntegerVector c=certificate;
  // To do: add truncation on c
  if(dotLong(c,*i)<=0)return false;

  for(IntegerVectorList::const_iterator j=g.begin();j!=g.end();j++)
    if(i!=j)
      if(dotLong(c,*j)>0)return false;
  return true;
}
/*	     for(IntegerVectorList::const_iterator j=g.begin();j!=g.end();j++)
	       {
		 /*		 double prod=0;
		 for(int i=0;i<work.nRows();i++)
		   {prod+=(*j)[i]*certificate[i];
		     //		 fprintf(stderr,"%f \n",prod);
		   }
		 int num,den;doubleToFraction(prod,num,den);
		 */
		 //		 fprintf(stderr,":%f: %i/%i\n",prod,num,den);
/*		 fprintf(stderr,":%i\n",dotLong(c,*j));
	       }
	*/       


static bool isFeasibleSolutionType2(IntegerVector const &s, IntegerVectorList const &inequalities, IntegerVectorList const &equations)
{
  // To do: do truncation
  if(s[0]<=0)return false;
  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
    if(dotLong(s,*i)<0)return false;
  for(IntegerVectorList::const_iterator i=equations.begin();i!=equations.end();i++)
    if(dotLong(s,*i)!=0)return false;
  return true;
}


static bool isInfeasibilityCertificateType2(IntegerVector const &c, int n, IntegerVectorList const &inequalities, IntegerVectorList const &equations)
{
  // To do: truncation
  int nInequalities=inequalities.size();
  if(!c.subvector(0,nInequalities).isNonNegative())return false;

  IntegerVector sum(n);
  int j=0;
  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++,j++)
    sum+=c[j]* *i;
  for(IntegerVectorList::const_iterator i=equations.begin();i!=equations.end();i++,j++)
    sum+=c[j]* *i;

  if(sum[0]>=0)return false;

  for(int i=1;i<n;i++)
    if(sum[i]<0)
      return false;

  return true;
}


bool LpSolverSoPlexCddGmp::isFacet(const IntegerVectorList &g, IntegerVectorList::const_iterator I)
{
   SPxSolver::Type           type           = SPxSolver::LEAVE;
   SPxSolver::Representation representation = SPxSolver::COLUMN;

   int lp_status=0;

   work.clear();

   work.build(g,I);

 retry:

   //   std::cerr<< work;
   work.solve();

   SPxSolver::Status stat = work.status();

   switch (stat)
     {
     case SPxSolver::OPTIMAL:
       {
         DVector objx(work.nCols());
         
         if( work.getPrimal(objx) != SPxSolver::ERROR )
	   {
	     vector<double> solution(work.nCols());
	     for(int i=0;i<work.nCols();i++)
	       solution[i]=objx[i];
	     
	     vector<int> solutionNum;
	     int denominator;
	     doubleVectorToFractions(solution,solutionNum,denominator);
	     IntegerVector s(solution.size());
	     for(int i=0;i<s.size();i++)s[i]=solutionNum[i];
	     
	     if(isFeasibleSolution(s,denominator,g,I))
	       {
		 log3 fprintf(Stderr,"Solution OK.\n");
		 return false;
	       }
	     log3 fprintf(Stderr,"Solution failed .\n");
	     goto fallBack;
	   }
       }
       break;
     case SPxSolver::UNBOUNDED:
       std::cerr << "LP is unbounded";
       lp_status=1;
       break;
     case SPxSolver::INFEASIBLE:
       {
	 DVector farkasx(work.nRows());
	 
	 if( work.getDualfarkas(farkasx) != SPxSolver::ERROR )
	   {
	     vector<double> certificate(work.nRows());
	     for(int i=0;i<work.nRows();i++)
	       certificate[i]=farkasx[i];
	     
	     vector<int> certificateNum;
	     int denominator;
	     doubleVectorToFractions(certificate,certificateNum,denominator);
	     IntegerVector c(certificate.size());
	     for(int i=0;i<c.size();i++)c[i]=certificateNum[i];
	     
	     if(isInfeasibilityCertificate(c,g,I))
	       {
		 log3 fprintf(Stderr,"Certificate for infeasibility OK.\n");
		 return true;
	       }
	     log3 fprintf(Stderr,"Certificate failed.\n");
	   }
	 else
	   {
	     log3 fprintf(Stderr,"Error while producing certificate\n");
	   }
	 goto fallBack;
       }
       break;
     case SPxSolver::ABORT_TIME:
       std::cout << "aborted due to time limit";
       lp_status=1;
       break;
     case SPxSolver::ABORT_ITER:
       std::cout << "aborted due to iteration limit";
       lp_status=1;
       break;
     case SPxSolver::ABORT_VALUE:
       std::cout << "aborted due to objective value limit";
       lp_status=1;
       break;
     default:
       std::cout << "An error occurred during the solution process";
       lp_status=1;
       break;
     }

 fallBack:
   log0 fprintf(Stderr,"Falling back on CddLib\n");
   return LpSolverCddGmp::isFacet(g,I);
}



bool LpSolverSoPlexCddGmp::hasHomogeneousSolution(int n, const IntegerVectorList &inequalities, const IntegerVectorList &equations)
{
   SPxSolver::Type           type           = SPxSolver::LEAVE;
   SPxSolver::Representation representation = SPxSolver::COLUMN;

   int lp_status=0;

   work.clear();

   work.buildType2(n,inequalities,equations);

 retry:

   //   std::cerr<< work;

   //   assert(0);
   work.solve();

   SPxSolver::Status stat = work.status();

   switch (stat)
     {
     case SPxSolver::OPTIMAL:
       {
         DVector objx(work.nCols());
         
         if( work.getPrimal(objx) != SPxSolver::ERROR )
	   {
	     vector<double> solution(work.nCols());
	     for(int i=0;i<work.nCols();i++)
	       solution[i]=objx[i];
	     
	     vector<int> solutionNum;
	     int denominator;
	     doubleVectorToFractions(solution,solutionNum,denominator);
	     IntegerVector s(solution.size());
	     for(int i=0;i<s.size();i++)s[i]=solutionNum[i];
	     
	     //  AsciiPrinter(Stderr).printVector(s);
	     if(isFeasibleSolutionType2(s,inequalities,equations))
	       {
		 log3 fprintf(Stderr,"Solution OK.\n");
		 return true;
	       }
	     log2 fprintf(Stderr,"Solution failed (Type2).\n");

	     /*	     for(int i=0;i<work.nCols();i++)
	       {
		 std::cerr<<solution[i]<<',';
	       }
	     std::cerr<< work;
	     AsciiPrinter(Stderr).printVector(s);
	     AsciiPrinter(Stderr).printVectorList(inequalities);
	     AsciiPrinter(Stderr).printVectorList(equations);
	     assert(0);
	     */
	     goto fallBack;
	   }
       }
       break;
     case SPxSolver::UNBOUNDED:
       std::cerr << "LP is unbounded";
       lp_status=1;
       break;
     case SPxSolver::INFEASIBLE:
       {
	 DVector farkasx(work.nRows());
	 
	 if( work.getDualfarkas(farkasx) != SPxSolver::ERROR )
	   {
	     vector<double> certificate(work.nRows());
	     for(int i=0;i<work.nRows();i++)
	       certificate[i]=farkasx[i];
	     
	     vector<int> certificateNum;
	     int denominator;
	     doubleVectorToFractions(certificate,certificateNum,denominator);
	     IntegerVector c(certificate.size());
	     for(int i=0;i<c.size();i++)c[i]=certificateNum[i];
	     
	     if(isInfeasibilityCertificateType2(c,n,inequalities,equations))
	       {
		 log3 fprintf(Stderr,"Certificate for infeasibility OK.\n");
		 return false;
	       }
	     
	       log2 fprintf(Stderr,"Certificate failed (Type2).\n");
	       /* std::cerr<< work;
	      std::cerr<< farkasx;
	     AsciiPrinter(Stderr).printVector(c);
	     AsciiPrinter(Stderr).printVectorList(inequalities);
	     AsciiPrinter(Stderr).printVectorList(equations);
	     assert(0);
	     */
	   }
	 else
	   {
	     log3 fprintf(Stderr,"Error while producing certificate\n");
	   }
	 goto fallBack;
	 }
       goto fallBack;
       break;
     case SPxSolver::ABORT_TIME:
       std::cout << "aborted due to time limit";
       lp_status=1;
       break;
     case SPxSolver::ABORT_ITER:
       std::cout << "aborted due to iteration limit";
       lp_status=1;
       break;
     case SPxSolver::ABORT_VALUE:
       std::cout << "aborted due to objective value limit";
       lp_status=1;
       break;
     default:
       std::cout << "An error occurred during the solution process";
       lp_status=1;
       break;
     }

 fallBack:
   log0 fprintf(Stderr,"Falling back on CddLib\n");
   return LpSolverCddGmp::hasHomogeneousSolution(n, inequalities,equations);
}


static LpSolverSoPlexCddGmp theLpSolverSoPlexCdd;
