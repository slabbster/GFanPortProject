#include "field.h"
#include "field_rationals.h"
#include <vector>
#include "vektor.h"
#include "linalg.h"

#include "printer.h"
using namespace std;

static void swapRow(vector<vector<FieldElement> > &m, int i, int j)
{
  assert(i!=j);

  vector<FieldElement> temp;
  temp=m[i];
  m[i]=m[j];
  m[j]=temp;
}

static void madd(vector<vector<FieldElement> > &m, int i, FieldElement a, int j, int width)
{
  assert(i!=j);

  //  if(a.isZero())return;
  for(int k=0;k<width;k++)
    m[j][k]=m[j][k]+m[i][k]*a;
}

static int findRowIndex(vector<vector<FieldElement> > &m, int column, int currentRow)
{
  for(int i=currentRow;i<m.size();i++)
    if(!m[i][column].isZero())return i;
  return -1;
}

static void printMatrix(vector<vector<FieldElement> > const &m, int height, int width)
{
  AsciiPrinter P(Stderr);
  for(int i=0;i<height;i++)
    {
      for(int j=0;j<width;j++)
	{
	  P.printFieldElement(m[i][j]);
	  P.printString(" ");
	}
      P.printNewLine();
    }
  P.printNewLine();
}

static int reduce(vector<vector<FieldElement> > &m, int height, int width, bool returnIfZeroDeterminant=false)
{
  int retSwaps=0;
  int currentRow=0;
  //  fprintf(Stderr,"Reducing\n");
  for(int i=0;i<width;i++)
    {
      //      printMatrix(m,height,width);

      int s=findRowIndex(m,i,currentRow);

      //      fprintf(Stderr,"rowIndex:%i currentRow: %i Column: %i\n",s,currentRow,i);
      if(s!=-1)
	{
	  if(s!=currentRow)
	    {
	      swapRow(m,currentRow,s);
	      retSwaps++;
	      
	      //	      fprintf(Stderr,"SWAP:\n");
	      //   printMatrix(m,height,width);
	    }
	  for(int j=currentRow+1;j<height;j++)
	    madd(m,currentRow,-m[j][i]*m[currentRow][i].inverse(),j,width);
	  currentRow++;
	}
      else
	if(returnIfZeroDeterminant)return -1;
    }
  //  fprintf(Stderr,"Done reducing\n");
  return retSwaps;
}

int determinantSign(IntegerVectorList const &l)
{
  /*  static Field* field;

  if(!field)field=Field::find("GmpRationals"); // this is a bit stupid. We should add a field_rationals header file instead
  */

  /*  int n=l.size();
  assert(n>0);
  //  if(n==0)return 1;
  vector<vector<FieldElement> >m;

  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      assert(i->size()==n);
      vector<FieldElement> v(n);
      for(int j=0;j<n;j++)
	v[j]=Q.zHomomorphism((*i)[j]);
      m.push_back(v);
    }
  int swaps=reduce(m,n,n,true);

  if(swaps==-1)return 0;

  FieldElement d=Q.zHomomorphism(1);
  for(int i=0;i<n;i++)
    d=d*m[i][i];
  */
  FieldMatrix M=integerMatrixToFieldMatrix(rowsToIntegerMatrix(l),Q);
  FieldElement d=M.reduceAndComputeDeterminant();
  int swaps=0;

  string s=d.toString(true,true);//very very stupid
  if(d.isZero())return 0;
  if(s[0]=='+')
    {
      if(swaps&1)return -1;
      else return 1;
    }
  if(s[0]=='-')
    {
      if(swaps&1)return 1;
      else return -1;
    }
  return 0;
}
