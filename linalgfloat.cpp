//#include "vektor.h"

#include "linalgfloat.h"

#include <iostream>
#include <fstream>
#include "cstdio"


using namespace std;

namespace linalgfloat
{
std::ostream& operator<<(std::ostream& s, const Vector &v)
{
  s<<"(";
  for(int j=0;j<v.size();j++)
    {
      if(j)s<<",";
      double temp=v[j];
      s<<temp;
    }
  s<<")";
  return s;
}


  Matrix Matrix::identity(int n)
  {
    Matrix ret(n,n);
    for(int i=0;i<n;i++)ret[i][i]=1;
    return ret;
  }
  Matrix combineOnTop(Matrix const &top, Matrix const &bottom)
  {
    if(top.getWidth()!=bottom.getWidth())
      {
	cerr<<top<<bottom;
      }
    assert(top.getWidth()==bottom.getWidth());
    Matrix ret(top.getHeight()+bottom.getHeight(),top.getWidth());
    for(int j=0;j<top.getWidth();j++)
      {
	for(int i=0;i<top.getHeight();i++)ret[i][j]=top[i][j];
	for(int i=0;i<bottom.getHeight();i++)ret[i+top.getHeight()][j]=bottom[i][j];
      }
    return ret;
  }
  Matrix Matrix::transposed()const
  {
    Matrix ret(getWidth(),getHeight());
    for(int i=0;i<getHeight();i++)
      for(int j=0;j<getWidth();j++)
	ret[j][i]=(*this)[i][j];
    return ret;
  }
  Matrix Matrix::operator-()const
  {
    Matrix ret(*this);
    for(int i=0;i<getHeight()*getWidth();i++)ret.data[i]=-ret.data[i];
    return ret;
  }

  std::ostream& operator<<(std::ostream& s, const Matrix &m)
{
  s<<m.height<<" "<<m.width<<endl;
  for(int i=0;i<m.height;i++)
    {
      for(int j=0;j<m.width;j++)
	{
	  if(j)s<<" ";
	  double temp=m[i][j];
	  //	  s<<m[i][j];
	  s<<temp;
	}
      s<<" hash:"<<(int)m.hashValue(i,m.width);
      s<<endl;
    }
  return s;
}

int Matrix::findRowIndex(int column, int currentRow)const
{
  int best=-1;
  int bestNumberOfNonZero;
  for(int i=currentRow;i<height;i++)
    if(!isZero((*this)[i][column]))
      {
	return i;//<-------------------------------------------- no need to find row with many zeros
	int nz=0;
	for(int k=column+1;k<width;k++)
	  if(!isZero((*this)[i][k]))nz++;
	if(best==-1)
	  {
	    best=i;
	    bestNumberOfNonZero=nz;
	  }
	else if(nz<bestNumberOfNonZero)
	  {
	    best=i;
	    bestNumberOfNonZero=nz;
	  }
      }
  return best;
}



int Matrix::reduce(bool returnIfZeroDeterminant) //bool reducedRowEcholonForm,
/* Run a Gauss reduction. Returns the number of swaps. The number of
   swaps is need if one wants to compute the determinant
   afterwards. In this case it is also a good idea to set the flag
   which make the routine terminate when a it is discovered the the
   determinant is zero. */
{
  //  if(width<=1)cerr<<height<<"x"<<width<<endl;
  int retSwaps=0;
  int currentRow=0;

  for(int i=0;i<width;i++)
    {
      int s=findRowIndex(i,currentRow);

      if(s!=-1)
	{
	  if(s!=currentRow)
	    {
	      swapRows(currentRow,s);
	      retSwaps++;
	    }
	  for(int j=currentRow+1;j<height;j++)
	      {
		multiplyAndAddRow(currentRow,-(*this)[j][i]/(*this)[currentRow][i],j);
	      }
	  currentRow++;
	}
      else
	if(returnIfZeroDeterminant)return -1;
    }

  return retSwaps;
}


int Matrix::numberOfPivots()const
{
  int ret=0;
  int pivotI=-1;
  int pivotJ=-1;
  while(nextPivot(pivotI,pivotJ))ret++;
  return ret;
}

int Matrix::reduceAndComputeRank()
{
  reduce(false);
  return numberOfPivots();
}

typ Matrix::reduceAndComputeDeterminant()
{
  assert(height==width);
  int swaps=reduce(false);

  int r=reduceAndComputeRank();
  //cerr<<*this;
  if(r!=height)return 0;

  typ ret=(1-2*(swaps&1));

  int pivotI=-1;
  int pivotJ=-1;
  while(nextPivot(pivotI,pivotJ))ret=ret*(*this)[pivotI][pivotJ];
  return ret;
}


Matrix Matrix::reduceAndComputeKernel()
{
  Matrix ret(width-reduceAndComputeRank(),width);

  REformToRREform();

  int k=0;
  int pivotI=-1;
  int pivotJ=-1;
  bool pivotExists=nextPivot(pivotI,pivotJ);
  for(int j=0;j<width;j++)
    {
      if(pivotExists && (pivotJ==j))
	{
	  pivotExists=nextPivot(pivotI,pivotJ);
	  continue;
	}
      int pivot2I=-1;
      int pivot2J=-1;
      while(nextPivot(pivot2I,pivot2J))
	{
	  ret[k][pivot2J]=(*this)[pivot2I][j]/(*this)[pivot2I][pivot2J];
	}
      ret[k][j]=-1;
      k++;
    }
  return ret;
}


Vector Matrix::normalForm(Vector v)const
{
  int pivotI=-1;
  int pivotJ=-1;
  int nonpivots=v.size();
  while(nextPivot(pivotI,pivotJ))
    {
      nonpivots--;
      v-=(v[pivotJ]/(*this)[pivotI][pivotJ])*(*this)[pivotI].toVector();
    }
  Vector ret(nonpivots);
  pivotI=-1;
  pivotJ=-1;
  int i=0;
  int last=-1;
  while(nextPivot(pivotI,pivotJ))
    {
      while(pivotJ-1>last)
	{
	  ret[i++]=v[++last];
	  //	    cerr<<"("<<(i-1)<<","<<last<<")";
	}
      last=pivotJ;
    }
  last++;
  while(last<width)
    ret[i++]=v[last++];
//  if(debug)cerr<<v<<":"<<ret<<endl;
  assert(i==nonpivots);
  return ret;
}

Matrix Matrix::normalForms(Matrix const &m)const
{
  //cerr<<*this;
  Matrix ret(m.height,width-numberOfPivots());
  for(int i=0;i<m.height;i++)
    ret[i].set(normalForm(m[i].toVector()));
  return ret;
}
void Matrix::cycleColumnsLeft(int offset)
{
  for(int i=0;i<height;i++)
    {
      Vector temp=(*this)[i].toVector();
      for(int j=0;j<width;j++)(*this)[i][j]=temp[(j+offset+width)%width];
    }
}
void Matrix::removeZeroRows()
{
  int n=0;
  for(int i=0;i<height;i++)
    {
      bool isZer=true;
      for(int j=0;j<width;j++)if(!isZero((*this)[i][j]))isZer=false;
      if(!isZer)n++;
    }
  Matrix ret(n,width);
  n=0;
  for(int i=0;i<height;i++)
    {
      bool isZer=true;
      for(int j=0;j<width;j++)if(!isZero((*this)[i][j]))isZer=false;
      if(!isZer)ret[n++].set((*this)[i].toVector());
    }
  *this=ret;
}

void Matrix::REformToRREform(bool scalePivotsToOne)
{
  int pivotI=-1;
  int pivotJ=-1;
  while(nextPivot(pivotI,pivotJ))
    {
      if(scalePivotsToOne)
	(*this)[pivotI].set((1/(*this)[pivotI][pivotJ])*(*this)[pivotI].toVector());
      //  cerr<<*this;
      for(int i=0;i<pivotI;i++)
	multiplyAndAddRow(pivotI,-((*this)[i][pivotJ])/((*this)[pivotI][pivotJ]),i);
    }
}

Matrix Matrix::inverse()const
{
  //  cerr<<"THIS"<<*this;
  assert(height==width);
  Matrix temp=combineOnTop(transposed(),identity(height)).transposed();
  // cerr<<"TEMP"<<temp;
  temp.reduce(false);
  //cerr<<"TEMP"<<temp;
  temp.REformToRREform(true);
  //cerr<<"TEMP"<<temp;
  Matrix ret=temp.submatrix(0,height,height,2*height);
  ///cerr<<"RET"<<ret;
  //cerr<<"PROD"<<ret*(*this);
  return ret;
}

Matrix Matrix::reduceDimension()const
{
  Matrix temp=*this;
  temp.reduce(false);

  //  cerr<<"IN"<<*this;

  int d=temp.numberOfPivots();
  Matrix ret(height,d);

  int pivotI=-1;
  int pivotJ=-1;
  int i=0;
  while(temp.nextPivot(pivotI,pivotJ))
    {
      for(int k=0;k<height;k++)ret[k][i]=(*this)[k][pivotJ];
      i++;
    }
  //  cerr<<"OUT"<<ret;
  return ret;
}
}
