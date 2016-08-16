#ifndef LINALGFLOAT_H_INCLUDED
#define LINALGFLOAT_H_INCLUDED

#include <assert.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <stdlib.h>

#define EPSILON 0.0001

#define typ double

namespace linalgfloat
{
class Vector{
    int n;
    std::vector<typ> data;
  public:
  Vector(int n_=0):
    n(n_),
      data(n_)
    {
      for(int i=0;i<n;i++)data[i]=0;
    }
    inline int size()const
    {
      return n;
    }
  inline typ& operator[](int n)
    {
      if(!(n>=0 && n<data.size()))assert(!"outOfRange(n,v.size())");
      return (data[n]);
    }
  inline typ const& operator[](int n)const
    {
      if(!(n>=0 && n<data.size()))assert(!"outOfRange(n,v.size())");
      return (data[n]);
    }
  friend Vector concatenation(Vector const &a, Vector const &b)
  {
    Vector ret(a.size()+b.size());
    for(int i=0;i<a.size();i++)ret[i]=a[i];
    for(int i=0;i<b.size();i++)ret[i+a.size()]=b[i];
    return ret;
  }
  friend inline typ dot(Vector const &a, Vector const &b)
  {
    assert(a.size()==b.size());
    typ ret=0;
    for(int i=0;i<a.size();i++)ret+=a[i]*b[i];
    return ret;
  }
  friend Vector operator*(typ s, const Vector& q)
  {
    Vector ret(q.size());
    for(int i=0;i<ret.size();i++)ret[i]=s*q[i];
    return ret;
  }
  Vector& operator+=(const Vector& q)
    {
      assert(q.size()==size());

      for(int i=0;i<size();i++)data[i]+=q[i];
      return *this;
    }
  Vector& operator-=(const Vector& q)
    {
      assert(q.size()==size());

      for(int i=0;i<size();i++)data[i]-=q[i];
      return *this;
    }
  friend std::ostream& operator<<(std::ostream& s, const Vector &v);
  Vector subvector(int begin, int end)const
    {
      assert(begin>=0);
      assert(end<=size());
      assert(end>=begin);
      Vector ret(end-begin);
      for(int i=0;i<end-begin;i++)
	ret[i]=data[begin+i];
      return ret;
    }
  Vector& madd(typ s, const Vector& q)
    {
      assert(q.size()==size());
      for(int i=0;i<q.size();i++)
	data[i]+=s*q.data[i];
      return *this;
    }
  bool isZero()const
  {
    for(int i=0;i<size();i++)
      if((data[i]>EPSILON)||(data[i]<-EPSILON))return false;
    return true;
  }

  friend class Matrix;
  };

class Matrix{
  int height,width;
  typ *data;
 public:
  static inline bool isZero(typ a){return (a<EPSILON)&&(a>-EPSILON);}
  class RowRef{
    int rowNum;
    Matrix &matrix;
  public:
  inline RowRef(Matrix &matrix_, int rowNum_):
    rowNum(rowNum_),
      matrix(matrix_)
      {
      }
    inline typ &operator[](int j)
      {
	assert(j>=0);
	assert(j<matrix.width);
	return matrix.data[matrix.width*rowNum+j];
      }
    Vector toVector()
    {
      Vector ret(matrix.width);
      for(int j=0;j<matrix.width;j++)
	ret[j]=matrix.data[matrix.width*rowNum+j];
      return ret;
    }
    void set(Vector const &v)
    {
      assert(v.size()==matrix.width);
      for(int j=0;j<matrix.width;j++)
	matrix.data[matrix.width*rowNum+j]=v[j];
    }
    bool isZero()const
    {
      for(int j=0;j<matrix.width;j++)if(!matrix.isZero(matrix.data[matrix.width*rowNum+j]))return false;
      return true;
    }
  };
  class const_RowRef{
    int rowNum;
    Matrix const &matrix;
  public:
  inline const_RowRef(const Matrix  &matrix_, int rowNum_):
    rowNum(rowNum_),
      matrix(matrix_)
      {
      }
    inline typ const &operator[](int j)
      {
	assert(j>=0);
	assert(j<matrix.width);
	return matrix.data[matrix.width*rowNum+j];
      }
    Vector toVector()
    {
      Vector ret(matrix.width);
      for(int j=0;j<matrix.width;j++)
	ret[j]=matrix.data[matrix.width*rowNum+j];
      return ret;
    }
    bool isZero()const
    {
      for(int j=0;j<matrix.width;j++)if(!matrix.isZero(matrix.data[matrix.width*rowNum+j]))return false;
      return true;
    }
  };
 Matrix():
  height(0),
    width(0)
  {
    data=0;
  }
 Matrix(int height_, int width_):
  height(height_),
    width(width_)
  {
    data=(typ*)malloc(height*width*sizeof(typ));
    //    for(int i=0;i<height*width;i++)data[i]=0;
    const int I=height*width;
    for(int i=0;i<I;i++)data[i]=0;
    //memset(data,0,sizeof(typ)*width*height);
  }
 Matrix(const Matrix &m):
  height(m.height),
    width(m.width)
    {
     data=(typ*)malloc(height*width*sizeof(typ));
     //for(int i=0;i<height*width;i++)data[i]=m.data[i];
     const int I=height*width;
     for(int i=0;i<I;i++)data[i]=m.data[i];
     //memcpy(data,m.data,sizeof(typ)*width*height);
    }
  Matrix& operator=(const Matrix& m)
    {
      if(this==&m)
	{
	  return *this;
	}
      else
	{
	  width=m.width;
	  height=m.height;
	  if(data)free(data);
	  data=(typ*)malloc(height*width*sizeof(typ));
	  //for(int i=0;i<height*width;i++)data[i]=m.data[i];
	  const int I=height*width;
	  for(int i=0;i<I;i++)data[i]=m.data[i];
	  //memcpy(data,m.data,sizeof(typ)*width*height);
	  return *this;
	}
    }
  ~Matrix()
    {
      free(data);
      data=0;
    }
  static Matrix identity(int n);
  /**
     Returns the number of rows of the matrix.
  */
  inline int getHeight()const
  {
    return height;
  }
  /**
     Returns the number of columns of the matrix.
  */
  inline int getWidth()const
  {
    return width;
  }
  inline RowRef operator[](int i)
  {
    assert(i>=0);
    assert(i<height);
    return RowRef(*this,i);
  }
  inline const_RowRef operator[](int i)const//should really return const_RowRef
  {
    assert(i>=0);
    assert(i<height);
    return const_RowRef(*this,i);
  }
  void multiplyAndAddColumn(int sourceColumn, typ scalar, int destinationColumn)
  {
    // for(int i=0;i<height;i++)
	  //  (*this)[i][destinationColumn]+=scalar* (*this)[i][sourceColumn];
      int width=this->width;int height=this->height;
          for(int i=0;i<height;i++)
      {
    	data[i*width+destinationColumn]+=scalar*data[i*width+sourceColumn];
	}
	  /*
    typ * __restrict  source=data+sourceColumn;
    typ * __restrict  dest=data+destinationColumn;
    for(int i=0;i<height;i++)
      {
    	(*(dest))+=scalar*(*(source));
	dest+=width;
	source+=width;
	}*/
  }
  void multiplyAndAddRow(int sourceRow, typ scalar, int destinationRow)
  {
    typ * __restrict source=data+width*sourceRow;
    typ * __restrict dest=data+width*destinationRow;
    for(int i=0;i<width;i++)
      dest[i]+=scalar*source[i];
	//    for(int i=0;i<width;i++)
	// (*this)[destinationRow][i]+=scalar* (*this)[sourceRow][i];
  }
  void scaleColumn(int column, typ scalar)
  {
    for(int i=0;i<height;i++)
      (*this)[i][column]*=scalar;
  }
  void scaleRow(int row, typ scalar)
  {
    for(int i=0;i<width;i++)
      (*this)[row][i]*=scalar;
  }

  friend std::ostream& operator<<(std::ostream& s, const Matrix &m);
  friend Matrix combineOnTop(Matrix const &top, Matrix const &bottom);
  Matrix transposed()const;
  Matrix operator-()const;
  Matrix operator*(Matrix const &b)const{assert(width==b.height);Matrix ret(height,b.width);for(int i=0;i<ret.height;i++)for(int j=0;j<ret.width;j++){typ s=0;for(int k=0;k<width;k++)s+=(*this)[i][k]*b[k][j];ret[i][j]=s;}return ret;}

  /**
   * Computes the dot product of the ith and jth row.
   */
  inline typ rowDot(int i, int j)
  {
	  assert(i>=0);
	  assert(i<getHeight());
	  assert(j>=0);
	  assert(j<getHeight());

	  typ ret=0;
	  typ *src1=data+width*i;
	  typ *src2=data+width*j;
	  for(int k=0;k<width;k++)
	  {
		  ret+=src1[k]*src2[k];
	  }
	  return ret;
  }

  /**
   * Computes the dot product of the ith and jth row.
   */
  inline typ rowDot(int i, Vector const &v)
  {
	  assert(i>=0);
	  assert(i<getHeight());
	  assert(width==v.size());

	  typ ret=0;
	  typ *src1=data+width*i;
	  for(int k=0;k<width;k++)
	  {
		  ret+=src1[k]*v[k];
	  }
	  return ret;
  }

  /*
    Computes the dot product of ith column of *this with jth row of m.
   */
  inline typ rowDotColumnOfOther(int i, Matrix const &m, int j)
  {
    assert(i>=0);
    assert(i<height);
    assert(j>=0);
    assert(j<m.width);
    assert(width==m.height);

    typ ret=0;
    typ *src1=data+width*i;
    typ *src2=m.data+j;
    for(int k=0;k<width;k++)
      {
	ret+=src1[k]*(*src2);
	src2+=m.width;
      }
    return ret;
  }

  int numberOfPivots()const;
  int reduceAndComputeRank();
  Matrix reduceAndComputeKernel();
  Vector normalForm(Vector v)const;//assume reduced
  Matrix normalForms(Matrix const &m)const;//assume reduced
  inline void swapRows(int a, int b)
  {for(int j=0;j<getWidth();j++){typ temp=(*this)[a][j];(*this)[a][j]=(*this)[b][j];(*this)[b][j]=temp;}}
  /*{
    assert(a>=0);
    assert(b>=0);
    assert(a<height);
    assert(b<height);
    typ *__restrict aRow=data+a*width;
    typ *__restrict bRow=data+b*width;
    for(int j=0;j<width;j++){typ temp=aRow[j];aRow[j]=bRow[j];bRow[j]=temp;}
    }*/
  int findRowIndex(int column, int currentRow)const;
  inline bool nextPivot(int &i, int &j)const//;
  //bool Matrix::nextPivot(int &i, int &j)const//iterates through the pivots in a matrix in reduced row echelon form. To find the first pivot put i=-1 and j=-1 and call this routine. When no more pivots are found the routine returns false.
{
  i++;
  if(i>=height)return false;
  while(++j<width)
    {
      if(!isZero((*this)[i][j])) return true;
    }
  return false;
}
  int reduce(bool returnIfZeroDeterminant);
  typ reduceAndComputeDeterminant();
  void cycleColumnsLeft(int offset);

  void REformToRREform(bool scalePivotsToOne=false);
  Matrix submatrix(int startRow, int startColumn, int endRow, int endColumn)const
  {
    assert(startRow>=0);
    assert(startColumn>=0);
    assert(endRow>=startRow);
    assert(endColumn>=startColumn);
    assert(endRow<=height);
    assert(endColumn<=width);
    Matrix ret(endRow-startRow,endColumn-startColumn);
    for(int i=startRow;i<endRow;i++)
      for(int j=startColumn;j<endColumn;j++)
	ret[i-startRow][j-startColumn]=(*this)[i][j];
    return ret;
  }
  void removeZeroRows();
  Matrix inverse()const;
  /*
    Returns a matrix with those columns removed, that would not conain a pivot in a row Echelon form.
   */
  Matrix reduceDimension()const;
  void maddRowToVector(int row, typ scalar, Vector &v)
  {
    assert(width==v.size());
    int offset=width*row;
    for(int i=0;i<width;i++)v.data[i]+=scalar*data[offset++];
    //    v.madd(scalar,(*this)[row].toVector());
  }
#if 0
  inline unsigned char hashValue(int row, int numberOfEntriesToConsider)const
  {
    unsigned char ret=0;
    typ *d=data+row*width;
    for(int i=0;i<numberOfEntriesToConsider;i++)
      ret+=((unsigned char*)(d+i))[6];
    return ret;
  }
#else
  inline unsigned char hashValue(int row, int numberOfEntriesToConsider)const
  {
    int ret=0;
    typ *d=data+row*width;
    for(int i=0;i<numberOfEntriesToConsider;i++)
      ret=((int*)(d+i))[1]+(ret>>7)+(ret<<25);
    return ret+(ret>>24)+(ret>>16)+(ret>>8);
  }
#endif
  bool rowsAreEqual(int row1, int row2, int numberOfEntriesToConsider)
  {
    typ *r1=data+row1*width;
    typ *r2=data+row2*width;
    for(int i=0;i<numberOfEntriesToConsider;i++)
      {
	if(!isZero(r1[i]-r2[i]))return false;
      }
    return true;
  }
  /*
   * This method transforms the set of row vectors into and orthogonal basis.
   * It returns a matrix describing the change
   */
  void orthogonalize()
  {
	  assert(0);
/*	  Matrix ret(height,height);
	  int retIndex=0;
	  list<int> toCheck;
	  for(int i=0;i<getHeight();i++)
	  {
		  Vector coef(height);
		  for(list<int>::const_iterator j=toCheck.begin();j!=toCheck.end();j++)
			  this->multiplyAndAddRow(*j,-rowDot(i,*j)/rowDot(*j,*j),i);//TODO: rowDot(j,j) is computed repeatedly - fix this
		  if(!(*this)[i].isZero())
			  {

				  retIndex++;
			  toCheck.push_back(i);
			  }
	  }
	  removeZeroRows();
*/
	  }
  /*
   * This method assumes that rows of the matrix are orthonogal and computes
   * the coefficients of the projection of v onto this basis.
   */
  Vector projectionCoefficients(Vector const &v)
  {
	  Vector ret(height);
	  for(int i=0;i<height;i++)ret[i]=rowDot(i,v)-rowDot(i,i);
	  return ret;
  }
};
};
#endif
