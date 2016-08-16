#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <vector>
#include <list>
#include <assert.h>
#include <algorithm>
#include "vektor.h"

using namespace std;

template <class typ> class Matrix{
  //public:
  int width,height;
  vector<Vektor<typ> > rows;
public:
  inline int getHeight()const{return height;};
  inline int getWidth()const{return width;};
  Matrix(const Matrix &a):rows(a.rows),width(a.getWidth()),height(a.getHeight()){
  }
  Matrix(int height_, int width_):rows(height_),height(height_),width(width_){
    assert(height>=0);
    assert(width>=0);
    for(int i=0;i<getHeight();i++)rows[i]=Vektor<typ>(width);
  };
  ~Matrix(){
  };
  Matrix():width(0),height(0){
  };
  Vektor<typ> column(int i)const
    {
      assert(i>=0);
      assert(i<getWidth());
      Vektor<typ> ret(getHeight());
      for(int j=0;j<getHeight();j++)ret[j]=rows[j][i];
      return ret;
    }
  Matrix transposed()const
    {
      Matrix ret(getWidth(),getHeight());
      for(int i=0;i<getWidth();i++)
	ret.rows[i]=column(i);
      return ret;
    }
  static Matrix identity(int n)
    {
      Matrix m(n,n);
      for(int i=0;i<n;i++)m.rows[i]=Vektor<typ>::standardVector(n,i);
      return m;
    }
  void append(Matrix const &m)
    {
      for(int i=0;i<m.height;i++)
	{
	  rows.push_back(m[i]);
	}
      height+=m.height;
    }
  IntegerVectorList getRows()const
    {
      IntegerVectorList ret;
      for(int i=0;i<height;i++)ret.push_back(rows[i]);
      return ret;
    }
  IntegerVector vectormultiply(IntegerVector const &v)const
    {
      assert(v.size()==width);
      IntegerVector ret(height);
      for(int i=0;i<height;i++)
	ret[i]=dot(rows[i],v);
      return ret;
    }
  /**
   * Decides if v is in the kernel of the matrix.
   */
  bool inKernel(IntegerVector const &v)const
    {
      assert(v.size()==width);
      for(int i=0;i<height;i++)
    	  if(dotLong(rows[i],v)!=0)return false;
      return true;
    }
  inline friend Matrix operator*(typ s, const Matrix& q)
    {
      Matrix p=q;
      for(int i=0;i<q.height;i++)p[i]=s*q[i];
      return p;
    }
  friend Matrix tropicalProduct(const Matrix& a, const Matrix& b)
  {
    int neutral=-888888888;
    assert(a.width==b.height);
    Matrix ret(a.height,b.width);
    for(int i=0;i<a.height;i++)
      for(int j=0;j<b.width;j++)
        {
          int sum=neutral;
          for(int k=0;k<b.height;k++)
            {
              if(a[i][k]+b[k][j]>sum)sum=a[i][k]+b[k][j];
            }
          ret[i][j]=sum;
        }
    return ret;
  }
  friend Matrix operator*(const Matrix& a, const Matrix& b)
    {
      assert(a.width==b.height);
      Matrix ret(b.width,a.height);
      for(int i=0;i<b.width;i++)
	ret[i]=a.vectormultiply(b.column(i));
      return ret.transposed();
    }
  /*  template<class T>
    Matrix<T>(const Matrix<T>& c):v(c.size()){
    for(int i=0;i<size();i++)v[i]=typ(c[i]);}
  */

  /**
     Returns the specified submatrix. The endRow and endColumn are not included.
   */
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
	ret[i-startRow][j-startColumn]=rows[i][j];
    return ret;
  }
  const Vektor<typ>& operator[](int n)const{assert(n>=0 && n<getHeight());return (rows[n]);}
//  Vektor<typ>& operator[](int n){assert(n>=0 && n<getHeight());return (rows[n]);}
// Bugfix for gcc4.5 (passing assertion to the above operator:
  Vektor<typ>& operator[](int n){if(!(n>=0 && n<getHeight()))(*(const Matrix<typ>*)(this))[n];return (rows[n]);}
};

typedef Matrix<int> IntegerMatrix;
typedef Matrix<double> FloatMatrix;

IntegerMatrix rowsToIntegerMatrix(IntegerVectorList const &rows, int width=-1);//width specifies the matrix width. If no width is specied the width is found by looking at the length of the rows. The function "asserts" if the length of the rows does not match the matrix size or if the width was not specified and could not be read off from the rows.


FloatMatrix integerToFloatMatrix(IntegerMatrix const &m);
IntegerVector flattenMatrix(IntegerMatrix const &m);
int rank(IntegerMatrix const &m);

#endif
