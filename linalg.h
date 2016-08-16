#ifndef LINALG_H_INCLUDED
#define LINALG_H_INCLUDED

#include <vector>
#include "field.h"
#include "field_rationals.h"
#include "printer.h"
#include "matrix.h"

/*
  Unfortunately it seems that these classes cannot share code with the
  Vektor and Matrix templates since the classes in this file carry
  typeinfo around. In particular it is not possible to use the Matrix
  template on the FieldElement type since a FieldElement has no
  default initialization from a double -- its initialization really
  depends on which Field we are using.
 */

using namespace std;

/**
   @brief A vector containing FieldElement s.

   When a FieldVector is constructed the Field of which the entries of
   the vector will be members must be specified either explicitly or
   implicitly.
*/
class FieldVector
{
  /**
     The field of which the entries of the vector are members.
   */
  Field theField;
  /**
     The length of the vector.
  */
  int n;
  /**
     The entries of the vector.
   */
  vector<FieldElement> v;
 public:
   /**
     Constructs a FieldVector with entries initialized to 0.

     @param _theField the Field of which the entries of the vector will be members.
     @param _n the number of entries in the vector.
   */
  FieldVector(Field const &_theField, int _n);
  /**
   * Produces a vector with all entries equal to 1.
   */
  static FieldVector allOnes(Field const &_theField, int _n);
  /**
      Returns the Field.
  */
  inline Field const&getField()const{return theField;}
  /**
     Provides access to i th entry of the vector.
   */
  FieldElement& operator[](int i);
  /**
     Provides access to i th entry of the vector.
   */
  FieldElement const& operator[](int i)const;
  /**
     Returns true if all entries of the vector equal the 0-element of the Field, false otherwise.
   */
  bool isZero()const;
  /**
     Computes the first integer vector on the half ray spanned by the
     vector. Assumes that the Field is Q, asserts if this is not the
     case. Asserts if the computed vector has entries that do not fit
     in the word size.
   */
  IntegerVector primitive()const;
  /**
     Returns the number of entries in the vector.
   */
  int size()const;
  /**
     Returns the subvector (v_begin,\dots,v_end-1).
   */
  FieldVector subvector(int begin, int end)const;
  /**
     Returns the subvector .
   */
  FieldVector subvector(list<int> const &l)const;
  /**
     Returns the concatenation of the two vectors.
   */
  friend FieldVector concatenation(FieldVector const &a, FieldVector const &b);
  /**
     Computes the dot product of the vectors a and b.
   */
  friend FieldElement dot(FieldVector const &a, FieldVector const &b);
  /**
     Computes the vector q scaled by a factor s.
   */
  friend FieldVector operator*(FieldElement s, const FieldVector& q);
  /**
     Computes the coordinate-wise subtraction a-b.
   */
  friend FieldVector operator-(FieldVector const &a, FieldVector const &b);
  /**
     Computes the coordinate-wise division of a by b.
   */
  friend FieldVector operator/(FieldVector const &a, FieldVector const &b);
  /**
     Adds q to the current vector coordinatewise. Asserts if the
     current vector and q do not have the same number of elements.
   */
  FieldVector& operator+=(const FieldVector& q);
  /**
     Adds s*q to the current vector coordinatewise. Asserts if the
     current vector and q do not have the same number of elements.
   */
  FieldVector& madd(const FieldElement &s, const FieldVector& q);
  /**
     Prints the vector using the Printer P.
   */
  void print(Printer &P)const;

  friend IntegerVector fieldVectorToIntegerVector(FieldVector const &v);//MUST CONTAIN INTEGER ENTRIES ONLY // should this be a member?
  bool operator==(FieldVector const &b)const;
};


/**
   @brief A matrix containing FieldElement s.

   When a FieldMatrix is constructed the Field of which the entries of
   the matrix will be members must be specified either explicitly or
   implicitly. The FieldMatrix class contains routines for performing
   Gauss reductions.
*/
class FieldMatrix
{
  Field theField;
  int width,height;
  vector<FieldVector> rows;
 public:
   /**
     Constructs a FieldMatrix with entries initialized to 0.

     @param _theField the Field of which the entries of the matrix will be members.
     @param _height the number of rows of the matrix.
     @param _width the number of columns of the matrix.
   */
  FieldMatrix(Field const &_theField, int _height, int _width);
  /**
     Returns the nxn identity matrix.
  */
  static FieldMatrix identity(Field const &_theField, int _n);
  /**
     Returns the Field.
  */
  inline Field const&getField()const{return theField;}
  /**
     Provides access to i th row of the matrix which is a FieldVector.
  */
  FieldVector& operator[](int i);
  /**
     Provides access to i th row of the matrix which is a FieldVector.
  */
  FieldVector const& operator[](int i)const;
  /**
     Returns the number of rows of the matrix.
  */
  int getHeight()const;
  /**
     Returns the number of columns of the matrix.
  */
  int getWidth()const;
  /**
     Swaps the i th and the j th row.
   */
  void swapRows(int i, int j);
  /**
     Returns the submatrix consisting of rows index by rowIndices.
   */
  FieldMatrix submatrixRows(list<int> const &rowIndices)const;
  /**
     Adds a times the i th row to the j th row.
  */
  void madd(int i, FieldElement a, int j);
  /**
     Returns the index of a row whose index is at least currentRow and
     which has a non-zero element on the column th entry. If no such
     row exists then -1 is returned. This routine is used in the Gauss
     reduction. To make the reduction more efficient the routine
     chooses its row with as few non-zero entries as possible.
   */
  int findRowIndex(int column, int currentRow)const;
  /**
     Prints the matrix using the Printer P.
   */
  void printMatrix(Printer &P)const;
  /**
     This method is used for iterating through the pivots in a matrix
     in row echelon form. To find the first pivot put i=-1 and
     j=-1 and call this routine. The (i,j) th entry of the matrix is a
     pivot. Call the routine again to get the next pivot. When no more
     pivots are found the routine returns false.
  */
  bool nextPivot(int &i, int &j)const;
  /**
     Performs a Gauss reduction and returns the number of row swaps
     done.  The result is a matrix in row echelon form. The pivots may
     not be all 1.  In terms of Groebner bases, what is computed is a
     minimal (not necessarily reduced) Groebner basis of the linear
     ideal generated by the rows.  The number of swaps is need if one
     wants to compute the determinant afterwards. In this case it is
     also a good idea to set the flag returnIfZeroDeterminant which
     make the routine terminate before completion if it discovers that
     the determinant is zero.
  */
  int reduce(bool returnIfZeroDeterminant=false, bool hermite=false); //bool reducedRowEcholonForm,
  /**
     Calls reduce() and returns the determinant of the matrix. If
     the matrix is not square the routine asserts.
   */
  FieldElement reduceAndComputeDeterminant();
  /**
     Calls reduce() and returns the rank of the matrix.
   */
  int reduceAndComputeRank();
  /**
     Scales every row of the matrix so that the first non-zero entry of
     each row is 1.
   */
  void scaleFirstNonZeroEntryToOne();
  /**
     Assumes that the matrix has a kernel of dimension 1.
     Calls reduce() and returns a non-zero vector in the kernel.
     If the matrix is an (n-1)x(n) matrix then the returned vector has
     the property that if appended as a row to the original matrix
     then the determinant of that matrix would be positive. Of course
     this property, as described here, only makes sense for ordered fields.
   */
  FieldVector reduceAndComputeVectorInKernel();
  /**
     Calls reduce() and constructs matrix whose rows forms a basis of
     the kernel of the linear map defined by the original matrix. The
     return value is the new matrix.
   */
  FieldMatrix reduceAndComputeKernel();
  /**
     Returns the transposed of the the matrix.
   */
  FieldMatrix transposed()const;
  /**
     Returns the matrix with rows in reversed order.
   */
  FieldMatrix flipped()const;
  /**
     Takes two matrices with the same number of columns and construct
     a new matrix which has the rows of the matrix top on the top and
     the rows of the matrix bottom at the bottom. The return value is
     the constructed matrix.
   */
  friend FieldMatrix combineOnTop(FieldMatrix const &top, FieldMatrix const &bottom);
  /**
     Takes two matrices with the same number of rows and construct
     a new matrix which has the columns of the matrix left on the left and
     the columns of the matrix right on the right. The return value is
     the constructed matrix.
   */
  friend FieldMatrix combineLeftRight(FieldMatrix const &left, FieldMatrix const &right);
  /**
     Computes a reduced row echelon form from a row echelon form. In
     Groebner basis terms this is the same as tranforming a minimal
     Groebner basis to a reduced one except that we do not force
     pivots to be 1 unless the scalePivotsToOne parameter is set.
   */
  int REformToRREform(bool scalePivotsToOne=false);
  /**
     This function may be called if the FieldMatrix is in Row Echelon
     Form. The input is a FieldVector which is rewritten modulo the
     rows of the matrix. The result is unique and is the same as the
     normal form of the input vector modulo the Groebner basis of the
     linear ideal generated by the rows of the matrix.
  */
  FieldVector canonicalize(FieldVector v)const;
  /**
     This routine removes the zero rows of the matrix.
   */
  void removeZeroRows();
  /**
     This routine produces a matrix in row Echelon form solving a
     system of the kind Ax=b, where A equals *this. The way this is
     done is by considering [A^t -I] and computing its RE form.  To
     solve the system afterwards simply call canonicalize of the new
     matrix on (b_1,\dots,b_m,0,\dots,0) to get
     (0,\dots,0,x_1,\dots,x_n).  If no solution exists then the first
     m entries of the returned vector are not all 0.
   */
  FieldMatrix solver()const;
  /**
     Computes the vector m scaled by a factor s.
   */
  friend FieldMatrix operator*(FieldElement s, const FieldMatrix& m);
  /**
     Returns the matrix product a*b.
   */
  friend FieldMatrix operator*(FieldMatrix const &a, const FieldMatrix &b);
  /**
     Returns the vector a*b, where a is considered as a row vector.
   */
  friend FieldVector operator*(FieldVector const &a, const FieldMatrix &b);
  /**
     Compares two matrices.
  */
  bool operator==(FieldMatrix const &b)const;
  /**
     Returns the rank of the matrix. This function is const and slow -
     will copy the matrix and perform Gauss reduction.
   */
  int rank()const;
  /**
     Returns the indices of the columns containing a pivot.
     The returned list is sorted.
     The matrix must be in row echelon form.
   */
  list<int> pivotColumns()const;
  /**
     Returns the indices of the columns not containing a pivot.
     The returned list is sorted.
     The matrix must be in row echelon form.
   */
  list<int> nonPivotColumns()const;
  /**
   * Returns the specified submatrix with size (endRow-startRow)x(endColumn-startColumn).
   */
  FieldMatrix submatrix(int startRow, int startColumn, int endRow, int endColumn)const;
  /**
   * Thinking of an upper triangular matrix as a Groebner basis this routine returns the normal form
   * of the linear polynomial represented by v. If coeff is non-zero, then that vector is assigned
   * the coefficients used in the reduction. That is, if the normal form is zero, then multiplying coeff
   * and *this gives v. The matrix must be in row-echelon form. coeff does not have to have the right length before the call.
   */
  FieldVector normalForm(FieldVector const &v, FieldVector *coeff)const;
};


FieldMatrix integerMatrixToFieldMatrix(IntegerMatrix const &m, Field f);
IntegerMatrix fieldMatrixToIntegerMatrix(FieldMatrix const &m);//MUST CONTAIN INTEGER ENTRIES ONLY
IntegerMatrix fieldMatrixToIntegerMatrixPrimitive(FieldMatrix const &m);//Scales to primitive vector
FieldVector integerVectorToFieldVector(IntegerVector const &v, Field f);

/**
   Computes a non-zero vector in the kernel of the given matrix. The
   dimension of the kernel is assumed to be 1. The returned vector is primitive.
*/
IntegerVector vectorInKernel(IntegerMatrix const &m);

/**
 * Using exact arithmetics this routine chooses a subset of the vectors in m, which generate the span of m.
 */
IntegerVectorList subsetBasis(IntegerVectorList const &m);

#endif
