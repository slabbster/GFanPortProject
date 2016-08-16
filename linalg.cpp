#include "linalg.h"
#include "field_rationals.h"
#include "linalgfloat.h"
#include "continuedfractions.h"
#include "log.h"

//--------------------------------------------------
// FieldVector
//--------------------------------------------------

FieldVector::FieldVector(Field const &_theField, int _n):
  n(_n),
  theField(_theField),
  v(_n,_theField.zHomomorphism(0))
{
}


FieldVector FieldVector::allOnes(Field const &_theField, int _n)
{
  FieldVector ret(_theField,_n);
  for(int i=0;i<_n;i++)ret[i]=_theField.zHomomorphism(1);
  return ret;
}


FieldElement& FieldVector::operator[](int i)
{
  assert(i>=0 && i<n);
  return (v[i]);
}


FieldElement const& FieldVector::operator[](int i)const
{
  assert(i>=0 && i<n);
  return (v[i]);
}


bool FieldVector::isZero()const
{
  for(int i=0;i<n;i++)if(!v[i].isZero())return false;
  return true;
}


IntegerVector FieldVector::primitive()const
{
//  debug<<*this<<"\n";
  return primitiveVector(v);
}


int FieldVector::size()const
{
  return n;
}


FieldVector FieldVector::subvector(int begin, int end)const
{
  assert(begin>=0);
  assert(end<=size());
  assert(end>=begin);

  FieldVector ret(theField,end-begin);
  for(int i=0;i<end-begin;i++)
    ret[i]=v[begin+i];

  return ret;
}


FieldVector FieldVector::subvector(list<int> const &l)const
{
  FieldVector ret(theField,l.size());
  int I=0;
  for(list<int>::const_iterator i=l.begin();i!=l.end();i++,I++)
    ret[I]=v[*i];

  return ret;
}


FieldVector concatenation(FieldVector const &a, FieldVector const &b)
{
  FieldVector ret(a.getField(),a.size()+b.size());
  int I=0;
  for(int i=0;i<a.size();i++,I++)ret[I]=a[i];
  for(int i=0;i<b.size();i++,I++)ret[I]=b[i];

  return ret;
}


FieldElement dot(FieldVector const &a, FieldVector const &b)
{
  FieldElement ret(a.theField);
  assert(a.size()==b.size());
  for(int i=0;i<a.size();i++)
    ret.madd(a[i],b[i]);
//  ret=ret+a[i]*b[i];

  return ret;
}

FieldVector operator*(FieldElement s, const FieldVector& q)
{
  FieldVector p=q;
  for(int i=0;i<q.size();i++)p[i]*=s;
  return p;
}


FieldVector operator-(FieldVector const &a, FieldVector const &b)
{
  FieldVector ret(a.getField(),a.size());
  for(int i=0;i<ret.size();i++)ret[i]=a[i]-b[i];
  return ret;
}


FieldVector operator/(FieldVector const &a, FieldVector const &b)
{
  FieldVector ret(a.getField(),a.size());
  for(int i=0;i<ret.size();i++)ret[i]=a[i]*(b[i].inverse());
  return ret;
}


FieldVector& FieldVector::operator+=(const FieldVector& q)
{
  assert(size()==q.size());
  for(int i=0;i<size();i++)v[i]=v[i]+q.v[i];
  return *this;
}


FieldVector& FieldVector::madd(const FieldElement &s, const FieldVector& q)
{
  assert(size()==q.size());
  for(int i=0;i<size();i++)
    v[i].madd(s,q.v[i]);
  //v[i]=v[i]+s*q.v[i];
  return *this;
}


void FieldVector::print(Printer &P)const
{
  for(int j=0;j<v.size();j++)
    {
      P.printFieldElement(v[j]);
      P.printString(" ");
    }
  P.printNewLine();
}


bool FieldVector::operator==(FieldVector const &b)const
{
  if(size()!=b.size())return false;

  FieldVector temp=*this;
  temp+=theField.zHomomorphism(-1)*b;

  return temp.isZero();
}


//--------------------------------------------------
// FieldMatrix
//--------------------------------------------------


FieldMatrix::FieldMatrix(Field const &_theField, int _height, int _width):
  theField(_theField),
  height(_height),
  width(_width),
  rows(_height,FieldVector(_theField,width))
{
}


FieldMatrix FieldMatrix::identity(Field const &_theField, int _n)
{
  FieldMatrix ret(_theField,_n,_n);
  for(int i=0;i<_n;i++)
    ret[i][i]=_theField.zHomomorphism(1);
  return ret;
}

FieldVector& FieldMatrix::operator[](int i)
{
  assert(i>=0 && i<height);
  return (rows[i]);
}


FieldVector const& FieldMatrix::operator[](int i)const
{
  assert(i>=0 && i<height);
  return (rows[i]);
}


int FieldMatrix::getHeight()const
{
  return height;
}


int FieldMatrix::getWidth()const
{
  return width;
}


void FieldMatrix::swapRows(int i, int j)
{
  assert(i!=j);

  FieldVector temp=rows[i];
  rows[i]=rows[j];
  rows[j]=temp;
}


FieldMatrix FieldMatrix::submatrixRows(list<int> const &rowIndices)const
{
  FieldMatrix ret(getField(),rowIndices.size(),getWidth());
  int I=0;
  for(list<int>::const_iterator i=rowIndices.begin();i!=rowIndices.end();i++,I++)
    ret[I]=rows[*i];
  return ret;
}


void FieldMatrix::madd(int i, FieldElement a, int j)
{
  assert(i!=j);
  assert(i>=0 && i<height);
  assert(j>=0 && j<height);

  if(!a.isZero())
  for(int k=0;k<width;k++)
    if(!rows[i][k].isZero())
      rows[j][k].madd(rows[i][k],a);
//  rows[j][k]=rows[j][k]+rows[i][k]*a;
}


/*Basic strategy
int FieldMatrix::findRowIndex(int column, int currentRow)const
{
  for(int i=currentRow;i<height;i++)
    if(!rows[i][column].isZero())return i;
  return -1;
}
*/

 //Advanced strategy
int FieldMatrix::findRowIndex(int column, int currentRow)const
{
  int best=-1;
  int bestNumberOfNonZero=0;
  for(int i=currentRow;i<height;i++)
    if(!rows[i][column].isZero())
      {
	int nz=0;
	for(int k=column+1;k<width;k++)
	  if(!rows[i][k].isZero())nz++;
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


void FieldMatrix::printMatrix(Printer &P)const
{
  if(0)
    {
      for(int i=0;i<height;i++)
	{
	  for(int j=0;j<width;j++)
	    {
	      P.printFieldElement(rows[i][j]);
	      P.printString(" ");
	    }
	  P.printNewLine();
	}
      P.printNewLine();
    }
  else
    {
      vector<int> widths(getWidth());
      for(int j=0;j<widths.size();j++)
	{
	  int maxWidth=0;
	  for(int i=0;i<getHeight();i++)
	    {
	      int l=rows[i][j].toString().size();
	      if(l>maxWidth)maxWidth=l;
	    }
	  widths[j]=maxWidth;
	}
      for(int i=0;i<height;i++)
	{
	  for(int j=0;j<width;j++)
	    {
	      if(j)P.printString(" ");
	      std::string s=rows[i][j].toString();
	      for(int k=s.size();k<widths[j];k++)P.printString(" ");
	      P.printString(s);
	    }
	  P.printNewLine();
	}
      P.printNewLine();
    }
}


bool FieldMatrix::nextPivot(int &i, int &j)const//iterates through the pivots in a matrix in reduced row echelon form. To find the first pivot put i=-1 and j=-1 and call this routine. When no more pivots are found the routine returns false.
{
  i++;
  if(i>=height)return false;
  while(++j<width)
    {
      if(!rows[i][j].isZero()) return true;
    }
  return false;
}


int FieldMatrix::REformToRREform(bool scalePivotsToOne)
{
  int ret=0;
  int pivotI=-1;
  int pivotJ=-1;
  while(nextPivot(pivotI,pivotJ))
    {
      if(scalePivotsToOne)
	rows[pivotI]=rows[pivotI][pivotJ].inverse()*rows[pivotI];
      for(int i=0;i<pivotI;i++)
        if(!rows[i][pivotJ].isZero())
          madd(pivotI,-rows[i][pivotJ]*rows[pivotI][pivotJ].inverse(),i);
    }
  return ret;
}

int FieldMatrix::reduce(bool returnIfZeroDeterminant, bool hermite) //bool reducedRowEcholonForm,
/* Run a Gauss reduction. Returns the number of swaps. The number of
   swaps is need if one wants to compute the determinant
   afterwards. In this case it is also a good idea to set the flag
   which make the routine terminate when a it is discovered the the
   determinant is zero. */
{
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
	    if(hermite)
	      {
		if(!rows[j][i].isZero())
		  {
		    FieldElement s(theField);
		    FieldElement t(theField);

		    FieldElement g=gcd(rows[currentRow][i],rows[j][i],s,t);
		    FieldElement u=-rows[j][i]*(g.inverse());
		    FieldElement v=rows[currentRow][i]*(g.inverse());
		    /* We want the (s,t) vector to be as small as possible.
		     * We are allowed to adjust by multiples of (u,v).
		     * The following computes the correct multiplier (in most cases).
		     */
		    {
		      FieldElement multiplier=(s*u+t*v)*((u*u+v*v).inverse());
		      double d=mpq_get_d(*(multiplier.getGmpRationalTemporaryPointer()));
		      multiplier=multiplier.getField()->zHomomorphism(-(((int)(d+10000.5))-10000));
                      //s=s+multiplier*u;
                      //t=t+multiplier*v;
                      s.madd(multiplier,u);
                      t.madd(multiplier,v);
		    }
/*                    if(((s-u)*(s-u)+(t-v)*(t-v)-(s*s+t*t)).sign()<0)//a simple test to avoid coefficient growth
                      {
                        s=s-u;
                        t=t-v;
                      }
                    if(((s+u)*(s+u)+(t+v)*(t+v)-(s*s+t*t)).sign()<0)//a simple test to avoid coefficient growth
                      {
                        s=s+u;
                        t=t+v;
                      }*/
//		    debug<<"i"<<i<<"s"<<s<<" t"<<t<<" g"<<g<<" u"<<u<<" v"<<v<<"\n";
		    for(int k=0;k<width;k++)
		      {
			FieldElement A=rows[currentRow][k];
			FieldElement B=rows[j][k];

			rows[currentRow][k]=s*A+t*B;
			rows[j][k]=u*A+v*B;
		      }
		  }
	      }
	    else
	      {
                if(!rows[j][i].isZero())
                  madd(currentRow,-rows[j][i]*rows[currentRow][i].inverse(),j);
	      }
	  currentRow++;
	}
      else
	if(returnIfZeroDeterminant)return -1;
//      debug<<*this;
    }

  return retSwaps;
}


int FieldMatrix::reduceAndComputeRank()
{
  reduce();
  int ret=0;
  /*      for(int i=0;i<height;i++)
	      if(!rows[i].isZero())ret++;*/
  int pivotI=-1;
  int pivotJ=-1;
  while(nextPivot(pivotI,pivotJ))ret++;
  return ret;
}


FieldElement FieldMatrix::reduceAndComputeDeterminant()
{
  assert(height==width);
  int swaps=reduce();

  int r=reduceAndComputeRank();

  if(r!=height)return theField.zHomomorphism(0);

  FieldElement ret=theField.zHomomorphism(1-2*(swaps&1));

  int pivotI=-1;
  int pivotJ=-1;
  while(nextPivot(pivotI,pivotJ))ret=ret*rows[pivotI][pivotJ];
  return ret;
}


FieldVector FieldMatrix::reduceAndComputeVectorInKernel()
{
  //FieldMatrix temp=*this;debug<<"comment this line\n";
  // TODO: (optimization) if the field is ordered, then it is better to just keep track of signs rather than
  // multiplying by sign*diagonalProduct*lastEntry at the end.
  int nswaps=this->reduce();
  FieldElement sign=theField.zHomomorphism(1-2*(nswaps&1));
  int rank=reduceAndComputeRank();
  assert(rank+1==width);

  REformToRREform();

  FieldVector ret(theField,width);

  FieldElement diagonalProduct=theField.zHomomorphism(1);
  {
    int pivot2I=-1;
    int pivot2J=-1;
    while(nextPivot(pivot2I,pivot2J))
      {
        diagonalProduct*=rows[pivot2I][pivot2J];
      }
  }
  {
    int j=nonPivotColumns().front();
    int pivot2I=-1;
    int pivot2J=-1;
    ret[j]=theField.zHomomorphism(-1);
    // Pretend that we are appending ret to the matrix, and reducing this
    // new row by the previous ones. The last entry of the resulting matrix
    // is lastEntry.
    FieldElement lastEntry=ret[j];
    while(nextPivot(pivot2I,pivot2J))
      {
        ret[pivot2J]=rows[pivot2I][j]*rows[pivot2I][pivot2J].inverse();
        //lastEntry+=ret[pivot2J]*ret[pivot2J];
        lastEntry+=-ret[pivot2J]*ret[pivot2J];
      }
    ret=(sign*diagonalProduct*lastEntry)*ret;
  }

/*  {
    int n=temp.getWidth();
    FieldMatrix temp2(temp.getField(),n,n);
    for(int i=0;i<n-1;i++)for(int j=0;j<n;j++)temp2[i][j]=temp[i][j];
    for(int j=0;j<n;j++)temp2[n-1][j]=ret[j];
    assert(temp2.reduceAndComputeDeterminant().sign()>0);
  }
*/
  return ret;
}


FieldMatrix FieldMatrix::reduceAndComputeKernel()
{
  FieldMatrix ret(theField,width-reduceAndComputeRank(),width);

  REformToRREform();

  //  AsciiPrinter P(Stderr);
  //  printMatrix(P);

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
      //      log0 fprintf(stderr,"%i\n",j);
      int pivot2I=-1;
      int pivot2J=-1;
      while(nextPivot(pivot2I,pivot2J))
	{
	  ret[k][pivot2J]=rows[pivot2I][j]*rows[pivot2I][pivot2J].inverse();
	}
      ret[k][j]=theField.zHomomorphism(-1);
      k++;
    }
  return ret;
}


void FieldMatrix::scaleFirstNonZeroEntryToOne()
{
  for(int i=0;i<height;i++)
    for(int j=0;j<width;j++)
      if(!rows[i][j].isZero())
	{
	  if(!rows[i][j].isOne())
	  {
	    FieldElement inv=rows[i][j].inverse();
	    for(int J=j;J<width;J++)
	      if(!rows[i][J].isZero())rows[i][J]*=inv;
	  }
//	  rows[i]=rows[i][j].inverse()*rows[i];
	  break;
	}
}


FieldMatrix FieldMatrix::transposed()const
{
  FieldMatrix ret(theField,width,height);
  for(int i=0;i<width;i++)
    for(int j=0;j<height;j++)
      ret[i][j]=rows[j][i];
  return ret;
}


FieldMatrix FieldMatrix::flipped()const
{
  FieldMatrix ret(theField,height,width);
  for(int i=0;i<height;i++)
    ret[i]=rows[height-1-i];
  return ret;
}


FieldVector FieldMatrix::canonicalize(FieldVector v)const
{
  assert(v.size()==getWidth());

  int pivotI=-1;
  int pivotJ=-1;

  while(nextPivot(pivotI,pivotJ))
    if(!v[pivotJ].isZero())
    {
      FieldElement s=-v[pivotJ]*rows[pivotI][pivotJ].inverse();

      for(int k=0;k<width;k++)
	if(!rows[pivotI][k].isZero())
          v[k].madd(rows[pivotI][k],s);
//      v[k]+=rows[pivotI][k]*s;
//      v[k]=v[k]+rows[pivotI][k]*s;
    }
  return v;
}

void FieldMatrix::removeZeroRows()
{
  int nonZeros=0;
  for(int i=0;i<height;i++)if(!(*this)[i].isZero())nonZeros++;
  if(nonZeros==height)return;

  FieldMatrix b(theField,nonZeros,width);

  int j=0;
  for(int i=0;i<height;i++)
    {
      if(!(*this)[i].isZero())
	{
	  b[j]=(*this)[i];
	  j++;
	}
    }
  *this=b;
}


FieldMatrix FieldMatrix::solver()const
{
  FieldMatrix ret=combineOnTop(*this,theField.zHomomorphism(-1)*identity(theField,getWidth())).transposed();

  ret.reduce();
  ret.REformToRREform();

  return ret;
}


FieldMatrix operator*(FieldElement s, const FieldMatrix& m)
{
  FieldMatrix ret=m;

  for(int i=0;i<ret.getHeight();i++)
    ret[i]=s*ret[i];

  return ret;
}

FieldMatrix operator*(FieldMatrix const &a, const FieldMatrix &b)
{
  assert(a.getWidth()==b.getHeight());

  FieldMatrix temp=b.transposed();
  FieldMatrix ret(Q,a.getHeight(),b.getWidth());

  for(int i=0;i<ret.getHeight();i++)
    for(int j=0;j<ret.getWidth();j++)
      ret[i][j]=dot(a[i],temp[j]);

  return ret;
}


FieldVector operator*(FieldVector const &a, const FieldMatrix &b)
{
  assert(a.size()==b.getHeight());
  FieldVector ret(b.getField(),b.getWidth());
  for(int i=0;i<a.size();i++)
    ret.madd(a[i],b[i]);
  return ret;
}


bool FieldMatrix::operator==(FieldMatrix const &b)const
{
  if(getHeight()!=b.getHeight())return false;
  if(getWidth()!=b.getWidth())return false;

  for(int i=0;i<getHeight();i++)
    if(!(rows[i]==b[i]))return false;

  return true;
}


int FieldMatrix::rank()const
{
  FieldMatrix temp=*this;
  return temp.reduceAndComputeRank();
}


list<int> FieldMatrix::pivotColumns()const
{
  list<int> ret;
  int pivotI=-1;
  int pivotJ=-1;
  while(nextPivot(pivotI,pivotJ))
    {
      ret.push_back(pivotJ);
    }
  return ret;
}


list<int> FieldMatrix::nonPivotColumns()const
{
  list<int> ret;
  int pivotI=-1;
  int pivotJ=-1;
  int firstPossiblePivot=0;
  while(nextPivot(pivotI,pivotJ))
    {
      for(int j=firstPossiblePivot;j<pivotJ;j++)
	ret.push_back(j);
      firstPossiblePivot=pivotJ+1;
    }
  for(int j=firstPossiblePivot;j<getWidth();j++)
    ret.push_back(j);

  return ret;
}


FieldMatrix FieldMatrix::submatrix(int startRow, int startColumn, int endRow, int endColumn)const
{
	assert(startRow<=endRow);
	assert(startColumn<=endColumn);
	assert(startRow>=0);
	assert(startColumn>=0);
	assert(endRow<=getHeight());
	assert(endColumn<=getWidth());
	FieldMatrix ret(theField,endRow-startRow,endColumn-startColumn);
	for(int i=startRow;i<endRow;i++)
		for(int j=startColumn;j<endColumn;j++)
			ret[i-startRow][j-startColumn]=rows[i][j];

	return ret;
}


FieldVector FieldMatrix::normalForm(FieldVector const &v, FieldVector *coeff)const
{
  if(coeff)*coeff=FieldVector(getField(),getHeight());
  FieldVector v2=v;
  int pivotI=-1;
  int pivotJ=-1;
  while(nextPivot(pivotI,pivotJ))
if(!v2[pivotJ].isZero())
  {
      FieldElement c=v2[pivotJ]*rows[pivotI][pivotJ].inverse();
      v2.madd(-c,(*this)[pivotI]);
      if(coeff)(*coeff)[pivotI]=c;
    }
  return v2;
}

//---------------------------------------


FieldMatrix integerMatrixToFieldMatrix(IntegerMatrix const &m, Field f)
{
  int height=m.getHeight();
  int width=m.getWidth();
  FieldMatrix ret(f,height,width);
  for(int i=0;i<height;i++)
    for(int j=0;j<width;j++)
      ret[i][j]=f.zHomomorphism(m[i][j]);
  return ret;
}


IntegerMatrix fieldMatrixToIntegerMatrix(FieldMatrix const &m)
{
  IntegerMatrix ret(m.getHeight(),m.getWidth());

  for(int i=0;i<m.getHeight();i++)
    ret[i]=fieldVectorToIntegerVector(m[i]);

  return ret;
}


IntegerMatrix fieldMatrixToIntegerMatrixPrimitive(FieldMatrix const &m)
{
  IntegerMatrix ret(m.getHeight(),m.getWidth());

  for(int i=0;i<m.getHeight();i++)
    ret[i]=m[i].primitive();

  return ret;
}


FieldVector integerVectorToFieldVector(IntegerVector const &v, Field f)
{
  FieldVector ret(f,v.size());
  for(int i=0;i<v.size();i++)
    ret[i]=f.zHomomorphism(v[i]);
  return ret;
}


IntegerVector fieldVectorToIntegerVector(FieldVector const &v)
{
  return toIntegerVector(v.v);
}


IntegerVector vectorInKernel(IntegerMatrix const &m) //computes a non-zero vector in the kernel of the given matrix. The dimension of the kernel is assumed to be 1.
{
	{
		// Try to use floating point operations and verify result.
		// Fall back on rational arithmetic in case of failure.
		linalgfloat::Matrix M(m.getHeight(),m.getWidth());
		for(int i=0;i<m.getHeight();i++)
			for(int j=0;j<m.getWidth();j++)
				M[i][j]=m[i][j];
		if(M.reduceAndComputeRank()!=M.getWidth()-1)goto fallback;
		linalgfloat::Matrix kernel=M.reduceAndComputeKernel();
		assert(kernel.getHeight()==1);
		vector<double> temp(kernel.getWidth());
		for(int i=0;i<kernel.getWidth();i++)temp[i]=kernel[0][i];
		vector<int> tempInt(kernel.getWidth());
		int denominator;
		doubleVectorToFractions(temp, tempInt, denominator);
		IntegerVector ret(kernel.getWidth());
		for(int i=0;i<kernel.getWidth();i++)ret[i]=tempInt[i];
		if(ret.max()>1000000)goto fallback;
		if(ret.min()<-1000000)goto fallback;
		if(ret.isZero())goto fallback;
		if(!m.inKernel(ret))goto fallback;
		return ret;
	}
	fallback:
	debug << "LINALG FALLBACK\n";


  FieldMatrix m2=integerMatrixToFieldMatrix(m,Q);
  FieldMatrix ker=m2.reduceAndComputeKernel();
  assert(ker.getHeight()==1);
  return ker[0].primitive();
}



FieldMatrix combineOnTop(FieldMatrix const &top, FieldMatrix const &bottom)
{
  assert(bottom.getWidth()==top.getWidth());
  FieldMatrix ret(top.theField,top.getHeight()+bottom.getHeight(),top.getWidth());
  for(int i=0;i<top.getHeight();i++)ret.rows[i]=top.rows[i];
  for(int i=0;i<bottom.getHeight();i++)ret.rows[i+top.getHeight()]=bottom.rows[i];

  return ret;
}


FieldMatrix combineLeftRight(FieldMatrix const &left, FieldMatrix const &right)
{
  assert(left.getHeight()==right.getHeight());
  FieldMatrix ret(left.theField,left.getHeight(),left.getWidth()+right.getWidth());
  for(int i=0;i<left.getHeight();i++)
    {
      for(int j=0;j<left.getWidth();j++)ret.rows[i][j]=left.rows[i][j];
      for(int j=0;j<right.getWidth();j++)ret.rows[i][j+left.getWidth()]=right.rows[i][j];
    }
  return ret;
}


IntegerVectorList subsetBasis(IntegerVectorList const &m)
{
  IntegerVectorList ret;
  int n=m.size();
  if(n==0)return ret;
  int d=m.begin()->size();
  FieldMatrix M(Q,d,n);
  int j=0;
  for(IntegerVectorList::const_iterator J=m.begin();J!=m.end();J++,j++)
    for(int i=0;i<d;i++)M[i][j]=Q.zHomomorphism((*J)[i]);
  M.reduce(false,false);

  int pivotI=-1;
  int pivotJ=-1;
  IntegerVectorList::const_iterator J=m.begin();
  j=0;
  while(M.nextPivot(pivotI,pivotJ))
    {
      while(pivotJ>j){j++;J++;}
      ret.push_back(*J);
    }

  return ret;
}
