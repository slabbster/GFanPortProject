#include "matrix.h"

IntegerMatrix rowsToIntegerMatrix(IntegerVectorList const &rows, int width)
{
  int height=rows.size();
  int a=-1;
  for(IntegerVectorList::const_iterator i=rows.begin();i!=rows.end();i++)
    {
      if(a==-1)
	a=i->size();
      else
	{
	  assert(a==i->size());
	}
    }
  if(a==-1)
    {
      assert(width!=-1);
    }
  else
    {
      if(width!=-1)
	{
	  assert(width==a);
	}
    }
  if(width==-1)width=a;
  IntegerMatrix ret(height,width);

  int j=0;
  for(IntegerVectorList::const_iterator i=rows.begin();i!=rows.end();i++)
    ret[j++]=*i;

  return ret;
}


FloatMatrix integerToFloatMatrix(IntegerMatrix const &m)
{
  FloatMatrix ret(m.getHeight(),m.getWidth());

  for(int i=0;i<m.getHeight();i++)
    for(int j=0;j<m.getWidth();j++)
      ret[i][j]=m[i][j];

  return ret;
}


IntegerVector flattenMatrix(IntegerMatrix const &m)
{
  IntegerVector ret(m.getHeight()*m.getWidth());

  for(int i=0;i<m.getHeight();i++)
    for(int j=0;j<m.getWidth();j++)
      ret[i*m.getWidth()+j]=m[i][j];

  return ret;
}


#include "linalg.h"

int rank(IntegerMatrix const &m)
{
  return integerMatrixToFieldMatrix(m,Q).rank();
}
