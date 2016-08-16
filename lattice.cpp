#include "lattice.h"

#include "linalg.h"
#include "field_rationals.h"


bool isPartOfAZBasis(IntegerVectorList const &l)
{
  if(l.size()==0)return true;
  IntegerMatrix A=rowsToIntegerMatrix(l).transposed();
  FieldMatrix Af=integerMatrixToFieldMatrix(A,Q);
  Af.reduce(false,true);

  int i=-1;
  int j=-1;
  int nPivots=0;
  while(Af.nextPivot(i,j))
    {
      if(!((Af[i][j]*Af[i][j])-Q.zHomomorphism(1)).isZero())return false;
      nPivots++;
    }
  return l.size()==nPivots;
}
