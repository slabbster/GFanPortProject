#ifndef LLL_H_INCLUDED
#define LLL_H_INCLUDED

#include "vektor.h"
#include "matrix.h"

using namespace std;

IntegerMatrix	mlll(IntegerMatrix &b, IntegerMatrix *inverseTransposedM=0);//LLL reduce the rows of b. The non-zero rows of b are put at the bottom rows. The returned matrix M satisfies: newb=M*oldb. The inverse transposed of M can also be computed.
IntegerMatrix latticeKernelOfTransposed(IntegerMatrix const &B);//The rows of the returned matrix is a lattice basis for the integer kernal of the transposed of B.
//int rank(basis_i B);

#endif
