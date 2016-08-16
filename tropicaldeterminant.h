#ifndef TROPICALDETERMINANT_H_INCLUDED
#define TROPICALDETERMINANT_H_INCLUDED

#include "matrix.h"

#define td_minusInfinity 0x80000000
int tropicalDeterminant(IntegerMatrix const &m);
void tropicalDeterminantTest();

#endif
