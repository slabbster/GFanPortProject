#ifndef MINKOWSKI_H_INCLUDED
#define MINKOWSKI_H_INCLUDED

#include "vektor.h"

typedef list<IntegerVectorList> IntegerVectorListList;

IntegerVectorList minkowski(IntegerVectorListList const &polytopes, IntegerVectorListList *sums);

#endif
