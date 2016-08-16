#ifndef REGULARSUBDIVISION_H_INCLUDED
#define REGULARSUBDIVISION_H_INCLUDED

#include <stdio.h>
#include "matrix.h"
#include <set>

set<set<int> > regularSubdivision(IntegerMatrix const &m, IntegerVector const &w);
void printSetSetInt(FILE *f, set<set<int> > const &faces);
//PolyhedralCone secondaryCone(IntegerMatrix const &m, set<set<int> > const &subdivision);

#endif
