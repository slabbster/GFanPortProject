#ifndef TROPICAL_H_INCLUDED
#define TROPICAL_H_INCLUDED

#include "polynomial.h"

//bool isFullColored(IntegerVectorList const &inequalityColors, IntegerVector const &v)
PolynomialSetList fullColoredIdeals(PolynomialSet const &g, bool skipColorTest=false);
bool containsMonomial(PolynomialSet const &ideal);//ideal must be homogeneous
Term computeTermInIdeal(PolynomialSet const &ideal, int m=0);//ideal must be homogeneous
PolynomialSet saturatedIdeal(PolynomialSet const &ideal); //ideal must be homogeneous

#endif
