#ifndef CONTINUEDFRACTIONS_H_INCLUDED
#define CONTINUEDFRACTIONS_H_INCLUDED

#include <vector>

using namespace std;

void doubleToFraction(double f, int &numerator, int &denominator, int maksIter=15);

void doubleVectorToFractions(vector<double> v, vector<int> &numerators, int &denominator);

#endif
