#ifndef POLYNOMIALRING_H_INCLUDED
#define POLYNOMIALRING_H_INCLUDED

#include <string>
#include <vector>
#include <set>

#include "field.h"
#include "vektor.h"

using namespace std;
class PolynomialRingImplementation
{
 public:
  int refCount;
  Field theField;
  int n;
  vector<string> variableNames;
  PolynomialRingImplementation(Field const &f, int numberOfVariables, vector<string> const &variableNames_):
    theField(f),
    n(numberOfVariables),
    refCount(0),
    variableNames(variableNames_)
    {
    }
};


class PolynomialRing
{
  PolynomialRingImplementation *implementingObject;
 public:
  inline int getNumberOfVariables()const{return implementingObject->n;}
  inline Field const&getField()const{return implementingObject->theField;}
/**
 * Returns a new polynomial ring, the same as *this, except that some variables have been added.
 * The variables to append are given in a comma separated string.
 */
  PolynomialRing withVariablesAppended(string variableNames)const;
  int variableIndex(string const &name)const;//returns -1 if no match
  string const &getVariableName(int i)const;
  vector<string> getVariableNames()const;
/**
 * Returns a string with comma separated of variable names used for printing.
 */
  string toStringVariableNames()const;

  class Polynomial one()const;
  class Polynomial zero()const;
  class Polynomial monomialFromExponentVector(IntegerVector const &v)const;
  class Polynomial ithVariable(int i)const;
  class Polynomial polynomialFromField(FieldElement const &c)const;

  //construtors
  PolynomialRing(Field const &f, int numberOfVariables);
  PolynomialRing(Field const &f, vector<string> const &variables);
  ~PolynomialRing();
  PolynomialRing(PolynomialRing const &a);//copy constructor
  PolynomialRing& operator=(const PolynomialRing& a);//assignment
};


vector<string> matrixVariableNames(string base, int height, int width);
vector<string> vectorVariableNames(string base, int n);
vector<string> subsetVariableNames(string base, int n, int choose, bool M2Convention);
int subsetToVariableIndex(set<int> const &s, int n, int choose, bool M2Convention);

#endif
