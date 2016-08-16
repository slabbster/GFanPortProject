#ifndef PARSER_H
#define PARSER_H

#include <string>

#include "vektor.h"
#include "monomial.h"
#include "polynomial.h"
#include "field.h"
#include "polynomialring.h"

using namespace std;

class Parser
{
 protected:
  void parserError(const char *expected, char c);
 public:
  virtual int parseInt()=0;
  virtual Monomial parseMonomial(PolynomialRing const &r)=0;
  virtual IntegerVector parseIntegerVector()=0;
};


class CharacterBasedParser : public Parser
{
 private:
  PolynomialRing azAZ(Field const &f, int n=52);
  bool isVariable(int c);
  bool isDigit(int c);
  int variableIndex(int c);
 protected:
  virtual int getChar()=0;
  virtual void ungetChar(int c)=0;
 public:
  int nextNonBlank();
  int nextNonBlankDoNotGet();
  bool isLeftBracket(int c);
  bool isRightBracket(int c);
  int parseChar();
  int parseInt();
  double parseFloat();
  ComplexNumber parseComplexNumber();
  FieldElement parseFieldElement(Field const &f);
  FieldElement parseFieldElementFromInteger(Field const &f);
  Monomial parseMonomial(PolynomialRing const &r);
  IntegerVector parseIntegerVector();
  FloatVector parseFloatVector();
  ComplexVector parseComplexVector();
  IntegerVectorList parseIntegerVectorList();
  IntegerVectorList parseIntegerVectorList4ti2();
  Term parseTerm(PolynomialRing const &r);
  Field parseField();
  string parseVariableName();
  vector<string> parseVariableList();
  PolynomialRing parsePolynomialRing();
  Polynomial parsePolynomial(PolynomialRing const &r);
  Polynomial parsePolynomialWithRing();
  PolynomialSet parsePolynomialSet(PolynomialRing const &r);
  PolynomialSet parsePolynomialSetWithRing();
  PolynomialSetList parsePolynomialSetList(PolynomialRing const &r);
  PolynomialSetList parsePolynomialSetListWithRing();
};


class FileParser : public CharacterBasedParser
{
  FILE *f;
 protected:
  virtual int getChar();
  virtual void ungetChar(int c);
 public:
  FileParser(FILE *f);
};


class StringParser : public CharacterBasedParser
{
  const char *s;
  int index;
  bool hasUngotten;
  char ungotten;
 protected:
  virtual int getChar();
  virtual void ungetChar(int c);
 public:
  StringParser(const char *s);
};

#endif
