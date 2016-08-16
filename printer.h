#ifndef Printer_INCLUDED
#define Printer_INCLUDED

class Printer;

#include <stdio.h>
#include <string>
#include "vektor.h"
#include "term.h"
#include "termorder.h"
#include "polynomial.h"

// We need to help the Mac compiler:
#define Stdin ((FILE*)stdin)
#define Stdout ((FILE*)stdout)
#define Stderr ((FILE*)stderr)

class Printer
{
  static bool assertOnPrinting;
 protected:
  static void printCheck();
  FILE *f;
  virtual string variableIndexToString(PolynomialRing const &r, int i);
 public:
  static void setAssertOnPrinting(bool b);
  Printer(FILE *f){this->f=f;};
  virtual void printVariable(PolynomialRing const &r, int i);
  virtual void printInteger(int i, int minimalFieldWidth=0)=0;
  virtual void printFloat(double i, int minimalFieldWidth=0)=0;
  virtual void printComplexNumber(ComplexNumber const &i, int minimalFieldWidth=0)=0;
//  virtual void printComplex(class complex const &c, int minimalFieldWidth=0)=0;
  virtual void printMonomial(const Monomial &m, bool alwaysWriteSign=false, bool writeIfOne=true)=0;
  virtual void printFieldElement(const FieldElement &e, bool writeIfOne=true, bool alwaysWriteSign=false)=0;
  virtual void printTerm(const Term &t)=0;
  virtual void printPolynomial(const Polynomial &p)=0;
  virtual void printPolynomialSet(const PolynomialSet &p, bool newLine=false)=0;
  virtual void printPolynomialSetList(const PolynomialSetList &s)=0;
  virtual void printVector(const IntegerVector &v, bool curly=false, int minimalFieldWidth=0)=0;
  virtual void printFloatVector(FloatVector const &v, bool curly=false, int minimalFieldWidth=0)=0;
  virtual void printComplexVector(ComplexVector const &v, bool curly=false, int minimalFieldWidth=0)=0;
  virtual void printVectorList(const IntegerVectorList &s, bool indexed=false)=0;
  virtual void printString(const string &s)=0;
  virtual void printNewLine()=0;
  virtual void printChar(int c){char s[2];s[0]=c;s[1]=0;printString(s);}
  virtual void printTermOrder(TermOrder const &t);
  virtual void printPolyhedralCone(class PolyhedralCone const &c, bool xml=false);
  virtual void printPolyhedralFan(class PolyhedralFan const &c);
  virtual void printField(class Field const &f);
  virtual void printPolynomialRing(class PolynomialRing const &r);
  Printer& operator<<(IntegerVector const &v)
  {
    printVector(v);
    return *this;
  }
/*  Printer& operator<<(FloatVector const &v)
  {
    printFloatVector(v);
    return *this;
  }*/
  Printer& operator<<(class FieldElement const &v);
  Printer& operator<<(class FieldVector const &v);
  Printer& operator<<(class FieldMatrix const &m);
  Printer& operator<<(IntegerVectorList const &l);
  Printer& operator<<(PolynomialRing const &r);
  Printer& operator<<(PolynomialSet const &l);
  Printer& operator<<(Polynomial const &p);
  Printer& operator<<(const string &s);
  Printer& operator<<(list<int> &l);
  Printer& operator<<(class PolyhedralCone const &c);
  Printer& operator<<(class PolyhedralFan const &f);
  Printer& operator<<(class TermOrder const &t);
  Printer& operator<<(int a);
  Printer& operator<<(double a);
};


class LatexPrinter:public Printer
{
  int mathModeLevel;
  void pushMathMode();
  void popMathMode();
 public:
  LatexPrinter(FILE *f):Printer(f){mathModeLevel=0;}
  virtual void printInteger(int i, int minimalFieldWidth=0);
  virtual void printFloat(double i, int minimalFieldWidth=0);
  virtual void printComplexNumber(ComplexNumber const &i, int minimalFieldWidth=0);
  virtual void printMonomial(const Monomial &m, bool alwaysWriteSign=false, bool writeIfOne=true);
  virtual void printFieldElement(const FieldElement &e, bool writeIfOne=true, bool alwaysWriteSign=false);
  virtual void printTerm(const Term &t);
  virtual void printPolynomial(const Polynomial &p);
  virtual void printPolynomialSet(const PolynomialSet &s, bool newLine=false);
  virtual void printPolynomialSetList(const PolynomialSetList &s);
  virtual void printVector(const IntegerVector &v, bool curly=false, int minimalFieldWidth=0);
  virtual void printFloatVector(FloatVector const &v, bool curly=false, int minimalFieldWidth=0);
  virtual void printComplexVector(ComplexVector const &v, bool curly=false, int minimalFieldWidth=0);
  virtual void printVectorList(const IntegerVectorList &s, bool indexed=false);
  virtual void printString(const string &s);
  virtual void printNewLine();
  void printLatexStart();
  void printLatexEnd();
};


class AsciiPrinter:public Printer
{
 public:
  AsciiPrinter(FILE *f):Printer(f){}
  virtual char vectorLeftBrackets(){return '(';}
  virtual char vectorRightBrackets(){return ')';}
  virtual char vectorListLeftBrackets(){return '{';}
  virtual char vectorListRightBrackets(){return '}';}

  virtual void printInteger(int i, int minimalFieldWidth=0);
  virtual void printFloat(double i, int minimalFieldWidth=0);
  virtual void printComplexNumber(ComplexNumber const &i, int minimalFieldWidth=0);
  virtual void printMonomial(const Monomial &m, bool alwaysWriteSign=false, bool writeIfOne=true);
  virtual void printFieldElement(const FieldElement &e, bool writeIfOne=true, bool alwaysWriteSign=false);
  virtual void printTerm(const Term &t);
  virtual void printPolynomial(const Polynomial &p);
  virtual void printPolynomialSet(const PolynomialSet &s, bool newLine=false);
  virtual void printPolynomialSetList(const PolynomialSetList &s);
  virtual void printVector(const IntegerVector &v, bool curly=false, int minimalFieldWidth=0);
  virtual void printFloatVector(FloatVector const &v, bool curly=false, int minimalFieldWidth=0);
  virtual void printComplexVector(ComplexVector const &v, bool curly=false, int minimalFieldWidth=0);
  virtual void printVectorList(const IntegerVectorList &s, bool indexed=false);
  virtual void printString(const string &s);
  virtual void printNewLine();
};


class XmlPrinter:public Printer
{
 public:
  XmlPrinter(FILE *f):Printer(f){}

  virtual void printInteger(int i, int minimalFieldWidth=0);
  virtual void printFloat(double i, int minimalFieldWidth=0);
  virtual void printComplexNumber(ComplexNumber const &i, int minimalFieldWidth=0);
  virtual void printMonomial(const Monomial &m, bool alwaysWriteSign=false, bool writeIfOne=true);
  virtual void printFieldElement(const FieldElement &e, bool writeIfOne=true, bool alwaysWriteSign=false);
  virtual void printTerm(const Term &t);
  virtual void printPolynomial(const Polynomial &p);
  virtual void printPolynomialSet(const PolynomialSet &s, bool newLine=false);
  virtual void printPolynomialSetList(const PolynomialSetList &s);
  virtual void printString(const string &s);
  virtual void printNewLine();

  virtual void printVector(const IntegerVector &v, bool curly=false, int minimalFieldWidth=0);
  virtual void printFloatVector(FloatVector const &v, bool curly=false, int minimalFieldWidth=0);
  virtual void printComplexVector(ComplexVector const &v, bool curly=false, int minimalFieldWidth=0);
  virtual void printVectorList(const IntegerVectorList &s, bool indexed=false);
};



class TopcomPrinter:public AsciiPrinter
{
 public:
  TopcomPrinter(FILE *f):AsciiPrinter(f){}
  virtual char vectorLeftBrackets(){return '[';}
  virtual char vectorRightBrackets(){return ']';}
  virtual char vectorListLeftBrackets(){return '[';}
  virtual char vectorListRightBrackets(){return ']';}
};


extern AsciiPrinter debug;
extern AsciiPrinter pout;

#endif
