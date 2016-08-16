#include "printer.h"
#include "termorder.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "polynomialring.h"
#include "linalg.h"
#include "log.h"


//--------------------------------------------------
// Printer
//--------------------------------------------------

bool Printer::assertOnPrinting;


void Printer::printCheck()
{
  if(assertOnPrinting)assert(0);
}

void Printer::setAssertOnPrinting(bool b)
{
  assertOnPrinting=true;
}


void Printer::printTermOrder(TermOrder const &t)
{
  printCheck();
  t.print(*this);
}

string Printer::variableIndexToString(PolynomialRing const &r, int i)
{
  return r.getVariableName(i);
  //  fprintf(Stderr,"vits%i\n",i);
    /*  static char s[4];
  if(i>=0 && i<26)sprintf(s,"%c",i+'a');
  else if(i>=26 && i<52)sprintf(s,"%c",i+'A'-26);
  else assert(0);
  return string(s);
    */
}


void Printer::printVariable(PolynomialRing const &r, int i)
{
  printCheck();
  printString(variableIndexToString(r,i));
}

void Printer::printPolyhedralCone(class PolyhedralCone const &c, bool xml)
{
  printCheck();
  PolyhedralCone C=c;
  C.print(this,xml);
}


void Printer::printPolyhedralFan(class PolyhedralFan const &c)
{
  printCheck();
  c.print(this);
}


void Printer::printField(class Field const &f)
{
  printCheck();
  printString(f.toString().c_str());
}

void Printer::printPolynomialRing(class PolynomialRing const &r)
{
  printCheck();
  printField(r.getField());
  printString("[");
  printString(r.toStringVariableNames());
  printString("]");
}



Printer& Printer::operator<<(class FieldElement const &v)
{
	*this<<v.toString();
  return *this;
}


Printer& Printer::operator<<(class FieldVector const &v)
{
  v.print(*this);
  return *this;
}


Printer& Printer::operator<<(class FieldMatrix const &m)
{
  m.printMatrix(*this);
  return *this;
}


Printer& Printer::operator<<(IntegerVectorList const &l)
{
  printVectorList(l);
  return *this;
}


Printer& Printer::operator<<(PolynomialRing const &r)
{
  printPolynomialRing(r);
  return *this;
}


Printer& Printer::operator<<(PolynomialSet const &l)
{
  printPolynomialSet(l);
  return *this;
}


Printer& Printer::operator<<(Polynomial const &p)
{
  printPolynomial(p);
  return *this;
}


Printer& Printer::operator<<(const string &s)
{
  printString(s);
  return *this;
}


Printer& Printer::operator<<(list<int> &l)
{
  printString("{");
  for(list<int>::const_iterator i=l.begin();i!=l.end();i++)
    {
      if(i!=l.begin())printString(",");
      printInteger(*i);
    }
  printString("}");
  return *this;
}


Printer& Printer::operator<<(PolyhedralCone const &c)
{
  printPolyhedralCone(c);
  return *this;
}


Printer& Printer::operator<<(PolyhedralFan const &f)
{
  printPolyhedralFan(f);
  return *this;
}


Printer& Printer::operator<<(TermOrder const &t)
{
  t.print(*this);
  return *this;
}

Printer& Printer::operator<<(int a)
{
  printInteger(a);
  return *this;
}


Printer& Printer::operator<<(double a)
{
  printFloat(a);
  return *this;
}

//--------------------------------------------------
// LatexPrinter
//--------------------------------------------------

void LatexPrinter::pushMathMode()
{
  printCheck();
  if(mathModeLevel==0)fprintf(f," $ ");
  mathModeLevel++;
}


void LatexPrinter::popMathMode()
{
  printCheck();
  mathModeLevel--;
  if(mathModeLevel==0)fprintf(f," $ ");
}


void LatexPrinter::printMonomial(const Monomial &m, bool alwaysWriteSign, bool writeIfOne)
{
  printCheck();
  const IntegerVector &exponent=m.exponent;
  const int sign=1;
  pushMathMode();

  bool variablePrinted=false;
  for(int i=0;i<exponent.size();i++)if(exponent[i]*sign>0)
    {
      fprintf(f,"%s",variableIndexToString(m.getRing(),i).c_str());
      if(int(exponent[i]*sign)!=1)
	{
	  fprintf(f,"^{%i}",int(exponent[i]*sign));
	}
      variablePrinted=true;
    }
  if(!variablePrinted && writeIfOne)fprintf(f,"1");

  popMathMode();
}


void LatexPrinter::printInteger(int i, int minimalFieldWidth)
{
  printCheck();
  fprintf(f,"%i",i);
}


void LatexPrinter::printFloat(double i, int minimalFieldWidth)
{
  printCheck();
  assert(0);
}


void LatexPrinter::printComplexNumber(ComplexNumber const &i, int minimalFieldWidth)
{
  printCheck();
  assert(0);
}


void LatexPrinter::printFieldElement(const FieldElement &e, bool writeIfOne, bool alwaysWriteSign)
{
  printCheck();
  printString(e.toString(writeIfOne,alwaysWriteSign,true));
}

void LatexPrinter::printTerm(const Term &t)
{
  printCheck();

}


void LatexPrinter::printPolynomial(const Polynomial &p)
{
  printCheck();
  pushMathMode();

  bool first=true;
  // If the polynomial has a marked term it is written first
  IntegerVector e=p.getMarked().m.exponent;
  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    if(e==i->first.exponent)
      {
	printFieldElement(i->second,i->first.exponent.isZero(),!first);
	printMonomial(i->first,false,false);
	first=false;
      }
  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    if(e!=i->first.exponent)
      {
	printFieldElement(i->second,i->first.exponent.isZero(),!first);
	printMonomial(i->first,false,false);
	first=false;
      }
  popMathMode();
}

void LatexPrinter::printPolynomialSet(const PolynomialSet &s, bool newLine)
{
  printCheck();
  pushMathMode();
  printString("\\{");
  if(newLine)
    {
      popMathMode();
      printNewLine();
      pushMathMode();
    }
  for(PolynomialSet::const_iterator i=s.begin();i!=s.end();i++)
    {
      if(i!=s.begin())
        {
          printString(",\n");
          if(newLine)
            {
              popMathMode();
              printNewLine();
              pushMathMode();
            }
        }
      printPolynomial(*i);
    }
  printString("\\}\n");
  popMathMode();
}


void LatexPrinter::printPolynomialSetList(const PolynomialSetList &s)
{
  printCheck();
  pushMathMode();
  printString("\\{");
  for(PolynomialSetList::const_iterator i=s.begin();i!=s.end();i++)
    {
      if(i!=s.begin())printString(",\n");
      printPolynomialSet(*i);
    }
  printString("\\}\n");
  popMathMode();
}


void LatexPrinter::printVectorList(const IntegerVectorList &s, bool indexed)
{
  printCheck();
  int index=0;
  pushMathMode();
  printString("\\{");
  for(IntegerVectorList::const_iterator i=s.begin();i!=s.end();i++)
    {
      if(i!=s.begin())printString(",\n");
      if(indexed)
	{
	  printInteger(index++);
	  printString(": ");
	}
      printVector(*i);
    }
  printString("\\}\n");
  popMathMode();
}


void LatexPrinter::printVector(const IntegerVector &v, bool curly, int minimalFieldWidth)
{
  printCheck();
  int s=v.size();
  int i=0;

  fprintf(f,"$(");
  if(i<s)fprintf(f,"%i",(v[i]));
  for(i++;i<s;i++)fprintf(f,",%i",(v[i]));
  fprintf(f,")$\n");
}


void LatexPrinter::printFloatVector(const FloatVector &v, bool curly, int minimalFieldWidth)
{
  printCheck();
  assert(0);
}


void LatexPrinter::printComplexVector(const ComplexVector &v, bool curly, int minimalFieldWidth)
{
  printCheck();
  assert(0);
}


void LatexPrinter::printString(const string &s)
{
  printCheck();
  fprintf(f,"%s",s.c_str());
}


void LatexPrinter::printNewLine()
{
  printCheck();
  fprintf(f,"\n\n");
}


void LatexPrinter::printLatexStart()
{
  printCheck();
  fprintf(f,"\\documentclass[12pt,a4paper,dvips]{article}"
          "\\usepackage[english]{babel}"
          "\\usepackage{t1enc}"
          "\\usepackage{a4}"
          "\\usepackage{epsfig}"
          "\\usepackage{amsfonts}"
          "\\usepackage{latexsym}"
          "\\begin{document}\n");
}


void LatexPrinter::printLatexEnd()
{
  printCheck();
  fprintf(f,"\\end{document}");
}

//--------------------------------------------------
// AsciiPrinter
//--------------------------------------------------

void AsciiPrinter::printInteger(int i, int minimalFieldWidth)
{
  printCheck();
  static char s[16];
  sprintf(s,"%%%ii",minimalFieldWidth);
  fprintf(f,s,i);
}


void AsciiPrinter::printFloat(double i, int minimalFieldWidth)
{
	  printCheck();
	  static char s[32];
	  sprintf(s,"%%%if",minimalFieldWidth);
	  fprintf(f,s,i);
}


void AsciiPrinter::printComplexNumber(ComplexNumber const &i, int minimalFieldWidth)
{
	  printCheck();
	  FloatVector t(2);
	  t[0]=i.real();
	  t[1]=i.imag();
	  printFloatVector(t);
}


void AsciiPrinter::printMonomial(const Monomial &m, bool alwaysWriteSign, bool writeIfOne)
{//changed from binomial printing to Laurent printing April 2010
  printCheck();
  const IntegerVector &exponent=m.exponent;

  //  fprintf(Stderr,"[%i]",m.exponent.size());

  bool variablePrinted=false;
  for(int i=0;i<exponent.size();i++)if(exponent[i]!=0)
    {
      if(variablePrinted)fprintf(f,"*");
      fprintf(f,"%s",variableIndexToString(m.getRing(),i).c_str());
      //      fprintf(f,"%c",i+'a');
      if(int(exponent[i])!=1)
	{
    	  if(int(exponent[i])>0)
    		  fprintf(f,"^%i",int(exponent[i]));
    	  if(int(exponent[i])<0)
    		  fprintf(f,"^(%i)",int(exponent[i]));
	}
      variablePrinted=true;
    }
  if(!variablePrinted && writeIfOne)fprintf(f,"1");
  //  fprintf(stderr,"(%i)",exponent.size());
}


void AsciiPrinter::printFieldElement(const FieldElement &e, bool writeIfOne, bool alwaysWriteSign)
{
  printCheck();
  printString(e.toString(writeIfOne,alwaysWriteSign));
}

void AsciiPrinter::printTerm(const Term &t)
{
  printCheck();

}


void AsciiPrinter::printPolynomial(const Polynomial &p)
{
  printCheck();
  bool first=true;

  if(p.terms.empty())
    {
      printString("0");
      return;
    }
  // If the polynomial has a marked term it is written first
  //   printString("_");
  IntegerVector e=p.getMarked().m.exponent;
  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    if(e==i->first.exponent)
      {
	printFieldElement(i->second,i->first.exponent.isZero(),!first);
	if((!i->first.exponent.isZero())&&(!i->second.isOne())&&(!(-(i->second)).isOne()))printString("*");
	printMonomial(i->first,false,false);
	first=false;
      }
  //    printString("_");
  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    if(e!=i->first.exponent)
      {
	printFieldElement(i->second,i->first.exponent.isZero(),!first);
	if((!i->first.exponent.isZero())&&(!i->second.isOne())&&(!(-(i->second)).isOne()))printString("*");
	printMonomial(i->first,false,false);
	first=false;
      }
}


void AsciiPrinter::printPolynomialSet(const PolynomialSet &s, bool newLine)
{
  printCheck();
  printString("{");
  printNewLine();
  for(PolynomialSet::const_iterator i=s.begin();i!=s.end();i++)
    {
      if(i!=s.begin())printString(",\n");
      printPolynomial(*i);
    }
  printString("}\n");
}


void AsciiPrinter::printPolynomialSetList(const PolynomialSetList &s)
{
  printCheck();
  printString("{");
  printNewLine();
  for(PolynomialSetList::const_iterator i=s.begin();i!=s.end();i++)
    {
      if(i!=s.begin())printString(",\n");
      printPolynomialSet(*i);
    }
  printString("}\n");
}


void AsciiPrinter::printVectorList(const IntegerVectorList &s, bool indexed)
{
  printCheck();
  int index=0;
  printChar(vectorListLeftBrackets());
  printNewLine();
  for(IntegerVectorList::const_iterator i=s.begin();i!=s.end();i++)
    {
      if(i!=s.begin())printString(",\n");
      if(indexed)
	{
	  printInteger(index++);
	  printString(": ");
	}
      printVector(*i);
      //      printNewLine();
    }
  printChar(vectorListRightBrackets());
  printNewLine();
}


void AsciiPrinter::printVector(const IntegerVector &v, bool curly, int minimalFieldWidth)
{
  printCheck();
  fprintf(f,"%c",curly?'{':vectorLeftBrackets());
  for(int i=0;i<v.size();i++)
    {
      if(i!=0)fprintf(f,",");
      printInteger(v[i],minimalFieldWidth);
    }
  fprintf(f,"%c",curly?'}':vectorRightBrackets());
}


void AsciiPrinter::printFloatVector(const FloatVector &v, bool curly, int minimalFieldWidth)
{
  printCheck();
  fprintf(f,"%c",curly?'{':vectorLeftBrackets());
  for(int i=0;i<v.size();i++)
    {
      if(i!=0)fprintf(f,",");
      printFloat(v[i],minimalFieldWidth);
    }
  fprintf(f,"%c",curly?'}':vectorRightBrackets());
}


void AsciiPrinter::printComplexVector(const ComplexVector &v, bool curly, int minimalFieldWidth)
{
  printCheck();
  fprintf(f,"%c",curly?'{':vectorLeftBrackets());
  for(int i=0;i<v.size();i++)
    {
      if(i!=0)fprintf(f,",");
      printComplexNumber(v[i],minimalFieldWidth);
    }
  fprintf(f,"%c",curly?'}':vectorRightBrackets());
}


void AsciiPrinter::printString(const string &s)
{
  printCheck();
  fprintf(f,"%s",s.c_str());
}


void AsciiPrinter::printNewLine()
{
  printCheck();
  fprintf(f,"\n");
}


//--------------------------------------------------
// XmlPrinter
//--------------------------------------------------
void XmlPrinter::printVector(const IntegerVector &v, bool curly, int minimalFieldWidth)
{
  printCheck();
  int s=v.size();
  int i=0;

  fprintf(f,"<v>");
  if(i<s)fprintf(f,"%i",(v[i]));
  for(i++;i<s;i++)fprintf(f," %i",(v[i]));
  fprintf(f,"</v>\n");
}


void XmlPrinter::printFloatVector(const FloatVector &v, bool curly, int minimalFieldWidth)
{
  printCheck();
  assert(0);
}


void XmlPrinter::printComplexVector(ComplexVector const &v, bool curly, int minimalFieldWidth)
{
  printCheck();
  assert(0);
}


void XmlPrinter::printVectorList(const IntegerVectorList &s, bool indexed)
{
  printCheck();
  fprintf(f,"<m>");
  for(IntegerVectorList::const_iterator i=s.begin();i!=s.end();i++)
    printVector(*i);
  fprintf(f,"</m>\n");
}


void XmlPrinter::printInteger(int i, int minimalFieldWidth)
{
  printCheck();
  assert(0);
}


void XmlPrinter::printFloat(double i, int minimalFieldWidth)
{
  printCheck();
  assert(0);
}


void XmlPrinter::printComplexNumber(ComplexNumber const &i, int minimalFieldWidth)
{
  printCheck();
  assert(0);
}


void XmlPrinter::printMonomial(const Monomial &m, bool alwaysWriteSign, bool writeIfOne)
{
  printCheck();
  assert(0);
}


void XmlPrinter::printFieldElement(const FieldElement &e, bool writeIfOne, bool alwaysWriteSign)
{
  printCheck();
  assert(0);
}


void XmlPrinter::printTerm(const Term &t)
{
  printCheck();
  assert(0);
}


void XmlPrinter::printPolynomial(const Polynomial &p)
{
  printCheck();
  assert(0);
}


void XmlPrinter::printPolynomialSet(const PolynomialSet &s, bool newLine)
{
  printCheck();
  assert(0);
}


void XmlPrinter::printPolynomialSetList(const PolynomialSetList &s)
{
  printCheck();
  assert(0);
}


void XmlPrinter::printString(const string &s)
{
  printCheck();
  assert(0);
}


void XmlPrinter::printNewLine()
{
  printCheck();
  assert(0);
}

/*
void XmlPrinter::printVector(const IntegerVector &v)
{
}


void XmlPrinter::printString(const string &s)
{
}
*/



AsciiPrinter debug(Stderr);
AsciiPrinter pout(Stdout);
