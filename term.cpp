#include "term.h"

#include "printer.h"

Term::Term(FieldElement const &c_, Monomial const &m_):m(m_),c(c_),theRing(m_.getRing())
{
  //  fprintf(Stderr,"TermConstructor :");
  //  AsciiPrinter(Stderr).printFieldElement(c);
  //  AsciiPrinter(Stderr).printMonomial(m);
}

void Term::operator*=(const Term &t)
{
  m.exponent+=t.m.exponent;
  c*=t.c;
}

Term Term::operator/(const Term &b)const
{
  return Term(c*b.c.inverse(),Monomial(theRing,m.exponent-b.m.exponent));
}

Term Term::operator/(const Monomial &b)const
{
  return Term(c,Monomial(theRing,m.exponent-b.exponent));
}

