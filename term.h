#ifndef TERM_H_INCLUDED
#define TERM_H_INCLUDED

#include "field.h"
#include "monomial.h"
#include "polynomialring.h"

class Term
{
  PolynomialRing theRing;
 public:
  Monomial m;
  FieldElement c;
  Term(PolynomialRing const &r):theRing(r),m(r),c(r.getField()){};
    Term(FieldElement const &c_, Monomial const &m_);
  void operator*=(const Term &t);
  PolynomialRing const &getRing()const{return theRing;}
  Term operator/(const Term &b)const;
  Term operator/(const Monomial &b)const;
};

#endif
