/*
 * padic.h
 *
 *  Created on: Dec 1, 2010
 *      Author: anders
 */

#ifndef PADIC_H_
#define PADIC_H_

#include "polynomial.h"
#include "symmetrictraversal.h"

Polynomial longDivision(Polynomial f, PolynomialSet const &G2, int prime, IntegerVector const &omega, TermOrder const &tieBreaker, PolynomialSet &H2, Polynomial &u2);
Polynomial SPolynomial(Polynomial const &a, Polynomial const &b, int prime, IntegerVector const &omega, TermOrder const &tieBreaker);
Polynomial pAdicInitialForm(PolynomialRing const &ZModPZRing, Polynomial const &f, int prime, IntegerVector const &omega);
PolynomialSet pAdicInitialForms(PolynomialRing const &ZModPZRing, PolynomialSet const &g, int prime, IntegerVector const &omega);
/**
 * Breaks ties.
 */
Polynomial pAdicInitialTerms(PolynomialRing const &ZModPZRing, Polynomial const &f, int prime, IntegerVector const &omega, TermOrder const &tieBreaker);
PolynomialSet pAdicInitialTerms(PolynomialRing const &ZModPZRing, PolynomialSet const &g, int prime, IntegerVector const &omega, TermOrder const &tieBreaker);
PolynomialRing residuePolynomialRing(PolynomialRing const &R, int prime);

IntegerVectorList normalPolyhedralInequalities(Polynomial const &f, int prime, IntegerVector const &omega, TermOrder const &tieBreaker);
IntegerVectorList normalPolyhedralInequalities(PolynomialSet const &g, int prime, IntegerVector const &omega, TermOrder const &tieBreaker);

void pAdicAutoReduce(PolynomialSet &g, int prime, IntegerVector const &omega, TermOrder const &tieBreaker);
void pAdicBuchberger(PolynomialSet &g, int prime, IntegerVector const &omega, TermOrder const &tieBreaker);

class PAdicGroebnerFanTraverser: public ConeTraverser
{
  PolynomialSet groebnerBasis;
  PolyhedralCone theCone;
  int n;//,d;
  int prime;
  void updatePolyhedralCone(IntegerVector const &omega, TermOrder const &tieBreaker);
public:
  PAdicGroebnerFanTraverser(PolynomialSet const &generators, int prime_);
  virtual void changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector);
  virtual IntegerVectorList link(IntegerVector const &ridgeVector);
  PolyhedralCone & refToPolyhedralCone();
//  PolynomialSet &refToGroebnerBasisRepresentation();
//  PolynomialSet initialIdeal()const;
};


#endif /* PADIC_H_ */
