/*
 * integergb.h
 *
 *  Created on: Dec 14, 2010
 *      Author: anders
 */

#ifndef INTEGERGB_H_
#define INTEGERGB_H_

/*
 * integergb.cpp
 *
 *  Created on: Dec 14, 2010
 *      Author: anders
 */

#include "polynomial.h"
#include "field_rationals.h"
#include "symmetrictraversal.h"

/*
 * Implemented according to [Becker, Weispfenning].
 */

Polynomial dDivision(Polynomial p, PolynomialSet const &l, TermOrder const &termOrder);
Polynomial spol(Polynomial const &g1, Polynomial const &g2);
Polynomial gpol(Polynomial const &g1, Polynomial const &g2);
void zAutoReduce(PolynomialSet *g, TermOrder const &termOrder);
void zBuchberger(PolynomialSet &F, TermOrder const &T);

class IntegerGroebnerFanTraverser: public ConeTraverser
{
  PolynomialSet groebnerBasis;
  PolyhedralCone theCone;
  int n;//,d;
  int prime;
  void updatePolyhedralCone();
public:
  IntegerGroebnerFanTraverser(PolynomialSet const &generators);
  virtual void changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector);
  virtual IntegerVectorList link(IntegerVector const &ridgeVector);
  PolyhedralCone & refToPolyhedralCone();
//  PolynomialSet &refToGroebnerBasisRepresentation();
//  PolynomialSet initialIdeal()const;
};

#endif /* INTEGERGB_H_ */
