#ifndef TROPICAL2_H_INCLUDED
#define TROPICAL2_H_INCLUDED

#include "polyhedralfan.h"
#include "polynomial.h"

/* This file contains tropical procedures implemented after April 1st 2005 */

Polynomial initialForm(Polynomial const &p, IntegerVector const &weight);
PolynomialSet initialForms(PolynomialSet const &groebnerBasis, IntegerVector const &weight);
PolynomialSet initialIdeal(PolynomialSet const &g, IntegerVector const &weight);//Assume homogeneous
Polynomial initialFormAssumeMarked(Polynomial const &p, IntegerVector const &weight);
PolynomialSet initialFormsAssumeMarked(PolynomialSet const &groebnerBasis, IntegerVector const &weight);
PolyhedralFan tropicalPrincipalIntersection(int n, PolynomialSet const &g, int linealitySpaceDimension=-1);
PolynomialSet guessInitialIdealWithoutMonomial(PolynomialSet const &groebnerBasis, PolynomialSet *fullNeighbourBasis, bool onlyCheckRays); //ideal must be homogeneous
PolynomialSet guessInitialIdealWithoutMonomialStably(PolynomialSet const &groebnerBasis, PolynomialSet *fullNeighbourBasis, bool onlyCheckRays);
  // fullNeighbourBasis is set to a Groebner basis of the full ideal. The returned basis and fullNeighbourBasis have at least one termorder in common

#endif
