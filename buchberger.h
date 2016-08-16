#ifndef BUCHBERGER_H_INCLUDED
#define BUCHBERGER_H_INCLUDED

#include "polynomial.h"

Polynomial sPolynomial(Polynomial a, Polynomial b);
void buchberger(PolynomialSet *g, TermOrder const &termOrder, bool allowSaturation=false);
void minimize(PolynomialSet *g);
/**
 * The Groebner basis g must be minimal before this function is called.
 */
void autoReduce(PolynomialSet *g, TermOrder const &termOrder);
static inline void autoReduce_(PolynomialSet *g, TermOrder const &termOrder){return autoReduce(g,termOrder);}//<----avoiding scoperules. Should start using namespaces.
bool isMarkedGroebnerBasis(PolynomialSet const &g);

/* For the autoReduction procedure the TermOrder argument is only used as an
argument for the division algorithm. This means that the input to autoReduce
should be marked and that the term order does not have to agree.
*/

/** This routine takes a marked Groebner basis and checks if it is minimal. TODO: there is room for improvement here. */
bool isMinimal(PolynomialSet const &g);
/** This routine takes a marekd Groebner basis and checks if it is autoreduced. This means that it also checks if the basis is minimal.*/
bool isReduced(PolynomialSet const &g);

#endif
