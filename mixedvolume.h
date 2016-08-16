#ifndef MIXEDVOLUME_H_INCLUDED
#define MIXEDVOLUME_H_INCLUDED

#include "polynomial.h"

/**
 * Computes the mixed volume of the Newton polytopes in g.
 */
int64 mixedVolume(PolynomialSet const &g);

/**
 * Returns true iff the mixed volume of the Newton polytopes in g is nonzero.
 */
bool mixedVolumePositive(PolynomialSet const &g);

#endif
