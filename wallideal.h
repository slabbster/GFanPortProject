#ifndef WALLIDEAL_H_INCLUDED
#define WALLIDEAL_H_INCLUDED

#include "polynomial.h"
#include "polyhedralcone.h"

/* The correct way to get a set of defining inequalities for a
   full-dimensional groebner cone is by calling wallInequalities() and
   possibly algebraicTest(). It is important not to use
   wallFlipableNormals() since this will not work for non-homogeneous
   ideals. The algebraicTest() however works fine in the
   non-homogeneous case. */

IntegerVectorList algebraicTest(IntegerVectorList const &l, PolynomialSet const &groebnerBasis);//run the list through the algebraic test
Polynomial wallPolynomial(Polynomial const &p, IntegerVector const &wallNormal);
PolynomialSet wallIdeal(PolynomialSet const &groebnerBasis, IntegerVector const &wallNormal);
PolynomialSet lowerDimensionalWallIdeal(PolynomialSet const &groebnerBasis, IntegerVectorList const &wallEqualities);
//IntegerVectorList wallNormals(PolynomialSet const &groebnerBasis);
IntegerVectorList wallRemoveScaledInequalities(IntegerVectorList const &l);
/**
 * This routine takes a list of marked polynomials. For each polynomial it takes the
 * exponent vector of the marked term and for every other term subtracts the exponent
 * vector of that term and stores the difference in a set which is output at the end.
 * Duplicates are removed.
 */
IntegerVectorList exponentDifferences(PolynomialSet const &g);
/**
 * This routine takes a list g of marked polynomials and computes a list of
 * strict inequalities that are as set are satisfied for a weight vector w
 * if and only if w induces the same marking on g.
 */
IntegerVectorList wallInequalities(PolynomialSet const &g);
bool wallContainsPositiveVector(IntegerVector const &wallNormal);
PolynomialSet flip(PolynomialSet const &groebnerBasis, IntegerVector const &wallNormal, TermOrder *autoReduceHint=0);
PolynomialSet flipMinkowski(PolynomialSet const &groebnerBasis, IntegerVector const &wallNormal);
void wallAddCoordinateWalls(IntegerVectorList &normals);
IntegerVectorList wallFlipableNormals(PolynomialSet const &groebnerBasis, bool isKnownToBeHomogeneous=false);
bool isIdealHomogeneous(PolynomialSet const &groebnerBasis); //wrt a positive vector
PolyhedralCone homogeneitySpace(PolynomialSet const &reducedGroebnerBasis);
PolyhedralCone groebnerCone(PolynomialSet const &reducedGroebnerBasis, bool useAlgebraicTest);
int dimensionOfHomogeneitySpace(PolynomialSet const &reducedGroebnerBasis);
PolynomialSet liftBasis(PolynomialSet const &toBeLifted, PolynomialSet const &originalBasisForFullIdeal);
IntegerVectorList normalizedWithSumsAndDuplicatesRemoved(IntegerVectorList const &a);
bool fastIsFacet(IntegerVectorList const &normals, IntegerVectorList::const_iterator i);//assumes that normals contains no duplicates (or scaled vectors) and that the cone is full-dimensional
IntegerVectorList fastNormals(IntegerVectorList const &inequalities);//using the two above routines. Assuming full-dimensionality

bool isMarkingConsistent(PolynomialSet const &g);

/**
   This routine takes a set of polynomials (not necessarily marked)
   and computes the smallest cone in the normal fan of the Minkowski
   sum of their Newton polytopes containing the vector w. In
   particular this routine is useful for computing Groebner cones.
   The returned cone is in canonical form.
 */
PolyhedralCone normalConeOfMinkowskiSum(PolynomialSet const &polynomials, IntegerVector const &w);

/**
   This routine returns a list of equations cutting out a space. This
   space is the smallest subspace of gradings for which all
   polynomials in the parameter set are homogeneous.
 */
IntegerVectorList commonHomogeneitySpaceEquations(PolynomialSet const &polynomials);

/**
   This routine returns a basis a space. This
   space is the smallest subspace of gradings for which all
   polynomials in the parameter set are homogeneous.
 */
IntegerVectorList commonHomogeneitySpaceGenerators(PolynomialSet const &polynomials);

#endif
