#ifndef SYMMETRY_H_INCLUDED
#define SYMMETRY_H_INCLUDED

#include <stdio.h>
#include <set>
#include "vektor.h"
#include "polynomial.h"
#include "termorder.h"

class SymmetryGroup{
  unsigned char *byteTable;
  int byteTableHeight;
  class Trie *trie;
 public:
  typedef set<IntegerVector,LexicographicTermOrder> ElementContainer;
  ElementContainer elements;
  int sizeOfBaseSet()const;
  /**
     Applies the permutation v to each vector in the list l.
   */
  static IntegerVectorList permuteIntegerVectorList(IntegerVectorList const &l, IntegerVector const &v);
  static Polynomial permutePolynomial(const Polynomial &p, IntegerVector const &v);
  static PolynomialSet permutePolynomialSet(PolynomialSet const &s, IntegerVector const &v);
  Polynomial computeUniqueRepresentative(Polynomial p);
  static IntegerVector compose(IntegerVector const &perm, IntegerVector const &b);
  static IntegerVector composeInverse(IntegerVector const &a, IntegerVector const &b);
  static IntegerVector identity(int numberOfElements);
  static IntegerVector inverse(IntegerVector const &a);
  static bool isPermutation(IntegerVector const &a);
/**
 * Takes a permutation of variables and sign changes of variables and combines them in to a single permutation of the variables and their inverses. That is, the permutation is a permutation on 2n elements.
 */
  static IntegerVector combinePermutationAndSignChanges(IntegerVector const &permutation, IntegerVector const &signChanges);
/**
 * Does the opposite transformation as combinePermutationAndSignChanges.
 */
  static void extractPermuationAndSignChanges(IntegerVector const &v, IntegerVector &permutation, IntegerVector &signChanges);
  /**
     The set of vectors which are not improved lexicographically when
     perm is applied to them is convex. Its closure is a
     halfspace. This routine returns the inner normal of this
     halfspace. The only exception is if perm is the identity then the
     zero vector is returned.
   */
  static IntegerVector fundamentalDomainInequality(IntegerVector const &perm);
  /**
     The set of vectors which cannot be improved lexicographically by
     applying an element from the group is a convex set. Its closure
     is a polyhedral cone. This routine returns a set of inequalities
     The returned list does not contain the zero vector.
   */
  IntegerVectorList fundamentalDomainInequalities()const;
  SymmetryGroup(int n);
  void computeClosure(IntegerVector const &v);
  void computeClosure(IntegerVectorList const &l);
  IntegerVectorList getUniqueGenerators()const;
  int orbitSize(IntegerVector const &stable)const;
  bool isTrivial()const;
  void print(FILE *f);
  /**
     The symmetry group acts on vectors by permuting the entries. The
     following routine returns a unique representative for the orbit
     containing v. This makes it easy to check if two elements are in
     the same orbit. The permutation used to get this representative
     is stored in *usedPermutation (if pointer not 0).
   */
  IntegerVector orbitRepresentative(IntegerVector const &v, IntegerVector *usedPermutation=0)const;
  /**
     This routine works as orbitRepresentative() except that the
     symmetry group considered is only the subgroup keeping the vector
     fixed fixed.
   */
  IntegerVector orbitRepresentativeFixing(IntegerVector const &v, IntegerVector const &fixed)const;

  // Methods for highly optimized symmetry group computations:
  void createByteTable();//Can only be called once. SymmetryGroup is not allowed to be changed afterwards or to be copied. Leaks memory at destruction.
  void createTrie();
  unsigned char *getByteTable()const;
  int getByteTableHeight()const;
};

/**
 * Sorts v and returns the number of swaps performed.
 */
int mergeSort(IntegerVector &v);

#endif
