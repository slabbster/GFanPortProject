#ifndef SUBSPACE_H_INCLUDED
#define SUBSPACE_H_INCLUDED

#include "vektor.h"
#include "polynomial.h"
#include "linalg.h"

/**
   @brief A subspace of Q^n.

   An object of this class represents a linear subspace of Q^n. All
   functionality is already available in the class FieldMatrix.
   However, the Subspace class is different in two ways: a Subspace
   cannot be changed after construction and the Subspace class
   provides operator< for sorting and distinguishing mathematical
   subspaces.
*/
class Subspace
{
  int n;//ambient dimension

  //  PolynomialSet basis;
  FieldMatrix basis2;
  IntegerVectorList integerRep;//REMOVE THIS
 public:
  /**
     Constructs a subspace given a generating set.

     @param generators generators for the subspace.
     @param ambientDimension the dimension of the ambient space.
     This number must equal the sizes of the vectors in generators.
   */
  Subspace(IntegerVectorList const &generators, int ambientDimension=-1);
  bool contains(IntegerVector const &v)const;
  int dimension();
  friend Subspace sum(Subspace const &a, Subspace const &b);
  int ambientDimension()const;
  IntegerVectorList getRepresentation()const;
  IntegerVector canonicalizeVector(IntegerVector const &v)const;
  bool operator<(Subspace const &b)const;
};

#endif
