#ifndef FIELDLP_H_INCLUDED
#define FIELDLP_H_INCLUDED

#include "linalg.h"
#include <set>

/**
   @brief A linear programming problem.

   This class represents linear programming problems of the form
   max(w,x) subject to Ax<=b. Or rather subject to
   A_ix<=b_i + epsilon^i. Notice that such perturbed system, if feassible,
   is fulldimensional. The problems can be stated over any ordered Field.

   The class has a state represented by the basis basis and vector x which first must be set by calling setCurrentBasis() or...
   A basis can only be set if the system is feassible.
   Afterwards step() can be called repeatedly until an optimal solution is found or the problem turns out to be unbounded.

	TODO: Anti-cycling rule goes wrong in this implementation. Fix this.
*/
class FieldLP
{
  Field theField;
  FieldMatrix A;
  FieldVector b;
  FieldVector w;
  FieldVector x;
  set<int> basis;

  static FieldVector subvector(FieldVector const &v, set<int> const &b);
  static FieldMatrix submatrix(FieldMatrix const &m, set<int> const &b);

  FieldVector edgeDirection(int i)const;
  bool isImprovingDirection(int i)const;
  FieldElement improvement(int i, int &newBasisMember)const;
 public:
  FieldLP(FieldMatrix const &A_, FieldVector const &b_);
  FieldVector basisToPoint(set<int> const &basis)const;
  void setObjectiveFunction(FieldVector const &w_);
  void setCurrentBasis(set<int> const &basis_);
  FieldMatrix computeLinealitySpace();
  FieldLP buildLPForFeasibilityCheck();
  void print(Printer &P)const;
  int step();
  bool findFeasibleBasis();
  FieldLP withNoLineality();
};

#endif
