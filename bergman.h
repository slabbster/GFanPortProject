#ifndef BERGMAN_H_INCLUDED
#define BERGMAN_H_INCLUDED

#include "polynomial.h"
#include "printer.h"
#include "symmetry.h"
#include "polyhedralfan.h"

/**
   Represents a tropical variety up to symmetry.

   New strategy implemented January 2009: We should only store the
   idealGroebnerBasis when actually needed. This means that for
   computing tropical varieties of curves it is stored, while for
   tropical varieties of prime ideals it should not be stored.
   Instead of storing the Groebner basis a polyhedral cone is stored.
 */

class BergmanFan
{
  bool simplicial;
 public:
  SymmetryGroup symmetryGroup;
  class CodimensionOneCone
    {
    public:
      list<int> incidenceList;
      IntegerVectorList incidencePermutationList;
      PolynomialSet idealGroebnerBasis;
      CodimensionOneCone(PolynomialSet const &idealGroebnerBasis_):
	idealGroebnerBasis(idealGroebnerBasis_)
	{
	}
    };
  class MaximalCone
  {
  public:
    bool idealBasisStored;
    int label;
    int numberOfFacets;
    PolynomialSet coneGroebnerBasis;
    PolynomialSet idealGroebnerBasis;
    PolyhedralCone theCone;
    int multiplicity;
    MaximalCone(PolynomialSet const &coneGroebnerBasis_, PolynomialSet const &idealGroebnerBasis_, bool storeIdealBasis, int label_, int numberOfFacets_=-1);
  };
  void print(Printer &p);
  bool contains(PolynomialSet const &g);
  typedef list<MaximalCone> MaximalConeList;
  typedef list<CodimensionOneCone> CodimensionOneConeList;

  BergmanFan():
    symmetryGroup(0),
    simplicial(false)
    {
    }
  MaximalConeList cones;
  CodimensionOneConeList codimensionOneCones;
  void setSymmetryGroup(SymmetryGroup const &s);
  PolyhedralFan toPolyhedralFan()const;
  int numberOfMaximalCones()const;
  void setSimplicial(bool b);
  bool isSimplicial()const;
  void computeMultiplicities();
};

BergmanFan bergmanRay(PolynomialSet const &idealGroebnerBasis);

/**
 * This routine computes a tropical curve. The input is a monomial-free
 * homogeneous (in some positive grading) ideal with dimension equal to
 * 1 plus the dimension of the homogeneity space.
 * The output is a BergmanFan object with one MaximalCone for each ray of the curve.
 */
BergmanFan bergmanRayIntersection(PolynomialSet const &idealGroebnerBasis);
BergmanFan bergman(PolynomialSet const &coneGroebnerBasis, PolynomialSet const &idealGroebnerBasis, SymmetryGroup const *symmetryGroup=0);

#endif
