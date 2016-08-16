#ifndef POLYHEDRALCONE_H_INCLUDED
#define POLYHEDRALCONE_H_INCLUDED

#include "vektor.h"
#include "printer.h"

/**
A PolyhedralCone is represented by halfSpaces and equations. The halfSpaces is a list of non-strict inequalities while equations is a list of equations.

A cone can be in one of the four states:
0) Nothing has been done to remove redundancies. This is the initial state.
1) A basis for the true, implied equations space has been computed. This means that the implied equations have been computed. In particular the dimension of the cone is known.
2) Redundant inequalities have been computed and have been eliminated. This means that the true set of facets is known - one for each element in halfSpaces.
3) The inequalities and equations from 2) have been transformed into a canonical form. Besides having a unique representation for the cone this also allows comparisons between cones with operator<().

It is not clear have this class is designed optimally. Should the user of this class know about the states?

Always putting the cone in state 1 after something has changed helps a lot. Then all operations can be performed except comparing and getting facets with out taking the cone to a special state.

Old explanation

 The idea is that equations is always updated by calling computeAndReduceLinearitySpace() when something is changed. If an element in halfSpaces is implied by equations it might be removed during this process. The set of halfspaces is not minimised until reduce() is called. The reason for this is that most functions do not require the heavy reduction computation. canonicalize() will (when implemented correctly) make sure that the defining set of equalities and inequalities depends deterministically on the actual set of points in the cone. The isCanonicalized variable will be set to true after this process, and polyhedral cones with this flag set can be ordered with operator<().

Under the assumption that computeAndReducedLinearitySpace is implemented correctly (and int is 32 and long long 64) getRelativeInteriorPoint() is guaranteed to produce correct output if it returns, otherwise it will assert(0).

The ordering will make high-dimensional cones largest.
 */

enum PolyhedralConePreassumptions{
  PCP_none=0,
  PCP_impliedEquationsKnown=1,
  PCP_facetsKnown=2
};

class PolyhedralCone
{
  int preassumptions;
  int state;
  int n;
  IntegerVector linearForm;//Associate a linear form to the cone. Useful for specifying tropical rational functions in terms of their breaking fan. The value does not have to be set
  IntegerVectorList inequalities;
  IntegerVectorList equations;
  mutable IntegerVectorList cachedExtremeRays;
/**
 * If this bool is true it means that cachedExtremeRays contains the extreme rays as found by extremeRays().
 */
  mutable bool haveExtremeRaysBeenCached;
  void ensureStateAsMinimum(int s);
  int multiplicity;//This really should not be a part of the PolyhedralCone class. But since Gfan passes parameters by value putting it here makes things easier
//  * (This version of the constructor should be called if the state of the generated cone is known a priori.)
//  PolyhedralCone(IntegerVectorList const &halfSpaces_, IntegerVectorList const &equations_, int ambientDimension, int state);
 public:
   /**
    * Starting as a test Aug 2010, the user can give information about what is known by specifying a PolyhedralConePreassumptions flag.
    */
     PolyhedralCone(IntegerVectorList const &halfSpaces_, IntegerVectorList const &equations_, int ambientDimension=-1, int preassumptions_=PCP_none);
	 static PolyhedralCone polyhedralConeWithKnownImpliedEquations(IntegerVectorList const &halfSpaces, IntegerVectorList const &equations, int ambientDimension)
	 {
//           return PolyhedralCone(halfSpaces,equations,ambientDimension,1);//CHANGEBACK
           return PolyhedralCone(halfSpaces,equations,ambientDimension,PCP_impliedEquationsKnown);
	 }

	 IntegerVector const &getLinearForm()const;
  void setLinearForm(IntegerVector const &linearForm_);
  bool isInStateMinimum(int s)const;
  int getState()const;
  IntegerVectorList getHalfSpaces()const;
  const IntegerVectorList &getEquations()const;
  IntegerVectorList generatorsOfSpan()const;
  IntegerVectorList generatorsOfLinealitySpace()const;

  bool areFacetsKnown()const{return (state>=2)||(preassumptions&PCP_facetsKnown);}
  bool areImpliedEquationsKnown()const{return (state>=1)||(preassumptions&PCP_impliedEquationsKnown);}

  void canonicalize();
  void findFacets();
  /**
   * After this function has been called a minimal set of implied equations for the cone is known and is
   * returned when calling getEquations(). The returned equations form a basis of the space of implied
   * equations.
   */
  void findImpliedEquations();

  PolyhedralCone(int ambientDimension=0);
//  PolyhedralCone(IntegerVectorList const &halfSpaces_, IntegerVectorList const &equations_, int ambientDimension=-1);
  IntegerVector getRelativeInteriorPoint()const;
  /**
     Assuming that this cone C is in state at least 3 (why not 2?), this routine returns a relative interior point v(C) of C with the following properties:
     1) v is a function, that is v(C) is found deterministically
     2) for any angle preserving, lattice preserving and lineality space preserving transformation T of R^n we have that v(T(C))=T(v(C)). This makes it easy to check if two cones in the same fan are equal up to symmetry. Here preserving the lineality space L just means T(L)=L.
  */
  IntegerVector getUniquePoint()const;
  /**
   * Takes a list of possible extreme rays and add up those actually contained in the cone.
   */
  IntegerVector getUniquePointFromExtremeRays(IntegerVectorList const &extremeRays)const;
  int ambientDimension()const;
  int dimension()const;
  int codimension()const;//dimension-ambientDimension
  int dimensionOfLinealitySpace()const;
  bool isZero()const;
  bool isFullSpace()const;
  friend PolyhedralCone intersection(const PolyhedralCone &a, const PolyhedralCone &b);
  friend PolyhedralCone product(const PolyhedralCone &a, const PolyhedralCone &b);
  static PolyhedralCone positiveOrthant(int dimension);
  static PolyhedralCone givenByRays(IntegerVectorList const &generators, IntegerVectorList const &linealitySpace, int n);
  void print(class Printer *p, bool xml=false)const;
  PolyhedralCone withLastCoordinateRemoved()const;
  /**
   * Largest dimensional cones are smallest.
   */
  friend bool operator<(PolyhedralCone const &a, PolyhedralCone const &b);
  friend bool operator!=(PolyhedralCone const &a, PolyhedralCone const &b);
  bool containsPositiveVector()const;
  bool contains(IntegerVector const &v)const;
  bool contains(IntegerVectorList const &l)const;
  bool contains(PolyhedralCone const &c)const;
  /**
   * Returns true if the PolyhedralCone contains v in its relative interior. False otherwise. The cone must be in state at least 1.
   */
  bool containsRelatively(IntegerVector const &v)const;
  bool isSimplicial()const;
  PolyhedralCone permuted(IntegerVector const &v)const;
  PolyhedralCone linealitySpace()const;
  PolyhedralCone dualCone()const;//only implemented for subspaces
  PolyhedralCone negated()const;
  /**
   * The returned list of extreme rays are orthogonal to the lineality space and primitive.
   * This makes them invariant under lattice preserving linear transformations.
   * If generators for the lineality space are known, they can be supplied. This can
   * speed up computations a lot.
   */
  IntegerVectorList extremeRays(IntegerVectorList const *generatorsOfLinealitySpace=0)const;
  bool checkDual(PolyhedralCone const &c)const;//check some necessary conditions
  int getMultiplicity()const;
  void setMultiplicity(int m);
  /**
     The cone defines two lattices, namely Z^n intersected with the
     span of the cone and Z^n intersected with the lineality space of
     the cone. Clearly the second is contained in the
     first. Furthermore, the second is a saturated lattice of the
     first. The quotient is torsion-free - hence a lattice. Generators
     of this lattice as vectors in the span of the cone are computed
     by this routine. The implied equations must be known when this
     function is called - if not the routine asserts.
   */
  IntegerVectorList quotientLatticeBasis()const;
  /**
     For a ray (dimension of lineality space +1 equals cone dimension)
     the quotent lattice described in quotientLatticeBasis() is
     isomorphic to Z. In fact the ray intersected with Z^n modulo the
     lineality space intersected with Z^n is a semigroup generated by
     just one element. This routine computes that element as an
     integer vector in the cone. Asserts if the cone is not a ray.
     Asserts if the implied equations have not been computed.
   */
  IntegerVector semiGroupGeneratorOfRay()const;

  /**
     Computes the link of the face containing v in its relative
     interior.
   */
  PolyhedralCone link(IntegerVector const &w)const;

  /**
     Computes the volume of the cone intersected with the cube with
     edge length two centered at the origin.
   */
  FieldElement volume()const;

  /**
     Tests if f is a face of the cone.
   */
  bool hasFace(PolyhedralCone const &f)const;
/**
 Computes the face of the cone containing v in its relative interior.
 */
  PolyhedralCone faceContaining(IntegerVector const &v)const;
  /**
   * Computes the projection of the cone to the first newn coordinates.
   * The ambient space of the returned cone has dimension newn.
   */
  PolyhedralCone projection(int newn)const;
};

#endif
