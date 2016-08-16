#ifndef POLYHEDRALFAN_H_INCLUDED
#define POLYHEDRALFAN_H_INCLUDED

#include <set>
#include <map>

#include "polyhedralcone.h"
#include "polynomial.h"
#include "printer.h"

class SymmetryGroup;

typedef set<PolyhedralCone> PolyhedralConeList;

typedef map<int,IntegerVectorList> IncidenceList;

enum FanPrintingFlags{
  FPF_conesCompressed=1,
  FPF_conesExpanded=2,
  FPF_cones=4,
  FPF_maximalCones=8,
  FPF_boundedInfo=16,
  FPF_values=32,
  FPF_group=64,
  FPF_multiplicities=128,
  FPF_xml=256,
  FPF_tPlaneSort=512,
  FPF_primitiveRays=1024,

  FPF_default=2+4+8
};


/** A PolyhedralFan is simply a collection of canonicalized PolyhedralCones.
 * It contains no combinatorial information in the sense of a polyhedral complex.
 * A cone being present in the PolyhedralFan corresponds to the cone and all its facets being present
 * in the mathematical object.
 * The intersection of cones in the fan must be a face of both.
 * In particular all cones in a PolyhedralFan have the same lineality space.*/
class PolyhedralFan
{
  int n;
  PolyhedralConeList cones;
 public:
  static IntegerVector stableRay(PolyhedralCone const &c, SymmetryGroup const *sym);
  static IntegerVector stableRay(IntegerVector const &v, IntegerVectorList const &equations, IntegerVectorList const &inequalities, SymmetryGroup const *sym);
  static class PolyhedralFan bergmanOfPrincipalIdeal(Polynomial const &p);
  static class PolyhedralFan normalFanOfNewtonPolytope(Polynomial const &p);
  static class PolyhedralFan fullSpace(int n);
  static class PolyhedralFan halfSpace(int n, int i);
  static class PolyhedralFan facetsOfCone(PolyhedralCone const &c);
  /**
     This routine computes a fan F with minimal support with the
     following property: the cone can be inserted to (is compatible
     with) F and this will make F complete. The support of F is the
     closure of the complement of C. If includec is true then C is
     indeed inserted into F before returning.
   */
  static class PolyhedralFan complementOfCone(PolyhedralCone const &c, bool includec=false);
  PolyhedralFan(int ambientDimension);
  void print(class Printer *p)const;
  //void printWithIndices(class Printer *p, SymmetryGroup *sym=0, const char *polymakeFileName=0)const; //fan must be pure
  vector<string> renamingStrings(IntegerVectorList const &theVectors, IntegerVectorList const &originalRays, IntegerVectorList const &linealitySpace, SymmetryGroup *sym)const;
  void printWithIndices(class Printer *p, int flags=FPF_default, SymmetryGroup *sym=0, vector<string> const *comments=0)const;
  //  void printWithIndices(class Printer *p, bool printMultiplicities, SymmetryGroup *sym, bool group=false, bool ignoreCones=false, bool xml=false, bool tPlaneSort=false, vector<string> const *comments=0)const;
  /* Read in a polyhedral fan, but with the cones containing w.  If
     present, only read in cones among coneIndices.  If sym is
     present, read COMPRESSED section and make w containment up to
     symmetry, taking all elements in the orbit that contains w into
     the fan.  If onlyMaximal is set then only maximal cones are read
     in.
   */
  static PolyhedralFan readFan(string const &filename, bool onlyMaximal=true, IntegerVector *w=0, set<int> const *conesIndice=0, SymmetryGroup const *sym=0, bool readCompressedIfNotSym=false);
  int getAmbientDimension()const;
  int getMaxDimension()const;
  int getMinDimension()const;
  friend PolyhedralFan refinement(const PolyhedralFan &a, const PolyhedralFan &b, int cutOffDimension=-1, bool allowASingleConeOfCutOffDimension=false);
  friend PolyhedralFan product(const PolyhedralFan &a, const PolyhedralFan &b);
  IntegerVectorList getRays(int dim=1);//This can be called for other dimensions than 1. The term "Rays" still makes sense modulo the common linearity space
  IntegerVectorList getRelativeInteriorPoints();
  PolyhedralCone const& highestDimensionalCone()const; //returns one of the cones in the fan with highest dimension
  void insert(PolyhedralCone const &c);
  void insertFacetsOfCone(PolyhedralCone const &c);
  void remove(PolyhedralCone const &c);
  void removeAllExcept(int a); //deletes all except the first a cones
  void removeAllLowerDimensional();
  /**
     Since the cones stored in a PolyhedralFan are cones of a
     polyhedral fan, it is possible to identify non maximal cones by
     just checking containment of relative interior points in other
     cones. This routine removes all non-maximal cones.
   */
  void removeNonMaximal();
  /**
     Returns the number of cones stored in the fan. This is not the number of cones in the fan in a mathematical sense.
   */
  int size()const;
  int dimensionOfLinealitySpace()const;
  void makePure();
  bool contains(PolyhedralCone const &c)const;
  /**
   * For a vector contained in the support of the fan represented by the fan object, this function
   * computes the cone that contains the vector in its relative interior.
   */
  PolyhedralCone coneContaining(IntegerVector const &v)const;
  PolyhedralFan facetComplex()const;
  /**
     This routine computes all cones in the fan and returns them as a new PolyhedralFan. This routine is extremely slow.
   */
  PolyhedralFan fullComplex()const;
  PolyhedralFan facetComplexSymmetry(SymmetryGroup const &sym, bool keepRays=false, bool dropLinealitySpace=false)const;
  PolyhedralFan rayComplexSymmetry(SymmetryGroup const &sym)const;
  IntegerVectorList getRaysInPrintingOrder(SymmetryGroup const *sym, bool upToSymmetry=false)const;
  IncidenceList getIncidenceList(SymmetryGroup *sym=0)const;
  bool isEmpty()const;
  bool isRefinementOf(PolyhedralFan const &f)const;
  /**
     Computes the link of the face containing w in its relative interior.
   */
  PolyhedralFan link(IntegerVector const &w)const;
  PolyhedralFan link(IntegerVector const &w, SymmetryGroup *sym)const;
  /**
     If the cones have been assigned linear forms then this specifies
     a piecewise linear function on the support of the fan. This
     routine evaluates that linear function.
   */
  int64 evaluatePiecewiseLinearFunction(IntegerVector const &x)const;

  typedef PolyhedralConeList::const_iterator coneIterator;
  PolyhedralFan::coneIterator conesBegin()const;
  PolyhedralFan::coneIterator conesEnd()const;
  /**
     Computes the volume of the d-dimensional cones stored in the
     fan. If the symmetry group is specified each volume of a cone is
     multiplied by the size of its orbit.
   */
  FieldElement volume(int d, SymmetryGroup *sym=0)const;
  /**
     Assuming that the fan is pure, this function tests whether the
     fan is connected in codimension 1.  IF SYMMETRY IS SPECIFIED THEN
     THE CHECK IS NOT COMPLETE AND IT WILL ONLY BE CHECKED IF THAT GIVEN
     TWO CONES THAT EXIST ELEMENTS IN THE TWO RESPECTIVE ORBITS
     WHICH ARE CONNECTED BY A RIDGE PATH.
   */
  bool isConnected(SymmetryGroup *sym=0)const;
  PolyhedralFan negated()const;


  /**
     Checks if c can be inserted in the fan without violating the fan
     axioms. Very slow.
   */
  bool isCompatible(PolyhedralCone const &c)const;
  /**
     Change the polyhedral structure of the fan and add in the cone
     c. Overlap of c and the cones are allowed. This routine is
     extremely slow.
   */
  void merge(PolyhedralCone const &c);
  /**
   * Converts a PolyhedralFan into a SymmetricComplex. This is used for homology computations, but not for printing yet.
   */
  class SymmetricComplex toSymmetricComplex(SymmetryGroup *sym);
};


void addFacesToSymmetricComplex(class SymmetricComplex &c, PolyhedralCone const &cone, IntegerVectorList const &facetCandidates, IntegerVectorList const &generatorsOfLinealitySpace);
void addFacesToSymmetricComplex(SymmetricComplex &c, set<int> const &indices, IntegerVectorList const &facetCandidates, int dimension, int multiplicity);


#endif
