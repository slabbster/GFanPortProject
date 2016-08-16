#ifndef SYMMETRICCOMPLEX_H_INCLUDED
#define SYMMETRICCOMPLEX_H_INCLUDED

#include <set>
#include <string>
#include <map>

using namespace std;

#include "symmetry.h"
#include "matrix.h"

class SymmetricComplex{
  int n;
  IntegerMatrix vertices;
  map<IntegerVector,int> indexMap;
  SymmetryGroup sym;
  IntegerVector dimensionsAtInfinity()const;
 public:
  class Cone
  {
    bool isKnownToBeNonMaximalFlag;
    //    bool ignoreSymmetry; //useful when computing faces and extracting symmetries
  public:
    vector<int> indices;//always sorted
    Cone(set<int> const &indices_, int dimension_, int multiplicity_, bool sortWithSymmetry, SymmetricComplex const &complex);
    set<int> indexSet()const;
    int dimension;
    int multiplicity;
    bool isKnownToBeNonMaximal()const{return isKnownToBeNonMaximalFlag;}
    void setKnownToBeNonMaximal(){isKnownToBeNonMaximalFlag=true;}
    bool isSubsetOf(Cone const &c)const;
    SymmetricComplex::Cone permuted(IntegerVector const &permutation, SymmetricComplex const &complex, bool withSymmetry)const;
    //    bool operator<(Cone const & b)const;
    /*    IntegerVector relativeInteriorPoint;
    IntegerVector smallestRepresentative;
    IntegerVector summary;*/
    IntegerVector sortKey;
    IntegerVector sortKeyPermutation;
    //    void computeRelativeInteriorPoint(SymmetricComplex const &complex);
    //    void computeSmallestRepresentative(SymmetricComplex const &complex);
    bool operator<(const Cone & b)const;
    //    void setIgnoreSymmetry(bool b){ignoreSymmetry=b;}
    bool isSimplicial(int linealityDim)const;
    void remap(SymmetricComplex &complex);
/**
 * This routine computes a basis for the orthogonal complement of the cone.
 * Notice that the lineality space, which is unknown at the time, is ignored.
 * This routine is deterministic and used for orienting the faces when computing homology.
 */
    IntegerVectorList orthogonalComplement(SymmetricComplex &complex)const;
  };
  typedef set<Cone> ConeContainer;
  ConeContainer cones;
  int dimension;
  SymmetricComplex(int n_, IntegerVectorList const &v, SymmetryGroup const &sym_);
  /**
   * Returns a reference to the matrix of vertices on which the complex is build.
   * The reference is valid as the Symmetric complex object exists.
   */
  IntegerMatrix const &getVertices()const{return vertices;}
  bool contains(Cone const &c)const;
  void insert(Cone const &c);
  int getMaxDim()const;
  int getMinDim()const;
  bool isMaximal(Cone const &c)const;
  bool isPure()const;
  IntegerVector fvector(bool boundedPart=false)const;
  string toString(int dimLow, int dimHigh, bool onlyMaximal, bool group, ostream *multiplicities=0, bool compressed=false, bool tPlaneSort=false, bool xml=false)const;
  bool isSimplicial()const;
  /**
     Calling this function will change the representative of each cone
     orbit by "applying" the permutation which will give the sortkey to
     the set of indices of the cone.
   */
  void remap();
  /**
   * Looks up the index of the vector among the vertices.
   */
  int indexOfVertex(IntegerVector const &v)const;
  int numberOfConesOfDimension(int d)const;
  /**
   * Given a cone this returns its index among all cones of that dimension.
   * Used for assigning "names" to cones.
   */
  int dimensionIndex(Cone const &c);
  /**
   * This routine is used for constructing the boundary map for homology computations.
   */
  void boundary(Cone const &c, vector<int> &indices, vector<int> &signs);
/**
 * This routine computes the ith boundary map for homology as a matrix.
 */
  IntegerMatrix boundaryMap(int i);
};

#endif
