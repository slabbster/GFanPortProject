#ifndef TRIANGULATION_H_INCLUDED
#define TRIANGULATION_H_INCLUDED

#include "matrix.h"
#include <list>

using namespace std;

class Triangulation{
 public:
  /**
     This class represents an oriented simplex. The simplex is stored
     as the list of indices of the vertices it contains together with
     an orientation flag which is either +1 or -1. The ordering of the
     indices is important. The list of indices can be sorted by
     calling coneSort(). This will change the orientation flag
     appropriately.
   */
  class Cone : public list<int>
  {
  public:
    int changeSign;
  Cone():
    changeSign(1)
      {
      }
  };
  /**
     The class Cone2 is just a short name for a list of integers.
   */
  typedef list<int> Cone2;
  // private:
  static int coneDim(Cone2 const &c, IntegerMatrix const &rays);
  static IntegerVectorList coneToVectorList(Cone2 const &c, IntegerMatrix const &rays);
  static void coneSort(Cone &c);
  static Cone firstSimplex(Cone const &c, IntegerMatrix const &rays);
  static IntegerVectorList coneComplement(Cone c, IntegerMatrix const &rays);//returns generators of orth. complement.
  static int signVisible(int v, Cone const &c, IntegerVectorList const &complement, IntegerMatrix const &rays);
  static list<Cone> triangulateRek(int d, Cone2 const &c, IntegerMatrix const &rays, bool revlex=false, bool ignoreContainedRays=false);

 public:
  /**
     Computes a triangulation of the rows of the matrix rays indexed
     by c. The tringulation is lexicographic unless revlex is set in
     which case a reverse lexicographic triangulation is
     computed. Notice that a reverselex triangulation does not always
     exist.
     CHOOSING REVLEX DOES NOT WORK YET.
   */
  static list<Cone> triangulate(Cone2 c, IntegerMatrix const &rays, bool revlex=false); //computes a lexicographic triangulation
  static list<Cone> triangulate(IntegerMatrix const &rays, bool revlex=false); //computes a lexicographic triangulation
  static list<Cone> boundary(list<Cone> cones);
  static IntegerVectorList normals(IntegerMatrix &rays);
  /**
     Takes a Cone list and turns it into a vector of lists of integers, thereby remove the orientation flag stored in the Cones.
   */
  static vector<list<int> > removeOrientation(list<Cone> const &triangulation);
  /**
     Takes a Cone vector and turns it into a vector of lists of integers, thereby remove the orientation flag stored in the Cones.
   */
  static vector<list<int> > removeOrientation(vector<Cone> const &triangulation);
};

#endif
