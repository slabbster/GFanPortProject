#ifndef __halfopencone_h
#define __halfopencone_h

#include "polyhedralcone.h"
#include "termorder.h"

#include "polyhedralfan.h"

class HalfOpenCone{
  static void appendList(IntegerVectorList &to, IntegerVectorList const &from, int appendValue);
  int liftedDimension;//ambient
 public: PolyhedralCone lifted;//remove public
  static IntegerVectorList shrink(const IntegerVectorList &l);
  HalfOpenCone(int dimension_, PolyhedralCone const &lifted_);
 public:
  int dimension;//ambient
  HalfOpenCone(PolyhedralCone C, TermOrder const &t);//only for full dimensional cones!
  HalfOpenCone(int dimension_, IntegerVectorList const &equations, IntegerVectorList const &nonstrict, IntegerVectorList const &strict, bool findFacets=false, bool canonicalize=false);
  HalfOpenCone(int ambientDimension);//full space
  bool isEmpty();
  friend HalfOpenCone intersection(const HalfOpenCone &a, const HalfOpenCone &b, bool findFacets=false);
  friend bool haveEmptyIntersection(const HalfOpenCone &a, const HalfOpenCone &b);
  PolyhedralCone closure();
  void splitIntoRelativelyOpenCones(list<HalfOpenCone> &l);
  void print(class Printer &p)const;
  bool contains(IntegerVector const &v)const;
  friend bool operator<(HalfOpenCone const &a, HalfOpenCone const &b);
  void canonicalize(){lifted.canonicalize();}

  /**
     Remove all coordinates from the space except those listed in chosen.
   */
  HalfOpenCone withChosenCoordinates(list<int> chosen)const;
  HalfOpenCone rewrite(FieldMatrix const &A, list<int> nonPivots)const;
  HalfOpenCone rewriteExpand(list<int> pivots, IntegerVectorList const &newEquations)const;
};

typedef list<HalfOpenCone> HalfOpenConeList;

HalfOpenConeList orientedBoundary(PolyhedralCone C, TermOrder const &t, HalfOpenCone *restrictingCone=0);

HalfOpenConeList splitIntoRelativelyOpenCones(HalfOpenConeList const &l);

class HalfOpenConeProcessor
{
public:
  bool savePartialResult;
  HalfOpenConeProcessor():savePartialResult(false)
  {
  }
  virtual void process(HalfOpenCone const &c)=0;
  void setSave()
  {
    savePartialResult=true;
  }
};

void tropicalHyperSurfaceIntersectionWithProcessor(int dimension, PolynomialSet const &g, HalfOpenConeProcessor &myProcessor, PolyhedralCone *restrictingCone=0, bool expand=false, int intervalLow=-1, int intervalHigh=-1);

HalfOpenConeList tropicalHyperSurfaceIntersection(int dimension, PolynomialSet const &g, HalfOpenCone *restrictingCone=0);
void tropicalHyperSurfaceIntersectionInSubspace(int dimension, PolynomialSet const &G, HalfOpenCone *restrictingCone, HalfOpenConeProcessor &processor, int intervalLow=-1, int intervalHigh=-1);
//HalfOpenConeList tropicalHyperSurfaceIntersectionInSubspace(int dimension, PolynomialSet const &g, HalfOpenCone *restrictingCone=0);
PolyhedralFan tropicalHyperSurfaceIntersectionClosed(int dimension, PolynomialSet const &g, PolyhedralCone *restrictingCone=0, bool expand=false, bool saveResult=false, int intervalLow=-1, int intervalHigh=-1);

void printHalfOpenConeList(HalfOpenConeList const &l, class Printer & p);//Assuming that the union is closed and the cones have a common linearity space


PolyhedralFan faceComplexOfCone(HalfOpenCone &c);

/**
 * This routine tricks the tropical hypersurface intersection
 * cone into checking whether the tropical stable intersection
 * of the list of polynomials is non-empty. This can be used
 * to check if a weight vector is contined in the stable
 * intersection - simply call this routine on the initial
 * forms of the generators.
 */
bool nonEmptyStableIntersection(PolynomialSet const &g);

#endif
