/*
 * app_traversetropicalintersection.cpp
 *
 *  Created on: Jan 4, 2011
 *      Author: anders
 */

#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "division.h"
#include "log.h"
#include "polyhedralcone.h"
#include "tropical2.h"
#include "wallideal.h"
#include "halfopencone.h"

#include "integergb.h"

#include <ostream>

using namespace std;





class TropicalIntersectionTraverser: public ConeTraverser
{
  PolynomialSet allPolys;
  PolynomialSet iForms;
  PolyhedralCone theCone;
  int n;//,d;
  int prime;
  bool collect;
  IntegerVectorList collection;
  void updatePolyhedralCone()
  {
    theCone=PolyhedralCone(fastNormals(wallInequalities(allPolys)),wallInequalities(iForms),n);
    theCone.canonicalize();
    if(collect)
      {
        collection.push_back(theCone.getRelativeInteriorPoint());
        int temp=collection.size();
        if((temp&(2*(temp-1)))==temp)pout<<"RELINT:\n"<<collection;
      }
  }
public:
  TropicalIntersectionTraverser(PolynomialSet const &generators, IntegerVector const &omega):
    allPolys(generators),
    iForms(generators.getRing()),
    n(generators.getRing().getNumberOfVariables()),
    collect(true)
  {
    IntegerVectorList l;
    l.push_back(omega);
    MatrixTermOrder T(l);
    iForms=initialForms(generators,omega);
    allPolys.mark(T);
    iForms.mark(T);
    updatePolyhedralCone();
  }
  virtual void changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)
  {
    IntegerVectorList l;
    l.push_back(ridgeVector);l.push_back(rayVector);
    MatrixTermOrder T(l);
    allPolys.mark(T);
    iForms=initialForms(initialForms(allPolys,ridgeVector),rayVector);
    iForms.mark(T);
    updatePolyhedralCone();
  }
  virtual IntegerVectorList link(IntegerVector const &ridgeVector)
  {
    PolyhedralFan l=tropicalHyperSurfaceIntersectionClosed(n, initialForms(allPolys,ridgeVector));
    assert(l.getMaxDimension()==l.dimensionOfLinealitySpace()+1);

    IntegerVectorList ret;
    for(set<PolyhedralCone>::const_iterator i=l.conesBegin();i!=l.conesEnd();i++)
      if(!(i->getUniquePoint().isZero()))ret.push_back(i->getUniquePoint());

    assert(ret.size()>1);

    return ret;
  }
  PolyhedralCone & refToPolyhedralCone()
  {
    return theCone;
  }
};





class TraverseTropicalIntersectionApplication : public GFanApplication
{
  SimpleOption optionSymmetry;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program \n";
  }
  TraverseTropicalIntersectionApplication():
    optionSymmetry("--symmetry","Do computations up to symmetry and group the output accordingly. If this option is used the program will read in a list of generators for a symmetry group after the other input.")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_traversetropicalintersection";
  }
  int main()
  {
    FileParser P(Stdin);
    PolynomialSet a=P.parsePolynomialSetWithRing();
    int n=a.getRing().getNumberOfVariables();
    IntegerVector omega=P.parseIntegerVector();


    {
      SymmetryGroup G(n);
      if(optionSymmetry.getValue())
        {
          IntegerVectorList generators=P.parseIntegerVectorList();
          G.computeClosure(generators);
          G.createTrie();
        }

      SymmetricTargetFanBuilder target(n,G);

      TropicalIntersectionTraverser traverser(a,omega);
      symmetricTraverse(traverser,target,&G);

      AsciiPrinter Q(Stdout);
      target.getFanRef().printWithIndices(&Q,
                                  FPF_default);

    }
    return 0;
  }
};

static TraverseTropicalIntersectionApplication theApplication;
