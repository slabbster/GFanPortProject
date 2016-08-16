/*
 * app_fancoarsening.cpp
 *
 *  Created on: Jul 19, 2010
 *      Author: anders
 */

#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "lp.h"
#include "gfanapplication.h"
#include "polyhedralcone.h"
#include "polyhedralfan.h"
#include "symmetry.h"

#include "polymakefile.h"

#include "traverser_sphere.h"

class FanCoarseningApplication : public GFanApplication
{

  StringOption inputOption;
//  SimpleOption resultantOption;

  multimap<IntegerVector,int> coneIndicesAtRidge;
  vector<PolyhedralCone> F1,F2;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }

  IntegerVectorList inequalitiesForCone(IntegerVector const &codim1interior, IntegerVector const &newDirection)
  {

    int i=0;

  }

  const char *helpText()
  {
    return "This program reads in a codimension-one fan (with support equal to the codimension-1 skeleton of the normal fan of a polytope) and constructs the normal fan of the polytope.\n";
  }
  FanCoarseningApplication():
    inputOption("-i","Specify the name of the input file.","polymake.out")
 //   resultantOption("--resultant","Take codim 1 skeleton and wipe out bad cones.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_fancoarsening";
  }

/*  bool zeroOrTwo(int v)
  {
    return (v==0) || (v==2);
  }*/

  class MyTarget : public SymmetricTarget
  {
    list<IntegerVector> normals;
    bool process(ConeTraverser &traverser)
    {
      SphereTraverser *t=dynamic_cast<SphereTraverser*> (&traverser);
      normals.push_back(t->currentNormal);
      return true;
    }
  public:
    PolyhedralCone toCone()
    {
      PolyhedralCone c(normals,IntegerVectorList(),normals.begin()->size());
      c.canonicalize();
      return c;
    }
  };


  class MyTraverser: public ConeTraverser
  {
    int currentConeIndex;
    vector<PolyhedralCone> const &cones;
    map<IntegerVector,list<int> > adjacency;
    PolyhedralCone theCone;
    //void updatePolyhedralCone();
    void moveToCone(IntegerVector const &v, IntegerVector const &normal)
    {
      SphereTraverser traverser(cones,adjacency,v,normal);

      MyTarget myTarget;
      SymmetryGroup s(normal.size());
      symmetricTraverse(traverser,myTarget);
      theCone=myTarget.toCone();
    }
  public:
//    IntegerVector currentNormal;
          MyTraverser(vector<PolyhedralCone> const &cones_, map<IntegerVector,list<int> > const &adjacency_, IntegerVector const &v, IntegerVector const &normal):
            cones(cones_),
            adjacency(adjacency_)
            {
              moveToCone(v,normal);
          }
          virtual void changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)
          {
            moveToCone(ridgeVector,rayVector);
          }
          virtual IntegerVectorList link(IntegerVector const &ridgeVector)
          {
            IntegerVector v=theCone.link(ridgeVector).getUniquePoint();
            IntegerVectorList ret;
            ret.push_back(v);
            ret.push_back(-v);
            return ret;
          }
          PolyhedralCone & refToPolyhedralCone()
          {
            return theCone;
          }
  };




  int main()
  {
    PolyhedralFan f1=PolyhedralFan::readFan(inputOption.getValue(),true,0,0,0);
    AsciiPrinter P(Stdout);

    for(PolyhedralFan::coneIterator i=f1.conesBegin();i!=f1.conesEnd();i++)F1.push_back(*i);
    map<IntegerVector,list<int> > adjacency;
    for(int i=0;i<F1.size();i++)
      {
        PolyhedralFan temp(F1[i].ambientDimension());
        temp.insert(F1[i]);
        PolyhedralFan facets=temp.facetComplex();
        for(PolyhedralFan::coneIterator j=facets.conesBegin();j!=facets.conesEnd();j++)
          adjacency[j->getUniquePoint()].push_back(i);
      }
//    SphereTraverser(vector<PolyhedralCone> const &cones_, map<IntegerVector,list<int> > const &ajacency_, IntegerVector const &startCone, IntegerVector const &normal);

    SphereTraverser traverser(F1,adjacency,F1[0].getUniquePoint(),*F1[0].getEquations().begin());

    MyTarget myTarget;
    SymmetryGroup s(f1.getAmbientDimension());
//    SymmetricTargetFanBuilder myTarget(f1.getAmbientDimension(),s);
    symmetricTraverse(traverser,myTarget);

//    myTarget.getFanRef().printWithIndices(&pout);
    myTarget.toCone().print(&debug);



    {
      MyTraverser traverser(F1,adjacency,F1[0].getUniquePoint(),*F1[0].getEquations().begin());

      SymmetryGroup s(f1.getAmbientDimension());
      SymmetricTargetFanBuilder myTarget(f1.getAmbientDimension(),s);
      symmetricTraverse(traverser,myTarget);
      myTarget.getFanRef().printWithIndices(&pout);
    }

    return 0;
  }
};

static FanCoarseningApplication theApplication;

