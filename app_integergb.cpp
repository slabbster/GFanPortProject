/*
 * app_integergb.cpp
 *
 *  Created on: Dec 14, 2010
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

#include "integergb.h"

#include <ostream>

using namespace std;

class IntegerGBApplication : public GFanApplication
{
  SimpleOption gbOption;
  SimpleOption initialOption;
  SimpleOption gFanOption;
  SimpleOption gConeOption;
  SimpleOption listOption;
  SimpleOption optionInputIsGroebnerBasis;
public:
  const char *helpText()
  {
    return "This program is an experimental implementation of Groebner bases for ideals in Z[x_1,...,x_n].\n"
    "Several operations are supported by specifying the appropriate option:\n"
    " (1) computation of the reduced Groebner basis with respect to a given vector (tiebroken lexicographically),\n"
    " (2) computation of an initial ideal,\n"
    " (3) computation of the Groebner fan,\n"
    " (4) computation of a single Groebner cone.\n"
    "Since Gfan only knows polynomial rings with coefficients being elements of a field, the ideal is specified by giving a set of polynomials in the polynomial ring Q[x_1,...,x_n]. That is, by using Q instead of Z when specifying the ring. The ideal MUST BE HOMOGENEOUS (in a positive grading) for computation of the Groebner fan. Non-homogeneous ideals are allowed for the other computations if the specified weight vectors are positive.\n"
    "NOTE: This program is experimental and expected to change behaviour in future releases, so don't write your SAGE and M2 interfaces just yet.\n";
  }
  IntegerGBApplication():
    gbOption("--groebnerBasis","Asks the program to compute a marked Groebner basis with respect to a weight vector tie-broken lexicographically.\n"
        "The input order is: Ring ideal vector.\n"),
    initialOption("--initialIdeal","Asks the program to compute an initial ideal with respect to a vector. "
        "The input order is: Ring ideal vector.\n"),
    gFanOption("--groebnerFan","Asks the program to compute the Groebner fan. \n "
        "The input order is: Ring ideal.\n"),
    gConeOption("--groebnerCone","Asks the program to compute a single Groebner cone containing the specified vector in its relative interior. The output is stored as a fan. "
        "The input order is: Ring ideal vector."),
    listOption("-m","For the operations taking a vector as input, read in a list of vectors instead, and perform the operation for each vector in the list."),
    optionInputIsGroebnerBasis("-g",
                               "Tells the program that the input is already a Groebner basis (with the initial term of each polynomial being "
                               "the first ones listed). Use this option if the usual --groebnerFan is too slow.\n")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_overintegers";
  }
  int main()
  {
    if(gbOption.getValue()+initialOption.getValue()+gFanOption.getValue()+gConeOption.getValue()!=1)
      {
        debug<<"WRONG COMBINATION OF COMMAND LINE OPTIONS\n";
        assert(0);
      }
    LexicographicTermOrder tieBreaker;
    FileParser P(Stdin);
    PolynomialSet a=P.parsePolynomialSetWithRing();
    int n=a.getRing().getNumberOfVariables();
//    IntegerVectorList omegas=P.parseIntegerVectorList();

    if(gFanOption.getValue())
      {
        SymmetryGroup G(n);
        SymmetricTargetFanBuilder target(n,G);

        if(optionInputIsGroebnerBasis.getValue())
          {
            IntegerVectorList empty;
            IntegerVector w=PolyhedralCone(wallInequalities(a),empty).getRelativeInteriorPoint();
            WeightReverseLexicographicTermOrder tieBreaker(w);
            zAutoReduce(&a,tieBreaker);
          }
        else
          {
            WeightReverseLexicographicTermOrder tieBreaker(IntegerVector::allOnes(n));
            zBuchberger(a,tieBreaker);
          }
        IntegerGroebnerFanTraverser traverser(a);
        symmetricTraverse(traverser,target);

        AsciiPrinter Q(Stdout);
        target.getFanRef().printWithIndices(&Q,
            FPF_default);
      }
    else
      {
        IntegerVectorList omegas;
        if(listOption.getValue())
          omegas=P.parseIntegerVectorList();
        else
          omegas.push_back(P.parseIntegerVector());

        for(IntegerVectorList::const_iterator i=omegas.begin();i!=omegas.end();i++)
          {
            if(i->size()!=a.getRing().getNumberOfVariables())
              {
                debug<<"ERROR: The number of entries of the weight vector is not equal to the number of variables in the ring.\n";
                assert(0);
              }
            if(gbOption.getValue())
              {
                WeightTermOrder T(*i);
                zBuchberger(a,T);
                pout<<a.getRing()<<a;
              }
            else if(initialOption.getValue())
              {
                WeightTermOrder T(*i);
                zBuchberger(a,T);

                pout<<a.getRing();
                pout<<initialForms(a,*i);
              }
            else if(gConeOption.getValue())
              {
                WeightTermOrder T(*i);
                zBuchberger(a,T);
                IntegerVectorList inequalities=wallInequalities(a);
                inequalities.push_back(IntegerVector::standardVector(n,0));
                IntegerVectorList empty;

                PolyhedralCone C(inequalities,empty,n);
                C=C.faceContaining(*i);
                C.canonicalize();

                PolyhedralFan F(C.ambientDimension());
                F.insert(C);
                F.printWithIndices(&pout,
                    FPF_default);
              }
      }
      }
/*
    for(IntegerVectorList::const_iterator i=omegas.begin();i!=omegas.end();i++)
      {
        debug<<"INTEGRAL GROEBNER BASIS:\n";
        WeightTermOrder T(*i);
        zBuchberger(a,T);
        debug<<"INTEGRAL GROEBNER BASIS:\n";
        debug<<a;

//        PolynomialRing ZModPZRing=residuePolynomialRing(a.getRing(), prime);

        debug<<"INTEGRAL INITIAL IDEAL:\n";
//        debug<<ZModPZRing;
        debug<<initialForms(a,*i);

//        debug<<"P-ADIC INITIAL TERMS:\n";
 //       debug<<ZModPZRing;
//        debug<<pAdicInitialTerms(ZModPZRing,a,prime,omega,tieBreaker);

//        debug<<"AUTOREDUCED:\n";
//        pAdicAutoReduce(a,prime,omega,tieBreaker);
//        debug<<a;

    debug<<"MAXIMAL GROEBNER CONE:\n";
    IntegerVectorList inequalities=wallInequalities(a);
    IntegerVectorList empty;
    PolyhedralCone C(inequalities,empty,n);
    debug<<C;
    debug<<"RAYS:\n";
    debug<<C.extremeRays();
      }

    {
      SymmetryGroup G(n);
      SymmetricTargetFanBuilder target(n,G);

      IntegerGroebnerFanTraverser traverser(a);
      symmetricTraverse(traverser,target);

      AsciiPrinter Q(Stdout);
      target.getFanRef().printWithIndices(&Q,
                                  FPF_default);

    }
*/    return 0;
  }
};

static IntegerGBApplication theApplication;
