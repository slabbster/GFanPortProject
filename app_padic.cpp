/*
 * app_padic.cpp
 *
 *  Created on: Dec 1, 2010
 *      Author: anders
 */

#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "division.h"
#include "log.h"
#include "polyhedralcone.h"

#include "padic.h"

#include <ostream>

using namespace std;

class PAdicApplication : public GFanApplication
{
  IntegerOption primeOption;
  SimpleOption gbOption;
  SimpleOption initialOption;
  SimpleOption gComplexOption;
  SimpleOption gPolyhedronOption;
  SimpleOption listOption;
public:
  const char *helpText()
  {
    return "This program is an experimental implementation of p-adic Groebner bases as proposed by Diane Maclagan.\n"
    "Several operations are supported by specifying the appropriate option:\n"
    " (1) computation of Groebner basis with respect to a given vector (tiebroken lexicographically),\n"
    " (2) computation of the p-adic initial ideal,\n"
    " (3) computation of the p-adic Groebner complex as defined by Maclagan and Sturmfels,\n"
    " (4) computation of a single polyhedron of the p-adic Groebner complex.\n"
    "The input ideal should be an ideal of the polynomial ring with coefficient field Q. The valuation is specified with the option -p. The ideal MUST BE HOMOGENEOUS (in a positive grading).\n"
    "Since gfan can only handle fans and not polyhedral complexes in general, what is computed as the Groebner complex is actually the \"fan over\" the complex - in other words, the first coordinate is supposed to be 1 in the output fan.\n"
    "Similarly, the weight vectors must be specified in an homogeneous way, for example by adding an additional 1 entry as first coordinate. (If fractions are needed, use the entry as a common denominator.) "
    "NOTE: This program is experimental and expected to change behaviour in future releases, so don't write your SAGE and M2 interfaces just yet. In particular this program uses the tropical minimum-convention!!\n";
  }
  PAdicApplication():
    primeOption("-p","Defines the prime used for the valuation.",2),
    gbOption("--groebnerBasis","Asks the program to compute a marked Groebner basis with respect to a weight vector (tie-broken lexicographically).\n"
        "The input order is: Ring ideal vector.\n"),
    initialOption("--initialIdeal","Asks the program to compute an initial ideal with respect to a vector. "
        "The input order is: Ring ideal vector.\n"),
    gComplexOption("--groebnerComplex","Asks the program to compute the p-adic Groebner complex. \n "
        "The input order is: Ring ideal.\n"),
    gPolyhedronOption("--groebnerPolyhedron","Asks the program to compute a single polyhedron of the Groebner complex containing the specified vector in its relative interior. The output is stored as a fan. "
        "The input order is: Ring ideal vector."),
    listOption("-m","For the operations taking a vector as input, read in a list of vectors instead, and perform the operation for each vector in the list.")
        {
    registerOptions();
  }
  const char *name()
  {
    return "_padic";
  }
  int main()
  {
    if(gbOption.getValue()+initialOption.getValue()+gComplexOption.getValue()+gPolyhedronOption.getValue()!=1)
      {
        debug<<"WRONG COMBINATION OF COMMAND LINE OPTIONS\n";
        assert(0);
      }
    LexicographicTermOrder tieBreaker;
    FileParser P(Stdin);
    int prime=primeOption.getValue();
    PolynomialSet a=P.parsePolynomialSetWithRing();
    int n=a.getRing().getNumberOfVariables();

    if(gComplexOption.getValue())
      {
        SymmetryGroup G(n+1);
        SymmetricTargetFanBuilder target(n+1,G);

        PAdicGroebnerFanTraverser traverser(a,prime);
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
            if(i->size()!=a.getRing().getNumberOfVariables()+1)
              {
                debug<<"ERROR: The number of entries of the weight vector is not one higher than the number of variables in the ring.\n";
                assert(0);
              }
            if(gbOption.getValue())
              {
                //debug<<"P-ADIC GROEBNER BASIS:\n";
                pAdicBuchberger(a,prime,*i,tieBreaker);
                pout<<a.getRing()<<a;
              }
            else if(initialOption.getValue())
              {
                pAdicBuchberger(a,prime,*i,tieBreaker);
                PolynomialRing ZModPZRing=residuePolynomialRing(a.getRing(), prime);

//                debug<<"P-ADIC INITIAL IDEAL:\n";
                pout<<ZModPZRing;
                pout<<pAdicInitialForms(ZModPZRing,a,prime,*i);
              }
//        PolynomialRing ZModPZRing=residuePolynomialRing(a.getRing(), prime);

//        debug<<"P-ADIC INITIAL IDEAL:\n";
//        debug<<ZModPZRing;
//        debug<<pAdicInitialForms(ZModPZRing,a,prime,*i);

//        debug<<"P-ADIC INITIAL TERMS:\n";
 //       debug<<ZModPZRing;
//        debug<<pAdicInitialTerms(ZModPZRing,a,prime,omega,tieBreaker);

//        debug<<"AUTOREDUCED:\n";
//        pAdicAutoReduce(a,prime,omega,tieBreaker);
//        debug<<a;
            else if(gPolyhedronOption.getValue())
              {
                pAdicBuchberger(a,prime,*i,tieBreaker);
                IntegerVectorList inequalities=normalPolyhedralInequalities(a,prime,*i,tieBreaker);
                inequalities.push_back(IntegerVector::standardVector(n+1,0));
                IntegerVectorList empty;

                PolyhedralCone C(inequalities,empty,n+1);
                C=C.faceContaining(*i);
                C.canonicalize();

                PolyhedralFan F(C.ambientDimension());
                F.insert(C);
                F.printWithIndices(&pout,
                    FPF_default);
//                pout<<C;
              }
/*    debug<<"GROEBNER POLYHEDRON:\n";
    IntegerVectorList inequalities=normalPolyhedralInequalities(a,prime,omega,tieBreaker);
    inequalities.push_back(IntegerVector::standardVector(n+1,0));
    IntegerVectorList empty;

    PolyhedralCone C(inequalities,empty,n+1);

    debug<<C;
    debug<<"RAYS:\n";
    debug<<C.extremeRays();
*/
      }
      }

    return 0;
  }
};

static PAdicApplication theApplication;
