#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "nbody.h"
#include "field_rationals.h"
#include <iostream>

class NBodyApplication : public GFanApplication
{
  IntegerOption NOption;
  SimpleOption optionWithMasses;
  SimpleOption optionSymmetric;
  SimpleOption optionDziobek;
  IntegerOption optionDeterminants;
  SimpleOption optionMM;
  SimpleOption optionSVariables;
  SimpleOption optionLaurent;
  //SimpleOption optionKaloshin;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program computes the AC equations for the nbody problem.\n";
  }
  NBodyApplication():
    NOption("-N","Specify number of bodies.",3),
    optionWithMasses("--masses","Include mass variables."),
    optionSymmetric("--symmetric","Produce the symmetric equations"),
    optionDziobek("--dziobek","Produce the Dziobek equations"),
    optionDeterminants("--codim","Produce the determinant equations",-1),
    optionMM("--mm","Add in additional mass equation"),
    optionSVariables("-s","Treat the Sij polynomials as variables"),
    optionLaurent("--laurent","Output Laurent polynomials instead of multiplying")
//    optionKaloshin("--kaloshin","Produce the polynomials from the Albouy Kaloshin paper")
  {
    registerOptions();
  }
  const char *name()
  {
    return "_nbody";
  }
/*
  int offset(int N, int i, int j)
  {
  }
  Polynomial varz(PolynomialRing const &r, int i)
  {
    int N=4;
    return r.ithVariable(i);
  }
  Polynomial varw(PolynomialRing const &r, int i)
  {
    int N=4;
    return r.ithVariable(N+i);
  }
  Polynomial varZ(PolynomialRing const &r, int i, int j)
  {
    int N=4;
    int sign=(i<j)?1:-1;
    if(j<i){i-=j;j-=i;i-=j;}
    return r.polynomialFromField(Q.zHomomorphism(sign)).ithVariable(offset(N,i,j));
  }
*/
  int main()
  {
    AsciiPrinter PP(stdout);

 /*   if(optionKaloshin.getValue())
      {
        PolynomialRing r=StringParser("Q[z1,z2,z3,z4,w1,w2,w3,w4,Z12,Z13,Z14,Z23,Z24,Z34,W12,W13,W14,W23,W24,W34]");

        PolynomialSet g(r);
        int N=4;
        for(int i=0;i<N;i++)
          for(int j=0;j<i;j++)
            {
              g.push_back((varz(r,i)-varz(r,j))*(varw(r,i)-varw(r,j))*(varw(r,i)-varw(r,j))*(varw(r,i)-varw(r,j))*varW(r,i,j)*varW(r,i,j)-r.one());
              g.push_back((varw(r,i)-varw(r,j))*(varz(r,i)-varz(r,j))*(varz(r,i)-varz(r,j))*(varz(r,i)-varz(r,j))*varZ(r,i,j)*varZ(r,i,j)-r.one());
            }
        for(int i=0;i<N;i++)
          {
            Polynomial A=-varz(r,i);
            Polynomial B=-varw(r,i);

            for(int j=0;j<N;j++)
            if(i!=j)
              {
                A+=varZ(r,i,j);//times m
                B+=varW(r,i,j);//times m
              }
            g.push_back(A);
            g.push_back(B);
          }
        PP<<g.getRing();
        PP<<g;
        return 0;
      }
   */

    PolynomialSet g=AlbouyChencinerEquations(NOption.getValue(),optionWithMasses.getValue(),optionSymmetric.getValue(),optionSVariables.getValue(),!optionLaurent.getValue());

    if(optionSVariables.getValue())
      {
	PolynomialSet s=SEquations(g.getRing(),NOption.getValue(),optionWithMasses.getValue());
	for(PolynomialSet::const_iterator i=s.begin();i!=s.end();i++)g.push_back(*i);
      }

    if(optionDziobek.getValue())
      {
	PolynomialSet dz=DziobekEquations(g.getRing(),NOption.getValue(),optionWithMasses.getValue(),optionSVariables.getValue(),!optionLaurent.getValue());
	for(PolynomialSet::const_iterator i=dz.begin();i!=dz.end();i++)g.push_back(*i);
      }

    if(optionDeterminants.getValue()!=-1)
      {
	PolynomialSet de=nbodyDeterminants(g.getRing(),NOption.getValue(),optionWithMasses.getValue(),1+NOption.getValue()-optionDeterminants.getValue());
	for(PolynomialSet::const_iterator i=de.begin();i!=de.end();i++)g.push_back(*i);
      }
    if(optionMM.getValue())g.push_back(massEquation(g.getRing(),NOption.getValue(),optionWithMasses.getValue()));

    PP<<g.getRing();
    PP<<g;

    return 0;
  }
};

static NBodyApplication theApplication;
