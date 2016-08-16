/*
 * app_librarytest.cpp
 *
 *  Created on: Sep 28, 2010
 *      Author: anders
 */

#include "gfanapplication.h"
#include "gfanlib.h"
#include "printer.h"
#include <iostream>
#include <fstream>
using namespace gfan;

class LibraryTestApplication : public GFanApplication
{
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program tests the gfan library.\n";
  }
  LibraryTestApplication()
  {
    registerOptions();
  }

  const char *name()
  {
    return "_librarytest";
  }

  int main()
  {
    int n=4;
    QVector s(n);
    QMatrix M(0,n);
    for(int i=0;i<n;i++)
      {
        QVector v=QVector::standardVector(n,i);
        std::cout << v;
        s-=v;
        M.append(QMatrix::rowVectorMatrix(s));
      }

    std::cout <<s<<M;
    M.reduce();
    std::cout <<s<<M;
    std::cerr<<"---------A"<<std::endl;





    ZMatrix A(3,2);

    A[0][0]=2;A[0][1]=2;
    A[1][0]=1;A[1][1]=2;
    A[2][0]=-2;A[2][1]=1;

    std::cout << ZCone::givenByRays(A,ZMatrix(0,2));


//    ZMatrix A(3,2);

    A[0][0]=2;A[0][1]=2;
    A[1][0]=1;A[1][1]=2;
    A[2][0]=-2;A[2][1]=1;

    ZMatrix temp(0,2);
    ZCone C(A,temp);

    std::cout<<C;
    C.canonicalize();
    std::cout<<C;

    std::cout<<"Relative interior point"<<endl<<C.getRelativeInteriorPoint()<<endl;
    std::cout<<"Extreme rays"<<endl<<C.extremeRays()<<endl;
    std::cout<<"Dual cone"<<endl<<C.dualCone()<<endl;
    std::cout<<"Unique point"<<endl<<C.getUniquePoint()<<endl;
    std::cout<<"Inequalities"<<endl<<C.getInequalities()<<endl;

    std::cout<<"Generators of span"<<endl<<C.generatorsOfSpan()<<endl;
    std::cout<<"Generators of lineality space"<<endl<<C.generatorsOfLinealitySpace()<<endl;

    Permutation a=Permutation::transposition(4,0,1);
    Permutation b=Permutation::cycle(4);
    SymmetryGroup G(4);
    G.computeClosure(a);
    G.computeClosure(b);

    std::cout<<G.size()<<":"<<G.orbitSize(ZVector::standardVector(4,0)+ZVector::standardVector(4,1))<<std::endl;



  /*  {
      std::cerr<<"TEST"<<std::endl;
      ZFan f(1);
//f.insert(ZCone::positiveOrthant(1));
      std::cerr<<f.toString();
      std::cerr<<"ENDTEST"<<std::endl;

    }
*/

    {
      SymmetryGroup sym(4);
      sym.computeClosure(Permutation::cycle(4));
      ZFan F(sym);

      F.insert(ZCone::positiveOrthant(4));
      std::cout<<F.toString();
      for(int i=0;i<10;i++)std::cerr<<F.numberOfConesOfDimension(i,false,false)<<std::endl;
      for(int i=0;i<4;i++)std::cerr<<F.getCone(3,i,false,false)<<std::endl;
std::cerr<<"AAA\n";
    }
    std::cerr<<"AAA\n";
    {
//      stringstream s;
//      std::string test="TEST";
 //     std::istringstream s(test);
      std::cerr<<"1AAA\n";
      std::fstream f;
      std::cerr<<"2AAA\n";
      f.open("fanfile");
      std::cerr<<"3AAA\n";
      ZFan G(f);
      std::cerr<<"4AAA\n";
      std::cout<<G.getAmbientDimension()<<std::endl;
      std::cerr<<"5AAA\n";
      std::cout<<G.toString();

      ZFan H=ZFan::fullFan(2);
//      std::cout<<ZFan::fullFan(2).toString();

      ZFan H2=H;
      std::cout<<H2.toString();
    }

    return 0;
  }
};

static LibraryTestApplication theApplication;
