#include "restrictedautoreduction.h"
#include "wallideal.h"
#include "tropical2.h"
#include "division.h"
#include <iostream>

void autoReduceRestricted(PolynomialSet *g, PolynomialSet const &initialIdeal)
{
  /*	  AsciiPrinter P(Stderr);
	  P.printPolynomialSet(*g);
	  P.printPolynomialSet(initialIdeal);
  */
	  assert(g->size()==initialIdeal.size());

  PolyhedralCone V=homogeneitySpace(initialIdeal);
  
  PolyhedralCone C=groebnerCone(*g,false);//no algebraic test - we hope that polyhedral methods are better
  PolyhedralCone C2=intersection(C,V);
  IntegerVector omega=C2.getRelativeInteriorPoint();

  while(1)
    {
      PolyhedralCone C=intersection(groebnerCone(*g,false),V);
      C.findFacets();
      IntegerVectorList normals=C.getHalfSpaces();
      
      bool allAreTrueFacets=true;
      for(IntegerVectorList::const_iterator i=normals.begin();i!=normals.end();i++)
	{
	  IntegerVectorList equations=C.getEquations();
	  equations.push_back(*i);
	  PolyhedralCone C2(C.getHalfSpaces(),equations,C.ambientDimension());
	  IntegerVector omega2=C2.getRelativeInteriorPoint();

	  //	  C2.print(&P);
	  bool isTrueFacet=false;
	  for(PolynomialSet::iterator j=g->begin();j!=g->end();j++)
	    {
	      Polynomial p=initialForm(*j,omega2)-initialForm(*j,omega);
	      if(!p.isZero())
		{
		  /*	  cerr << "poly1:";
		  P.printPolynomial(*j);
		  cerr << endl << "in_omega";
		  P.printPolynomial(initialForm(*j,omega));
		  cerr << endl << "in_omega2";
		  P.printPolynomial(initialForm(*j,omega2));
		  */
		  Polynomial q=divisionLift(p,initialIdeal,*g,LexicographicTermOrder(),true);
		  *j-=q;
		  /*cerr << endl << "q";
		  P.printPolynomial(q);
		  cerr << endl << "r";
		  P.printPolynomial(division(p,initialIdeal,LexicographicTermOrder()));
		  cerr << endl << "poly2:";
		  P.printPolynomial(*j);
		  cerr << endl;
		  */
		  if(!division(p,initialIdeal,LexicographicTermOrder()).isZero())isTrueFacet=true;
		}
	    }
	  if(!isTrueFacet)allAreTrueFacets=false;
	}
      
      if(allAreTrueFacets)break;
    }
}
