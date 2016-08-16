#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "breadthfirstsearch.h"
#include "termorder.h"
#include "ep_standard.h"
#include "ep_xfig.h"
#include "gfanapplication.h"
#include "timer.h"
#include "log.h"
#include "matrix.h"
#include "lll.h"
#include "polyhedralfan.h"
#include "linalg.h"
#include "determinant.h"
#include "triangulation.h"
#include "intsinpolytope.h"
#include "graph.h"

#include "triangulation2.h"

#include "traverser_secondaryfan.h"
#include "symmetrictraversal.h"

#include <iostream>
#include <algorithm>

class SecondaryFanApplication : public GFanApplication
{
  SimpleOption hirschOption;
  SimpleOption searchOption;
  IntegerOption scaleOption;
  StringOption optionRestrictingFan;
  SimpleOption symmetryOption;
  SimpleOption optionIgnoreCones;
public:
  const char *helpText()
  {
    return "This program computes the secondary fan of a vector configuration. The configuration is given as an ordered list of vectors. In order to compute the secondary fan of a point configuration an additional coordinate of ones must be added. For example {(1,0),(1,1),(1,2),(1,3)}.\n";
  }
  SecondaryFanApplication():
    searchOption("--unimodular","Use heuristics to search for unimodular triangulation rather than computing the complete secondary fan"),
    scaleOption("--scale","Assuming that the first coordinate of each vector is 1, this option will take the polytope in the 1 plane and scale it. The point configuration will be all lattice points in that scaled polytope. The polytope must have maximal dimension. When this option is used the vector configuration must have full rank. This option may be removed in the future."),
    symmetryOption("--symmetry","Tells the program to read in generators for a group of symmetries (subgroup of $S_n$) after having read in the vector configuration. The program checks that the configuration stays fixed when permuting the variables with respect to elements in the group. The output is grouped according to the symmetry.\n"),
    optionRestrictingFan("--restrictingfan","Specify the name of a file containing a polyhedral fan in Polymake format. The computation of the Secondary fan will be restricted to this fan. If the --symmetry option is used then this restricting fan must be invariant under the symmetry and the orbits in the file must be with respect to the specified group of symmetries. The orbits of maximal cones of the file are then read in rather than the maximal cones.\n",0),
    optionIgnoreCones("--nocones","Tells the program not to output the CONES and MAXIMAL_CONES sections, but still output CONES_COMPRESSED and MAXIMAL_CONES_COMPRESSED if --symmetry is used."),
    hirschOption("--hirsch","")
  {
    hirschOption.hide();
    registerOptions();
  }

  const char *name()
  {
    return "_secondaryfan";
  }

  PolyhedralFan enumerate(Triangulation2 const &t)
  {
    PolyhedralFan ret(t.getN());
    list<Triangulation2> active;
    active.push_back(t);
    IntegerVectorList interiorPoints;
    interiorPoints.push_back(t.interiorPoint());
    while(!active.empty())
      {
	Triangulation2 a=active.front();
	PolyhedralCone C=a.secondaryCone();


	//	if(active.size()>100)break;//SLETMIGGGGG
	//log0 fprintf(stderr,"a\n");
	/*	{
	  PolyhedralCone C2=C;
	  C2.canonicalize();
	  }*/
	C.canonicalize();
	//log0 fprintf(stderr,"b\n");
	ret.insert(C);
	AsciiPrinter P(Stderr);
	//	C.print(&P);
	active.pop_front();
	//	fprintf(stderr,"pop\n");
	IntegerVectorList flips=a.facets();
	for(IntegerVectorList::const_iterator i=flips.begin();i!=flips.end();i++)
	  {
	    {
	      IntegerVectorList t=C.getEquations();
	      t.push_back(*i);
	      PolyhedralCone CF(C.getHalfSpaces(),t);
	      CF.findFacets();
	      //	      CF.canonicalize();
	    }

	    if(!i->isNonNegative()) //is this the right condition or should i be negated?
	    //   if(!(-*i).isNonNegative()) //is this the right condition or should i be negated?
	      {
		Triangulation2 b=a;
		log1 AsciiPrinter(Stderr)<<*i;
		/*fprintf(stderr,"Before:");
		  b.print(P);*/
		//		b.flip(*i);
		b.flipNew(-*i);
		/*fprintf(stderr,"After:");
		  b.print(P);*/
		if(!b.isEmpty())
		  {
		    IntegerVectorList inequalities=b.inequalities();
		    bool isKnown=false;
		    for(IntegerVectorList::const_iterator j=interiorPoints.begin();j!=interiorPoints.end();j++)
		      {
			bool match=true;
			for(IntegerVectorList::const_iterator k=inequalities.begin();k!=inequalities.end();k++)
			  {
			    if(dotLong(-*k,*j)<0)
			      {
				match=false;
				break;
			      }
			  }
			if(match)isKnown=true;
		      }
		    if(!isKnown)
		      {
			active.push_back(b);
			interiorPoints.push_back(b.interiorPoint());
		      }
		  }
	      }
	  }
      }
    return ret;
  }
  PolyhedralFan interactive(Triangulation2 const &t)
  {
    Triangulation2 a=t;

    while(1)
      {
	fprintf(stdout,"Triangles in current triangulation:%i\n",a.bases.size());
	PolyhedralCone C=a.secondaryCone();

	C.canonicalize();


	AsciiPrinter Pstd(Stderr);
	IntegerVectorList flips=a.facets();
	int I=0;
	for(IntegerVectorList::const_iterator i=flips.begin();i!=flips.end();i++,I++)
	  {
	    fprintf(stdout,"%i:\n",I);
	    Pstd.printVector(*i);


	    if(!i->isNonNegative())
	      {
		Triangulation2 b=a;
		//fprintf(stderr,"Before:");
		//b.print(P);
		//		b.flip(*i);
		/*		b.flipNew(*i);
		//fprintf(stderr,"After:");
		//b.print(P);
		fprintf(stdout,"Triangles in new triangulation:%i\n",b.bases.size());
		*/
	      }

	    fprintf(stdout,"\n");
	  }
	int s;
	//cin >> s;
	fscanf(stdin,"%i",&s);
	if((s>=0)&&(s<I))
	  {
	    IntegerVectorList::const_iterator i=flips.begin();
	    for(int I=0;I<s;I++)i++;
	    a.flipNew(*i);
	  }
	else
	  {
	    a.print(Pstd);
	  }
      }
    PolyhedralFan ret(0);

    return ret;
  }
  PolyhedralFan automatic(Triangulation2 const &t, int abortAtSize)
  {
    Triangulation2 a=t;

    cout << "Looking for triangulation with "<<abortAtSize<<" simplices."<<endl;

    while(1)
      {
	fprintf(stdout,"Triangles in current triangulation:%i\n",a.bases.size());
	//	PolyhedralCone C=a.secondaryCone();

	//	C.canonicalize();


	AsciiPrinter Pstd(Stderr);
	IntegerVectorList flips=a.facets();
	int I=0;
	cerr << "Number of facets of secondary cone: " << flips.size() <<endl;
	for(IntegerVectorList::const_iterator i=flips.begin();i!=flips.end();i++,I++)
	  {
	    if(!i->isNonNegative())
	      {
		Triangulation2 b=a;
		b.flipNew(-*i);
		fprintf(stdout,"Triangles in new triangulation:%i\n",b.bases.size());

		if(b.bases.size()==abortAtSize)
		  {
		    b.print(Pstd);

		    exit(0);
		  }

		if((b.bases.size()>a.bases.size())||((rand()&127)==0))
		  {
		    a=b;
		    break;
		  }
	      }
	  }
      }
    PolyhedralFan ret(0);

    return ret;
  }
  PolyhedralFan automaticHirsch(Triangulation2 const &t)
  {
    Triangulation2 a=t;
    while(1)
      {
	int nVertices=a.bases.size();
	int nEdges=a.coDimensionOneTriangles().size();
	int diameter=a.edgeGraph().diameter();
	int dimension=a.getD();
	int nFacets=a.usedRays().size();
	fprintf(stdout,"NVER: %i NEDGES: %i DIAMETER:%i DIMENSION:%i NFACETS:%i\n",nVertices,nEdges,diameter,dimension,nFacets);

	AsciiPrinter Pstd(Stderr);
	IntegerVectorList flips=a.facets();
	int I=0;
	float currentScore=a.hirschScore();
	cerr << "Current score: " << currentScore <<endl;
	for(IntegerVectorList::const_iterator i=flips.begin();i!=flips.end();i++,I++)
	  {
	    if(!i->isNonNegative())
	      {
		Triangulation2 b=a;
		b.flipNew(-*i);
		float bScore=b.hirschScore();
		fprintf(stdout,"New score:%f\n",bScore);

		if((bScore>currentScore)||((rand()&31)==0))
		  {
		    a=b;
		    break;
		  }
	      }
	  }
      }
    PolyhedralFan ret(0);

    return ret;
  }
  int main()
  {
    IntegerMatrix A=rowsToIntegerMatrix(FileParser(Stdin).parseIntegerVectorList()).transposed();

    int n=A.getWidth();

    SymmetryGroup s(n);
    if(symmetryOption.getValue())
      {
	IntegerVectorList generators=FileParser(Stdin).parseIntegerVectorList();
	for(IntegerVectorList::const_iterator i=generators.begin();i!=generators.end();i++)
	  {
	    assert(i->size()==n);
	    FieldMatrix M1=integerMatrixToFieldMatrix(A,Q);
	    FieldMatrix M2=integerMatrixToFieldMatrix(rowsToIntegerMatrix(SymmetryGroup::permuteIntegerVectorList(A.getRows(),*i)),Q);
	    M1.reduce();
	    M1.REformToRREform(true);
	    M2.reduce();
	    M2.REformToRREform(true);

	    if(!(M1==M2))
	      {
		AsciiPrinter(Stderr) << "Permutation "<< *i <<
		  " does not keep the configuration fixed.\n";
		assert(0);
	      }
	  }
	s.computeClosure(generators);
        s.createTrie();
      }


    if(scaleOption.getValue())
      {
	if(rank(A)!=A.getHeight())
	  {
	    cerr << "The vector configuration must have full rank in order to use the scale option.\n";
	    assert(0);
	  }

	int s=scaleOption.getValue();

	cout << "Input configuration:" << endl;
	AsciiPrinter(Stdout)<<A.transposed().getRows();

	IntegerVectorList empty;
	PolyhedralCone dual(A.transposed().getRows(),empty,n);
	dual.canonicalize();
	assert(dual.dimensionOfLinealitySpace()==0);
	IntegerVectorList inequalities=dual.extremeRays();

	IntegerMatrix M1=rowsToIntegerMatrix(inequalities);
	IntegerMatrix M2=-1*M1.submatrix(0,1,M1.getHeight(),M1.getWidth());
	IntegerVector rightHandSide=s*M1.transposed().submatrix(0,0,1,M1.getWidth())[0];
	IntegerVector v=s*A.submatrix(1,0,A.getHeight(),1).transposed()[0];

	IntegerVectorList l=intsInPolytopeGivenIneqAndPt(M2,rightHandSide,v);

	IntegerVectorList lT=rowsToIntegerMatrix(l).transposed().getRows();
	lT.push_front(IntegerVector::allOnes(l.size()));
	l=rowsToIntegerMatrix(lT).transposed().getRows();

	cout << "New configuration:" << endl;
	AsciiPrinter(Stdout)<<l;
	A=rowsToIntegerMatrix(l).transposed();
      }

    /* If the vector configuration A does not have full rank then
       change coordinates. */
    if(rank(A)!=A.getHeight())
      {
	FieldMatrix M=integerMatrixToFieldMatrix(A,Q);
	M.reduce(false,true);//force integer operations - preserving volume
	M.removeZeroRows();
	A=fieldMatrixToIntegerMatrix(M);
      }


    Triangulation2 t(A);

    /* Convert a Triangulation to a Triangulation2 */
    {
      list<Triangulation::Cone> T=Triangulation::triangulate(A.transposed());
      for(list<Triangulation::Cone>::const_iterator i=T.begin();i!=T.end();i++)
	{
	  IntegerVector v=i->size();
	  int J=0;
	  for(Triangulation::Cone::const_iterator j=i->begin();j!=i->end();j++,J++)
	    v[J]=*j;
	    t.bases.insert(v);
	}
    }

    if(searchOption.getValue())
      {
	PolyhedralFan f=automatic(t,t.totalVolume());
      }
    else if(hirschOption.getValue())
      {
	PolyhedralFan f=automaticHirsch(t);
      }
    else
{
		SymmetricTargetFanBuilder target(n,s);

		if(!optionRestrictingFan.getValue())
		{
			SecondaryFanTraverser traverser(t);
			symmetricTraverse(traverser,target,&s);
		}
		else
		{
		    PolyhedralFan f1=PolyhedralFan::readFan(optionRestrictingFan.getValue(),true,0,0,/*optionSymmetry.getValue()?&s:0*/0,false/*true*/);

		    for(PolyhedralFan::coneIterator i=f1.conesBegin();i!=f1.conesEnd();i++)
			{
		    	static int a;
		    	log2 cerr<<"Processing Cone "<<a++<<" which has dimension "<<i->dimension()<<endl;
		    	SecondaryFanTraverser traverser(triangulationWithFullDimensionalIntersection(t,*i),*i);
				symmetricTraverse(traverser,target,&s);
			}
		}

	        target.getFanRef().printWithIndices(&pout,
	                                    (symmetryOption.getValue()?FPF_group|FPF_conesCompressed:0)|
	                                    (optionIgnoreCones.getValue()?0:FPF_conesExpanded)|
	                                    FPF_maximalCones|FPF_cones,
	                                    &s);
/*		target.getFanRef().printWithIndices(&p,
											FPF_default|
											(symmetryOption.getValue()?FPF_group|FPF_conesCompressed:0),
											&s);
*/
/*    	PolyhedralFan f=enumerate(t).negated();//Changing sign

        f.printWithIndices(&p,
			   FPF_default|
			   (symmetryOption.getValue()?FPF_group|FPF_conesCompressed:0),
			   &s);
  */
		}
    return 0;
  }
};

static SecondaryFanApplication theApplication;

