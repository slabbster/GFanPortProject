#include "polyhedralcone.h"

#include "lp.h"
#include "subspace.h"
#include "symmetry.h"
#include "polymakefile.h"
#include "wallideal.h"
#include <sstream>
#include "triangulation.h"
#include "linalg.h"
#include "linalgfloat.h"
#include "continuedfractions.h"
#include "log.h"

//--------------------------------
// PolyhedralCone
//--------------------------------

static bool compareIntegerLists(IntegerVectorList const &a, IntegerVectorList const &b)
{
  assert(a.size()==b.size());

  IntegerVectorList::const_iterator B=b.begin();
  for(IntegerVectorList::const_iterator A=a.begin();A!=a.end();A++)
    {
      if(LexicographicTermOrder()(*A,*B))return true;
      if(LexicographicTermOrder()(*B,*A))return false;
      B++;
    }
  return false;
}

bool operator<(PolyhedralCone const &a, PolyhedralCone const &b)
{
  assert(a.state>=3);
  assert(b.state>=3);

  if(a.dimension()>b.dimension())return true;//INVERTED
  if(a.dimension()<b.dimension())return false;//INVERTED

  if(a.equations.size()<b.equations.size())return true;
  if(a.equations.size()>b.equations.size())return false;

  if(a.n<b.n)return true;
  if(a.n>b.n)return false;

  if(a.inequalities.size()<b.inequalities.size())return true;
  if(a.inequalities.size()>b.inequalities.size())return false;

  if(compareIntegerLists(a.equations,b.equations))return true;
  if(compareIntegerLists(b.equations,a.equations))return false;

  if(compareIntegerLists(a.inequalities,b.inequalities))return true;
  if(compareIntegerLists(b.inequalities,a.inequalities))return false;

  return false;
}

bool operator!=(PolyhedralCone const &a, PolyhedralCone const &b)
{
  return (a<b)||(b<a);
}


void PolyhedralCone::ensureStateAsMinimum(int s)
{
  if((state<1) && (s==1))
    {
      //log0 cerr<<"Number of inequalities: "<<halfSpaces.size()<<" Number of equations: "<<equations.size()<<endl;

/*      if(inequalities.size()==0)
        {
          equations=subsetBasis(equations);
        }
      else*/
      {
	IntegerMatrix temp=rowsToIntegerMatrix(equations,n);
	FieldMatrix m=integerMatrixToFieldMatrix(temp,Q);
	m.reduce();
	m.removeZeroRows();

	IntegerVectorList newInequalities;
	for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
	  {
	    FieldVector w=integerVectorToFieldVector(*i,Q);
	    w=m.canonicalize(w);
	    if(!w.isZero())
	      newInequalities.push_back(w.primitive());
	  }

        inequalities=newInequalities;
        removeDuplicates(inequalities);
        equations=fieldMatrixToIntegerMatrixPrimitive(m).getRows();
      }

      //      log0 cerr<<"Number of inequalities: "<<halfSpaces.size()<<" Number of equations: "<<equations.size()<<endl;

      //      log0 fprintf(Stderr,"+\n");
      if(!(preassumptions&PCP_impliedEquationsKnown))
      if(inequalities.size()>1)//there can be no implied equation unless we have at least two inequalities
        removeRedundantRows(&inequalities,&equations,false);
      //log0 cerr<<"Number of inequalities: "<<halfSpaces.size()<<" Number of equations: "<<equations.size()<<endl;
      //log0 fprintf(Stderr,"+\n");
    }
  if((state<2) && (s>=2) && !(preassumptions&PCP_facetsKnown))
    {
      //      halfSpaces.sort();
      //      AsciiPrinter(Stderr).printVectorList(halfSpaces);

       if(inequalities.size()>25)
	 //	if(0)
	 {
	  IntegerVectorList h1;
	  IntegerVectorList h2;
	  bool a=false;
	  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
	    {
	      if(a)
		h1.push_back(*i);
	      else
		h2.push_back(*i);
	      a=!a;
	    }
	  PolyhedralCone c1(h1,equations);
	  PolyhedralCone c2(h2,equations);
	  c1.ensureStateAsMinimum(2);
	  c2.ensureStateAsMinimum(2);
	  inequalities=c1.inequalities;
	  for(IntegerVectorList::const_iterator i=c2.inequalities.begin();i!=c2.inequalities.end();i++)
	    inequalities.push_back(*i);
	}


       //      fprintf(Stderr,"Number half spaces: %i, number of equations: %i\n",halfSpaces.size(),equations.size());
      if(equations.size())
	{
	  FieldMatrix equationSpace=integerMatrixToFieldMatrix(rowsToIntegerMatrix(equations,n),Q);
	  equationSpace.reduce();
	  IntegerVectorList halfSpaces2;
	  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
	    {
	      halfSpaces2.push_back(equationSpace.canonicalize(integerVectorToFieldVector(*i,Q)).primitive());
	    }
	  for(IntegerVectorList::const_iterator i=halfSpaces2.begin();i!=halfSpaces2.end();i++)
	    if((i->max()>1000) || (i->min()<-1000))goto fallBack;//more overflows caught in lp_cdd.
	  inequalities=fastNormals(halfSpaces2);
	  goto noFallBack;
	fallBack:
	  removeRedundantRows(&inequalities,&equations,true);
	noFallBack:;
	}
      else
	inequalities=fastNormals(inequalities);
      //   fprintf(Stderr,"done\n");
    }
  if((state<3) && (s>=3))
    {
      //      fprintf(Stderr,"Number half spaces: %i, number of equations: %i\n",halfSpaces.size(),equations.size());
      //      log1 fprintf(Stderr,"Computing subspace...\n");
      Subspace v(equations,n);

      equations=v.getRepresentation();
      //      log1 fprintf(Stderr,"...done computing subspace.\n");


      for(IntegerVectorList::iterator i=inequalities.begin();i!=inequalities.end();i++)
	{
	  *i=v.canonicalizeVector(*i);
	}
      inequalities.sort();
      //fprintf(Stderr,"done\n");
    }
  state=s;
}

void PolyhedralCone::canonicalize()
{
  ensureStateAsMinimum(3);
}


void PolyhedralCone::findFacets()
{
  ensureStateAsMinimum(2);
}

void PolyhedralCone::findImpliedEquations()
{
  ensureStateAsMinimum(1);
}

/*PolyhedralCone::PolyhedralCone(IntegerVectorList const &halfSpaces_, IntegerVectorList const &equations_, int ambientDimension, int state):
	inequalities(halfSpaces_),
	equations(equations_),
	state(0),
	multiplicity(1),
	n(ambientDimension)
{
	this->state=state;
}
*/

PolyhedralCone::PolyhedralCone(int ambientDimension):
  n(ambientDimension),
  state(1),
  preassumptions(PCP_impliedEquationsKnown|PCP_facetsKnown),
  multiplicity(1),
  haveExtremeRaysBeenCached(false)
{
}


PolyhedralCone::PolyhedralCone(IntegerVectorList const &halfSpaces_, IntegerVectorList const &equations_, int ambientDimension, int preassumptions_):
//PolyhedralCone::PolyhedralCone(IntegerVectorList const &halfSpaces_, IntegerVectorList const &equations_, int ambientDimension):
  inequalities(halfSpaces_),
  equations(equations_),
//  equations(subsetBasis(equations_)),
  state(0),
  preassumptions(preassumptions_),
  multiplicity(1),
  haveExtremeRaysBeenCached(false)
{
  n=ambientDimension;
  if(n==-1)
    {
      if(!halfSpaces_.empty())
	n=halfSpaces_.begin()->size();
      else if(!equations_.empty())
	n=equations_.begin()->size();
      else
	{
	  assert(0);
	}
    }
  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
    {
      assert(i->size()==n);
    }
  for(IntegerVectorList::const_iterator i=equations.begin();i!=equations.end();i++)
    {
      /*AsciiPrinter P(Stderr);
      P.printString("FJDSKFJA\n");
      P.printVector(*i);
      */
      assert(i->size()==n);
    }
  ensureStateAsMinimum(1);
  //  computeAndReduceLinearitySpace();
}


IntegerVector PolyhedralCone::getRelativeInteriorPoint()const
{
  //  ensureStateAsMinimum(1);
  assert(state>=1);
  IntegerVector ret;
  IntegerVectorList g=equations;
  int numberOfEqualities=g.size();
  g.insert(g.end(),inequalities.begin(),inequalities.end());
  IntegerVector equalitySet(g.size());
  for(int i=0;i<numberOfEqualities;i++)equalitySet[i]=1;

  if(!g.empty())
    ret=relativeInteriorPoint(ambientDimension(),g,&equalitySet);
  else
    ret=IntegerVector(n);//cone is the full space, lp code would fail, since the dimension is unknown

  for(IntegerVectorList::const_iterator i=equations.begin();i!=equations.end();i++)
    {
      assert(dotLong(*i,ret)==0);
    }

  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
    {
      if(!(dotLong(*i,ret)>0))
	{
	  fprintf(Stderr,"PolyhedralCone::relativeInteriorPoint() : halfSpaces not reduced or mistake in cdd interface!!!\n");
	}
    }

  return ret;
}


IntegerVectorList PolyhedralCone::getHalfSpaces()const
{
  return inequalities;
}


const IntegerVectorList &PolyhedralCone::getEquations()const
{
  assert(state>=1);
  return equations;
}


IntegerVectorList PolyhedralCone::generatorsOfSpan()const
{
	assert(isInStateMinimum(1));
	IntegerVectorList empty;
	PolyhedralCone temp(empty,getEquations(),n);
	return temp.dualCone().getEquations();
}


IntegerVectorList PolyhedralCone::generatorsOfLinealitySpace()const
{
  IntegerVectorList l=equations;
  l.insert(l.end(),inequalities.begin(),inequalities.end());
  FieldMatrix L=integerMatrixToFieldMatrix(rowsToIntegerMatrix(l,ambientDimension()),Q);
  return fieldMatrixToIntegerMatrixPrimitive(L.reduceAndComputeKernel()).getRows();
//  return linealitySpace().generatorsOfSpan();
}


int PolyhedralCone::ambientDimension()const
{
  return n;
}


int PolyhedralCone::codimension()const
{
  return ambientDimension()-dimension();
  //  return getEquations().size();
}


int PolyhedralCone::dimension()const
{
  assert(state>=1);
  //  ensureStateAsMinimum(1);
  return ambientDimension()-equations.size();
}


bool PolyhedralCone::isZero()const
{
  return dimension()==0;
}


bool PolyhedralCone::isFullSpace()const
{
	for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
		if(!i->isZero())return false;
	for(IntegerVectorList::const_iterator i=equations.begin();i!=equations.end();i++)
		if(!i->isZero())return false;
	return true;
}


PolyhedralCone intersection(const PolyhedralCone &a, const PolyhedralCone &b)
{
  assert(a.ambientDimension()==b.ambientDimension());
  IntegerVectorList inequalities=a.inequalities;
  inequalities.insert(inequalities.end(),b.inequalities.begin(),b.inequalities.end());
  IntegerVectorList equations=a.equations;


  equations.insert(equations.end(),b.equations.begin(),b.equations.end());
  {
    removeDuplicates(equations);
    removeDuplicates(inequalities);
    IntegerVectorList Aequations=a.equations;
    IntegerVectorList Ainequalities=a.inequalities;
    removeDuplicates(Aequations);
    removeDuplicates(Ainequalities);
    if((Ainequalities.size()==inequalities.size()) && (Aequations.size()==equations.size()))return a;
    IntegerVectorList Bequations=b.equations;
    IntegerVectorList Binequalities=b.inequalities;
    removeDuplicates(Bequations);
    removeDuplicates(Binequalities);
    if((Binequalities.size()==inequalities.size()) && (Bequations.size()==equations.size()))return b;
  }

  return PolyhedralCone(inequalities,equations,a.ambientDimension());
}

PolyhedralCone product(const PolyhedralCone &a, const PolyhedralCone &b)
{
  IntegerVectorList equations2;
  IntegerVectorList inequalities2;

  int n1=a.n;
  int n2=b.n;

  for(IntegerVectorList::const_iterator i=a.equations.begin();i!=a.equations.end();i++)
    equations2.push_back(concatenation(*i,IntegerVector(n2)));
  for(IntegerVectorList::const_iterator i=b.equations.begin();i!=b.equations.end();i++)
    equations2.push_back(concatenation(IntegerVector(n1),*i));
  for(IntegerVectorList::const_iterator i=a.inequalities.begin();i!=a.inequalities.end();i++)
    inequalities2.push_back(concatenation(*i,IntegerVector(n2)));
  for(IntegerVectorList::const_iterator i=b.inequalities.begin();i!=b.inequalities.end();i++)
    inequalities2.push_back(concatenation(IntegerVector(n1),*i));

  PolyhedralCone ret(inequalities2,equations2,n1+n2);
  ret.setMultiplicity(a.getMultiplicity()*b.getMultiplicity());
  ret.setLinearForm(concatenation(a.getLinearForm(),b.getLinearForm()));

  ret.ensureStateAsMinimum(a.state);
  ret.ensureStateAsMinimum(b.state);

  return ret;
}

PolyhedralCone PolyhedralCone::positiveOrthant(int dimension)
{
  IntegerVectorList halfSpaces;

  for(int i=0;i<dimension;i++)halfSpaces.push_back(IntegerVector::standardVector(dimension,i));

  IntegerVectorList empty;
  return PolyhedralCone(halfSpaces,empty,dimension);
}


bool PolyhedralCone::isInStateMinimum(int s)const
{
  return state>=s;
}

int PolyhedralCone::getState()const
{
  return state;
}


void PolyhedralCone::print(class Printer *p, bool xml)const
{
  if(0)
    {
      p->printString("Printing PolyhedralCone");
      p->printNewLine();
      p->printString("Ambient dimension: ");
      p->printInteger(n);
      p->printNewLine();
      if(isInStateMinimum(1))
	{
	  p->printString("Dimension: ");
	  p->printInteger(dimension());
	  p->printNewLine();
	}
      p->printString("Linearity space:");
      //  p->printNewLine();
      p->printVectorList(equations);
      p->printString("Inequalities:");
      p->printVectorList(inequalities);
      p->printString("Relative interior point:\n");
      p->printVector(getRelativeInteriorPoint());
      p->printNewLine();
      p->printString("Done printing PolyhedralCone.");
      p->printNewLine();
    }
  else
    {
      PolymakeFile polymakeFile;
      polymakeFile.create("NONAME","PolyhedralCone","PolyhedralCone",xml);
      polymakeFile.writeCardinalProperty("AMBIENT_DIM",n);
      if(isInStateMinimum(1))
	{
	  polymakeFile.writeCardinalProperty("DIM",dimension());
	  //need to check that the following is done correctly
	  //       	  polymakeFile.writeCardinalProperty("LINEALITY_DIM",dimensionOfLinealitySpace());

	  polymakeFile.writeMatrixProperty("IMPLIED_EQUATIONS",rowsToIntegerMatrix(equations,n));
	}
      polymakeFile.writeCardinalProperty("LINEALITY_DIM",dimensionOfLinealitySpace());
      polymakeFile.writeMatrixProperty("LINEALITY_SPACE",rowsToIntegerMatrix(linealitySpace().dualCone().getEquations(),n));


      if(isInStateMinimum(2))
	polymakeFile.writeMatrixProperty("FACETS",rowsToIntegerMatrix(inequalities,n));
      else
	polymakeFile.writeMatrixProperty("INEQUALITIES",rowsToIntegerMatrix(inequalities,n));

      polymakeFile.writeCardinalVectorProperty("RELATIVE_INTERIOR_POINT",getRelativeInteriorPoint());


      stringstream s;
      polymakeFile.writeStream(s);
      string S=s.str();
      p->printString(S.c_str());
    }
}


static IntegerVector dehomogenize(IntegerVector const &v)
{
  assert(v.size()>0);

  IntegerVector ret(v.size()-1);

  for(int i=0;i<v.size()-1;i++)ret[i]=v[i];
  return ret;
}

static IntegerVectorList dehomogenize(IntegerVectorList const &l)
{
  IntegerVectorList ret;

  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      ret.push_back(dehomogenize(*i));
    }
  return ret;
}

PolyhedralCone PolyhedralCone::withLastCoordinateRemoved()const
{
  assert(n>0);

  return PolyhedralCone(dehomogenize(inequalities),dehomogenize(equations));
}

bool PolyhedralCone::containsPositiveVector()const
{
  PolyhedralCone temp=intersection(*this,PolyhedralCone::positiveOrthant(n));
  IntegerVector v=temp.getRelativeInteriorPoint();
  return v.isPositive();
}


int PolyhedralCone::dimensionOfLinealitySpace()const
{
  if(inequalities.empty())return dimension();
  IntegerVectorList a;
  PolyhedralCone temp(a,inequalities);
  temp=intersection(temp,*this);

  return temp.dimension();
}


bool PolyhedralCone::contains(IntegerVector const &v)const
{
  for(IntegerVectorList::const_iterator i=equations.begin();i!=equations.end();i++)
    {
      if(dotLong(*i,v)!=0)return false;
    }

  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
    {
      if(dotLong(*i,v)<0)return false;
    }

  return true;
}


bool PolyhedralCone::contains(IntegerVectorList const &l)const
{
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    if(!contains(*i))return false;
  return true;
}

bool PolyhedralCone::contains(PolyhedralCone const &c)const
{
  PolyhedralCone c2=intersection(*this,c);
  PolyhedralCone c3=c;
  c2.canonicalize();
  c3.canonicalize();
  return !(c2!=c3);
}


bool PolyhedralCone::containsRelatively(IntegerVector const &v)const
{
  assert(state>=1);
  for(IntegerVectorList::const_iterator i=equations.begin();i!=equations.end();i++)
    {
      if(dotLong(*i,v)!=0)return false;
    }

  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
    {
      if(dotLong(*i,v)<=0)return false;
    }

  return true;
}


PolyhedralCone PolyhedralCone::permuted(IntegerVector const &v)const
{
  PolyhedralCone ret(SymmetryGroup::permuteIntegerVectorList(inequalities,v),SymmetryGroup::permuteIntegerVectorList(equations,v),n);
  if(state>=1)ret.state=1;
  if(state>=2)ret.state=2;

  ret.ensureStateAsMinimum(state);
  ret.setMultiplicity(multiplicity);
  return ret;
}

IntegerVector PolyhedralCone::getUniquePoint()const
{
  IntegerVectorList rays=extremeRays();
  IntegerVector ret(n);
  for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
    ret+=*i;

  assert(containsRelatively(ret));//remove this check later
  return ret;
}

IntegerVector PolyhedralCone::getUniquePointFromExtremeRays(IntegerVectorList const &extremeRays)const
{
	IntegerVector ret(n);
	for(IntegerVectorList::const_iterator i=extremeRays.begin();i!=extremeRays.end();i++)
		if(contains(*i))ret+=*i;
	return ret;
}

/**
 * Returns a primitive vector parallel to the projection of v onto the orthogonal complement of E.
 */
static IntegerVector primitiveProjection(IntegerVectorList const &E, IntegerVector &v, bool useFloat)
{
	int n=v.size();
if(useFloat)
	{
	IntegerMatrix E2=rowsToIntegerMatrix(E,n);
	linalgfloat::Matrix E3(E.size(),n);
	for(int i=0;i<E3.getHeight();i++)
		for(int j=0;j<E3.getWidth();j++)
			E3[i][j]=E2[i][j];
	cerr<<E3;
	E3.orthogonalize();
	linalgfloat::Vector v2(n);
	for(int i=0;i<n;i++)v2[i]=v[i];
	cerr<<E3;
	linalgfloat::Vector coefficients=E3.projectionCoefficients(v2);

	cerr<<coefficients<<"\n";

	vector<double> temp(coefficients.size());
	for(int i=0;i<coefficients.size();i++)temp[i]=coefficients[i];
	vector<int> tempInt(coefficients.size());
	int denominator;
	doubleVectorToFractions(temp, tempInt, denominator);
	IntegerVector tempInt2(tempInt.size());for(int i=0;i<tempInt.size();i++)tempInt2[i]=tempInt[i];
	if(tempInt2.max()>1000)goto fallback;
	if(tempInt2.min()<-1000)goto fallback;
	if(denominator>1000)goto fallback;
	if(denominator<-1000)goto fallback;
	IntegerVector ret=denominator*v;
	int I=0;
	for(IntegerVectorList::const_iterator i=E.begin();i!=E.end();i++,I++)ret-=coefficients[I]*(*i);
	if(!E2.inKernel(ret))goto fallback;
	return normalized(ret);
	}
	fallback:
	debug<<"FLOATCONELINALG FALLBACK\n";

	IntegerVectorList inequalities;
	inequalities.push_back(v);
	FieldMatrix linealitySpaceOrth=combineOnTop(integerMatrixToFieldMatrix(rowsToIntegerMatrix(E,n),Q),integerMatrixToFieldMatrix(rowsToIntegerMatrix(inequalities,n),Q));
	FieldMatrix temp=combineOnTop(linealitySpaceOrth.reduceAndComputeKernel(),integerMatrixToFieldMatrix(rowsToIntegerMatrix(E,n),Q));
	FieldMatrix temp2=temp.reduceAndComputeKernel();
	assert(temp2.getHeight()==1);
	  return temp2[0].primitive();
}

IntegerVectorList PolyhedralCone::extremeRays(IntegerVectorList const *generatorsOfLinealitySpace)const
{
  assert((dimension()==ambientDimension()) || (state>=3));

  if(haveExtremeRaysBeenCached)return cachedExtremeRays;
  //	 this->print(&debug);
  /* This code actually works even for lower dimensional cones if they have been canonicalized. */

  //  AsciiPrinter(Stderr).printVectorList(halfSpaces);
  //  AsciiPrinter(Stderr).printVectorList(equations);

  IntegerVectorList ret;

  // log0 fprintf(Stderr,"calling double description (cddlib)\n");
  IntegerVectorList indices=extremeRaysInequalityIndices(inequalities);
  //  log0 fprintf(Stderr,"returning\n");
  // log0 fprintf(Stderr,"computing relative interior points\n");

//  debug<<"INDICES"<<indices;
  for(IntegerVectorList::const_iterator i=indices.begin();i!=indices.end();i++)
    {
//        log0 AsciiPrinter(Stderr)<<*i;
      if(1)
	{
      /* At this point we know lineality space, implied equations and
	 also inequalities for the ray. To construct a vector on the
	 ray which is stable under (or indendent of) angle and
	 linarity preserving transformation we find the dimension 1
	 subspace orthorgonal to the implied equations and the
	 lineality space and pick a suitable primitive generator */

    	  /* To be more precise,
    	   * let E be the set of equations, and v the inequality defining a  ray R.
    	   * We wish to find a vector satisfying these, but it must also be orthogonal
    	   * to the lineality space of the cone, that is, in the span of {E,v}.
    	   * One way to get such a vector is to project v to E an get a vector p.
    	   * Then v-p is in the span of {E,v} by construction.
    	   * The vector v-p is also in the orthogonal complement to E by construction,
    	   * that is, the span of R.
    	   * We wish to argue that it is not zero.
    	   * That would imply that v=p, meaning that v is in the span of the equations.
    	   * However, that would contradict that R is a ray.
    	   * In case v-p does not satisfy the inequality v (is this possible?), we change the sign.
    	   *
    	   * As a consequence we need the following procedure
    	   * primitiveProjection():
    	   *    Input: E,v
    	   *    Output: A primitive representation of the vector v-p, where p is the projection of v onto E
    	   *
    	   * Notice that the output is a Q linear combination of the input and that p is
    	   * a linear combination of E. The check that p has been computed correctly,
    	   * it suffices to check that v-p satisfies the equations E.
    	   * The routine will actually first compute a multiple of v-p.
    	   * It will do this using floating point arithmetics. It will then transform
    	   * the coefficients to get the multiple of v-p into integers. Then it
    	   * verifies in exact arithmetics, that with these coefficients we get a point
    	   * satisfying E. It then returns the primitive vector on the ray v-p.
    	   * In case of a failure it falls back to an implementation using rational arithmetics.
    	   */


	  IntegerVector asVector(inequalities.size());
    //  log0 AsciiPrinter(Stderr)<<asVector;
	  for(int j=0;j<i->size();j++)asVector[(*i)[j]]=1;
     // log0 AsciiPrinter(Stderr)<<asVector;

	  IntegerVectorList equations=this->equations;
	  IntegerVectorList inequalities;

	  IntegerVector theInequality;

    //  log0 AsciiPrinter(Stderr)<<inequalities;

      IntegerVectorList::const_iterator a=this->inequalities.begin();
	  for(int j=0;j<asVector.size();j++,a++)
	    if(asVector[j])
	      equations.push_back(*a);
	    else
	    	theInequality=*a;

	  assert(!theInequality.isZero());

	  IntegerVector thePrimitiveVector;
	  if(generatorsOfLinealitySpace)
	  {
		IntegerMatrix temp=rowsToIntegerMatrix(equations,n);
		temp.append(rowsToIntegerMatrix(*generatorsOfLinealitySpace,n));

//		debug<<*generatorsOfLinealitySpace;
//			debug.printVectorList(temp.getRows());
		thePrimitiveVector=vectorInKernel(temp);
	  }
	  else
	  {
	  //log0  AsciiPrinter(Stderr)<<equations;
/*		  {
			  IntegerVectorList equations2=this->equations;
			  for(IntegerVectorList::const_iterator i=equations.begin();i!=equations.end();i++)equations2.push_back(*i);
		  debug<<primitiveProjection(equations2,theInequality,true)<<"\n";
		  debug<<primitiveProjection(equations2,theInequality,false)<<"\n";
		  }
*/
		  /** TODO: These calls are slow, but used often in the symmetric traversal. Maybe they should be done in floating point somehow.*/
	  FieldMatrix linealitySpaceOrth=combineOnTop(integerMatrixToFieldMatrix(rowsToIntegerMatrix(this->equations,n),Q),integerMatrixToFieldMatrix(rowsToIntegerMatrix(this->inequalities,n),Q));
	  FieldMatrix temp=combineOnTop(linealitySpaceOrth.reduceAndComputeKernel(),integerMatrixToFieldMatrix(rowsToIntegerMatrix(equations,n),Q));
	  FieldMatrix temp2=temp.reduceAndComputeKernel();

	  assert(temp2.getHeight()==1);
	  thePrimitiveVector=temp2[0].primitive();

	//  debug<<thePrimitiveVector<<"\n";
	  }
	  if(!contains(thePrimitiveVector))thePrimitiveVector=-thePrimitiveVector;
	  ret.push_back(thePrimitiveVector);
	}
      else
	{
  /*      IntegerVectorList equations;
      for(int j=0;j<i->size();j++)
	{
	  IntegerVectorList::const_iterator a=halfSpaces.begin();
	  for(int k=0;k<(*i)[j];k++)
	    {
	      assert(a!=halfSpaces.end());
	      	      a++;
	    }
	  assert(a!=halfSpaces.end());
	  equations.push_back(*a);
	  }*/

	  IntegerVector asVector(inequalities.size());
	  for(int j=0;j<i->size();j++)asVector[(*i)[j]]=1;

	  IntegerVectorList equations=this->equations;
	  IntegerVectorList inequalities;

	  IntegerVectorList::const_iterator a=inequalities.begin();
	  for(int j=0;j<asVector.size();j++,a++)
	    if(asVector[j])
	      equations.push_back(*a);
	    else
	      inequalities.push_back(*a);


	  //log0 fprintf(Stderr,"Equations %i, HalfSpaces: %i\n",equations.size(),inequalities.size());
	  IntegerVector u=PolyhedralCone(inequalities,equations).getRelativeInteriorPoint();
	  if(!u.isZero())
	    ret.push_back(u);
	  else
	    {
	      log2 fprintf(Stderr,"Remember to fix cdd double description interface\n");
	    }
	}
    }
  // log0 fprintf(Stderr,"done computing relative interior points\n");


  /*   //triangulation method. Keep this code.
  {
    IntegerMatrix temp=rowsToIntegerMatrix(halfSpaces);
 log0 fprintf(Stderr,"calling double description (triangulation)\n");
    IntegerVectorList ret2=Triangulation::normals(temp);
 log0 fprintf(Stderr,"returning\n");
 return ret2;

    AsciiPrinter(Stderr).printVectorList(halfSpaces);
    AsciiPrinter(Stderr).printVectorList(equations);
    fprintf(Stderr,"dim:%i\n",dimension());
    //ret.sort();
    AsciiPrinter(Stderr).printVectorList(ret);
    //  AsciiPrinter(Stderr).printVectorList(fastNormals(ret2));
    AsciiPrinter(Stderr).printVectorList(ret2);

    fprintf(stderr,"-----------------------\n");
  }
  */

  cachedExtremeRays=ret;
  haveExtremeRaysBeenCached=true;

  return ret;
}


bool PolyhedralCone::isSimplicial()const
{
  assert(state>=2);

  //  ensureStateAsMinimum(2);
  //  AsciiPrinter P(Stderr);
  //  print(&P);
  return codimension()+getHalfSpaces().size()+dimensionOfLinealitySpace()==n;
}


bool PolyhedralCone::checkDual(PolyhedralCone const &c)const
{
  assert(dimensionOfLinealitySpace()+c.dimension()==ambientDimension());

  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
    {
      assert(c.contains(*i));
    }
  for(IntegerVectorList::const_iterator i=equations.begin();i!=equations.end();i++)
    {
      assert(c.contains(*i));
    }
  return true;
}


PolyhedralCone PolyhedralCone::dualCone()const
{
  assert(state>=1);

  IntegerVectorList dualInequalities,dualEquations;

  dual(ambientDimension(),inequalities,equations,&dualInequalities,&dualEquations);

  PolyhedralCone ret(dualInequalities,dualEquations);

  ret.ensureStateAsMinimum(state);
  //  ret.canonicalize();


  assert(checkDual(ret));
  assert(ret.checkDual(*this));

  return ret;
}


PolyhedralCone PolyhedralCone::negated()const
{
  IntegerVectorList inequalities2;
  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)inequalities2.push_back(-*i);
//  PolyhedralCone ret(inequalities2,equations,n);
  PolyhedralCone ret(inequalities2,equations,n,(areFacetsKnown()?PCP_facetsKnown:0)|(areImpliedEquationsKnown()?PCP_impliedEquationsKnown:0));
  ret.ensureStateAsMinimum(state);
  return ret;
}

PolyhedralCone PolyhedralCone::linealitySpace()const
{
  IntegerVectorList l1=getEquations();
  IntegerVectorList l2=getHalfSpaces();

  l1.splice(l1.begin(),l2);

  IntegerVectorList temp;
  PolyhedralCone ret(temp,l1,n);
  ret.ensureStateAsMinimum(state);
  return ret;
}


int PolyhedralCone::getMultiplicity()const
{
  return multiplicity;
}


void PolyhedralCone::setMultiplicity(int m)
{
  multiplicity=m;
}


IntegerVectorList PolyhedralCone::quotientLatticeBasis()const
{
  assert(isInStateMinimum(1));// Implied equations must have been computed in order to know the span of the cone

  int a=equations.size();
  int b=inequalities.size();

  // Implementation below could be moved to nonLP part of code.

  // small vector space defined by a+b equations.... big by a equations.

  FieldMatrix M=combineLeftRight(combineLeftRight(
						  integerMatrixToFieldMatrix(rowsToIntegerMatrix(equations,n),Q).transposed(),
						  integerMatrixToFieldMatrix(rowsToIntegerMatrix(inequalities,n),Q).transposed()),
				 FieldMatrix::identity(Q,n));
  M.reduce(false,true);
  /*
    [A|B|I] is reduced to [A'|B'|C'] meaning [A'|B']=C'[A|B] and A'=C'A.

    [A'|B'] is in row echelon form, implying that the rows of C' corresponding to zero rows
    of [A'|B'] generate the lattice cokernel of [A|B] - that is the linealityspace intersected with Z^n.

    [A'] is in row echelon form, implying that the rows of C' corresponding to zero rows of [A'] generate
    the lattice cokernel of [A] - that is the span of the cone intersected with Z^n.

    It is clear that the second row set is a superset of the first. Their difference is a basis for the quotient.
   */
  IntegerVectorList ret;

  for(int i=0;i<M.getHeight();i++)
    if(M[i].subvector(0,a).isZero()&&!M[i].subvector(a,a+b).isZero())
      {
        ret.push_back(fieldVectorToIntegerVector(M[i].subvector(a+b,a+b+n)));
      }

  return ret;
}


IntegerVector PolyhedralCone::semiGroupGeneratorOfRay()const
{
  IntegerVectorList temp=quotientLatticeBasis();
  assert(temp.size()==1);
  for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
    if(dotLong(temp.front(),*i)<0)
      {
	temp.front()=-temp.front();
	break;
      }
  return temp.front();
}


IntegerVector const &PolyhedralCone::getLinearForm()const
{
  return linearForm;
}


void PolyhedralCone::setLinearForm(IntegerVector const &linearForm_)
{
  linearForm=linearForm_;
}


PolyhedralCone PolyhedralCone::link(IntegerVector const &w)const
{
	/*
	 * Observe that the inequalities giving rise to facets
	 * also give facets in the link, if they are kept as
	 * inequalities. This means that the state cannot decrease when taking links.
	 *
	 */
//  assert(state>=3);
  IntegerVectorList inequalities2;
  for(IntegerVectorList::const_iterator j=inequalities.begin();j!=inequalities.end();j++)
    if(dotLong(w,*j)==0)inequalities2.push_back(*j);
//  PolyhedralCone C(inequalities2,getEquations(),n);
//  C.canonicalize();
//  PolyhedralCone C(inequalities2,getEquations(),n,state);//STATE-----------------------------------------------------
  PolyhedralCone C(inequalities2,getEquations(),n,(areImpliedEquationsKnown()?PCP_impliedEquationsKnown:0)|(areFacetsKnown()?PCP_facetsKnown:0));
  C.ensureStateAsMinimum(state);

  C.setLinearForm(getLinearForm());
  C.setMultiplicity(getMultiplicity());

  return C;
}


PolyhedralCone PolyhedralCone::givenByRays(IntegerVectorList const &generators, IntegerVectorList const &linealitySpace, int n)
{
  //rewrite modulo lineality space
  IntegerVectorList newGenerators=generators;
  {
    Subspace l(linealitySpace,n);
    for(IntegerVectorList::const_iterator i=generators.begin();i!=generators.end();i++)
      newGenerators.push_back(l.canonicalizeVector(*i));
  }

  PolyhedralCone dual(newGenerators,linealitySpace,n);
  dual.findFacets();
  dual.canonicalize();
  IntegerVectorList inequalities=dual.extremeRays();

  IntegerVectorList span=generators;
  for(IntegerVectorList::const_iterator i=linealitySpace.begin();i!=linealitySpace.end();i++)span.push_back(*i);
  FieldMatrix m2Q=integerMatrixToFieldMatrix(rowsToIntegerMatrix(span,n),Q);
  IntegerVectorList equations=fieldMatrixToIntegerMatrixPrimitive(m2Q.reduceAndComputeKernel()).getRows();

  return PolyhedralCone(inequalities,equations,n);
}


FieldElement PolyhedralCone::volume()const
{
  AsciiPrinter P(Stderr);

  PolyhedralCone Ctemp=intersection(*this,this->linealitySpace().dualCone());

  Ctemp.canonicalize();
  cerr<<"testestests";
  IntegerVectorList eq2=rowsToIntegerMatrix(Ctemp.equations,n).transposed().getRows();
  eq2.push_front(IntegerVector(Ctemp.equations.size()));
  eq2=rowsToIntegerMatrix(eq2).transposed().getRows();
  cerr<<"testestests2";

  IntegerVectorList in2=rowsToIntegerMatrix(Ctemp.inequalities,n).transposed().getRows();
  in2.push_front(IntegerVector(Ctemp.inequalities.size()));
  in2=rowsToIntegerMatrix(in2).transposed().getRows();

  for(int i=0;i<n;i++)
    {
      in2.push_back(IntegerVector::standardVector(n+1,i+1)+IntegerVector::standardVector(n+1,0));
      //in2.push_back(-IntegerVector::standardVector(n+1,i+1)+IntegerVector::standardVector(n+1,0));
    }

  PolyhedralCone lifted(in2,eq2,n+1);

  lifted.canonicalize();
  IntegerMatrix A=rowsToIntegerMatrix(lifted.extremeRays()).transposed();//points are columns

  cerr << "height " << A.getHeight() << " width" <<A.getWidth()<<endl;

  P<<A.transposed().getRows();

  FieldMatrix A2(Q,A.getHeight()-1,A.getWidth());
  for(int i=0;i<A.getHeight()-1;i++)
    {
      A2[i]=integerVectorToFieldVector(A[i+1],Q)/integerVectorToFieldVector(A[0],Q);
    }
  A2=A2.transposed();//Now points are rows

  P<<"Triangulating\n";
  list<Triangulation::Cone> T=Triangulation::triangulate(A.transposed(),true);//revlex
  P<<"Done triangulating\n";

  FieldElement ret(Q);

  for(list<Triangulation::Cone>::const_iterator i=T.begin();i!=T.end();i++)
    {
      FieldMatrix S=A2.submatrixRows(*i);

      FieldMatrix S1(Q,S.getHeight()-1,S.getWidth());
      for(int j=0;j<S1.getHeight();j++)
	S1[j]=S[j+1]-S[0];

      FieldMatrix S2=S1*(S1.transposed());
      FieldElement square=S2.reduceAndComputeDeterminant();

      //      ret+=square.squareroot();
    }

  return ret;
}


#include "halfopencone.h"
bool PolyhedralCone::hasFace(PolyhedralCone const &f)const
{
  if(!contains(f))return false;
  IntegerVectorList linealitySpace=this->linealitySpace().dualCone().getEquations();
  IntegerVectorList rays=extremeRays();

  for(IntegerVectorList::const_iterator i=linealitySpace.begin();i!=linealitySpace.end();i++)
    if(!f.contains(*i))return false;

  IntegerVectorList strict;
  for(IntegerVectorList::const_iterator i=rays.begin();i!=rays.end();i++)
    if(!f.contains(*i))
      strict.push_back(*i);

  IntegerVectorList linealitySpace2=f.linealitySpace().dualCone().getEquations();
  IntegerVectorList rays2=f.extremeRays();

  for(IntegerVectorList::const_iterator i=rays2.begin();i!=rays2.end();i++)
    {
      linealitySpace2.push_back(*i);
    }
  IntegerVectorList empty;
  HalfOpenCone C(n,linealitySpace2, empty, strict);
  return !C.isEmpty();
}

PolyhedralCone PolyhedralCone::faceContaining(IntegerVector const &v)const
{
	assert(n==v.size());
	assert(contains(v));
	IntegerVectorList newEquations=equations;
	IntegerVectorList newInequalities;
	for(IntegerVectorList::const_iterator i=inequalities.begin();i!=inequalities.end();i++)
		if(dotLong(*i,v))
			newInequalities.push_back(*i);
		else
			newEquations.push_back(*i);

	PolyhedralCone ret(newInequalities,newEquations,n,(state>=1)?PCP_impliedEquationsKnown:0);
	ret.ensureStateAsMinimum(state);
	return ret;
}

PolyhedralCone PolyhedralCone::projection(int newn)const
{
	assert(newn<=n);
	assert(newn>=0);
	IntegerVectorList rays=extremeRays();
	IntegerVectorList lines=linealitySpace().generatorsOfSpan();
	rays=rowsToIntegerMatrix(rays,n).submatrix(0,0,rays.size(),newn).getRows();
	lines=rowsToIntegerMatrix(lines,n).submatrix(0,0,lines.size(),newn).getRows();
	return givenByRays(rays,lines,newn);
}
